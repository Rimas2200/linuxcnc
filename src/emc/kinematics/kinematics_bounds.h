/* Conservative interval helpers for inverse kinematics.
 *
 * License: GPL Version 2
 */
#ifndef KINEMATICS_BOUNDS_H
#define KINEMATICS_BOUNDS_H

#include "emcpos.h"
#include "rtapi_math.h"

typedef struct {
    double lower, upper;
} KinematicsInterval;

static inline KinematicsInterval kinematicsInterval(double lower, double upper)
{
    KinematicsInterval result = {lower, upper};
    return result;
}

/* Allow for rounding in arithmetic and libm without changing the thread's
 * rounding mode. Nonfinite bounds are rejected by the caller. */
static inline KinematicsInterval kinematicsIntervalWiden(double lower,
                                                        double upper)
{
    double scale = fmax(DBL_MIN, fmax(fabs(lower), fabs(upper)));
    double padding = 16 * DBL_EPSILON * scale;
    return kinematicsInterval(lower - padding, upper + padding);
}

static inline int kinematicsIntervalValid(KinematicsInterval value)
{
    return isfinite(value.lower) && isfinite(value.upper)
        && value.lower <= value.upper;
}

static inline KinematicsInterval kinematicsIntervalAdd(KinematicsInterval a,
                                                       KinematicsInterval b)
{
    return kinematicsIntervalWiden(a.lower + b.lower, a.upper + b.upper);
}

static inline KinematicsInterval kinematicsIntervalSubtract(KinematicsInterval a,
                                                            KinematicsInterval b)
{
    return kinematicsIntervalWiden(a.lower - b.upper, a.upper - b.lower);
}

static inline KinematicsInterval kinematicsIntervalScale(KinematicsInterval a,
                                                         double scale)
{
    if (scale >= 0)
        return kinematicsIntervalWiden(a.lower * scale, a.upper * scale);
    return kinematicsIntervalWiden(a.upper * scale, a.lower * scale);
}

static inline KinematicsInterval kinematicsIntervalSquare(KinematicsInterval a)
{
    double left = a.lower * a.lower;
    double right = a.upper * a.upper;
    double lower = a.lower <= 0 && a.upper >= 0 ? 0 : fmin(left, right);
    KinematicsInterval result = kinematicsIntervalWiden(lower, fmax(left, right));
    result.lower = fmax(0, result.lower);
    return result;
}

/* Include endpoint values and interior extrema. For large angles, use
 * [-1,1] to avoid errors in reduction by the double approximation of 2*pi. */
static inline KinematicsInterval kinematicsIntervalTrig(KinematicsInterval a,
                                                        int cosine)
{
    const double period = 2 * PM_PI;
    const double phase = cosine ? 0 : PM_PI / 2;
    double left, right;
    KinematicsInterval result;
    if (!kinematicsIntervalValid(a) || a.upper - a.lower >= period
            || fmax(fabs(a.lower), fabs(a.upper)) > 1048576)
        return kinematicsInterval(-1, 1);

    left = cosine ? cos(a.lower) : sin(a.lower);
    right = cosine ? cos(a.upper) : sin(a.upper);
    result = kinematicsIntervalWiden(fmin(left, right), fmax(left, right));
    a = kinematicsIntervalWiden(a.lower, a.upper);
    if (ceil((a.lower - phase) / period) <= floor((a.upper - phase) / period))
        result.upper = 1;
    if (ceil((a.lower - phase - PM_PI) / period)
            <= floor((a.upper - phase - PM_PI) / period))
        result.lower = -1;
    result.lower = fmax(-1, result.lower);
    result.upper = fmin(1, result.upper);
    return result;
}

/* atan2 has no stationary value inside a quadrant: corner values suffice
 * unless the box touches the negative-X branch cut or contains the origin.
 * Include both sides of that discontinuity, including signed zero. */
static inline KinematicsInterval kinematicsIntervalAtan2(KinematicsInterval y,
                                                         KinematicsInterval x)
{
    double a, b, c, d;
    if (x.lower <= 0 && y.lower <= 0 && y.upper >= 0)
        return kinematicsInterval(-PM_PI, PM_PI);
    a = atan2(y.lower, x.lower);
    b = atan2(y.lower, x.upper);
    c = atan2(y.upper, x.lower);
    d = atan2(y.upper, x.upper);
    return kinematicsIntervalWiden(fmin(fmin(a, b), fmin(c, d)),
                                    fmax(fmax(a, b), fmax(c, d)));
}

/* Bound scarakins' inverse, including its elbow-cosine clamp for unreachable
 * radii. Subdivision tightens the bounds; boxes crossing an angle branch cut
 * can still give wide joint ranges. This does not check reachability.
 *
 * dimensions[] holds D1..D6. Writes joints 0..5 on success and leaves the
 * output unchanged on error. */
static inline int scaraInverseBounds(const EmcPose *lower, const EmcPose *upper,
                                    const double dimensions[6],
                                    unsigned long iflags,
                                    double joint_lower[], double joint_upper[])
{
    KinematicsInterval pose[9], xt, yt, radius2, cosine, elbow, shoulder;
    KinematicsInterval tool_angle, link_x, link_y, result[6], inner2, outer2;
    double denominator;
    int i;
    if (!lower || !upper || !dimensions || !joint_lower || !joint_upper)
        return -1;

    pose[0] = kinematicsInterval(lower->tran.x, upper->tran.x);
    pose[1] = kinematicsInterval(lower->tran.y, upper->tran.y);
    pose[2] = kinematicsInterval(lower->tran.z, upper->tran.z);
    pose[3] = kinematicsInterval(lower->a, upper->a);
    pose[4] = kinematicsInterval(lower->b, upper->b);
    pose[5] = kinematicsInterval(lower->c, upper->c);
    pose[6] = kinematicsInterval(lower->u, upper->u);
    pose[7] = kinematicsInterval(lower->v, upper->v);
    pose[8] = kinematicsInterval(lower->w, upper->w);
    for (i = 0; i < 9; i++) {
        if (!kinematicsIntervalValid(pose[i]))
            return -1;
    }
    for (i = 0; i < 6; i++) {
        if (!isfinite(dimensions[i]))
            return -1;
    }
    if (dimensions[1] <= 0 || dimensions[3] <= 0)
        return -1;

    tool_angle = kinematicsIntervalScale(pose[5], PM_PI / 180);
    xt = kinematicsIntervalSubtract(pose[0],
        kinematicsIntervalScale(kinematicsIntervalTrig(tool_angle, 1), dimensions[5]));
    yt = kinematicsIntervalSubtract(pose[1],
        kinematicsIntervalScale(kinematicsIntervalTrig(tool_angle, 0), dimensions[5]));
    radius2 = kinematicsIntervalAdd(kinematicsIntervalSquare(xt),
                                     kinematicsIntervalSquare(yt));
    /* Use the same rounded constants as the point inverse. */
    inner2 = kinematicsInterval(dimensions[1] * dimensions[1],
                                 dimensions[1] * dimensions[1]);
    outer2 = kinematicsInterval(dimensions[3] * dimensions[3],
                                 dimensions[3] * dimensions[3]);
    denominator = 2 * dimensions[1] * dimensions[3];
    if (!kinematicsIntervalValid(radius2) || !kinematicsIntervalValid(inner2)
            || !kinematicsIntervalValid(outer2) || !isfinite(denominator)
            || denominator <= 0 || !isfinite(1 / denominator))
        return -1;
    cosine = kinematicsIntervalScale(kinematicsIntervalSubtract(
        kinematicsIntervalSubtract(radius2, inner2), outer2), 1 / denominator);
    if (!kinematicsIntervalValid(cosine))
        return -1;
    cosine.lower = fmax(-1, fmin(1, cosine.lower));
    cosine.upper = fmax(-1, fmin(1, cosine.upper));
    elbow = kinematicsIntervalWiden(acos(cosine.upper), acos(cosine.lower));
    if (iflags)
        elbow = kinematicsInterval(-elbow.upper, -elbow.lower);

    link_x = kinematicsIntervalAdd(kinematicsInterval(dimensions[1], dimensions[1]),
        kinematicsIntervalScale(kinematicsIntervalTrig(elbow, 1), dimensions[3]));
    link_y = kinematicsIntervalScale(kinematicsIntervalTrig(elbow, 0), dimensions[3]);
    shoulder = kinematicsIntervalSubtract(kinematicsIntervalAtan2(yt, xt),
                                            kinematicsIntervalAtan2(link_y, link_x));
    result[0] = kinematicsIntervalScale(shoulder, 180 / PM_PI);
    result[1] = kinematicsIntervalScale(elbow, 180 / PM_PI);
    /* Preserve the point inverse's order of operations. */
    denominator = dimensions[0] + dimensions[2] - dimensions[4];
    result[2] = kinematicsIntervalSubtract(kinematicsInterval(denominator, denominator), pose[2]);
    result[3] = kinematicsIntervalSubtract(pose[5],
        kinematicsIntervalAdd(result[0], result[1]));
    result[4] = pose[3];
    result[5] = pose[4];
    for (i = 0; i < 6; i++) {
        if (!kinematicsIntervalValid(result[i]))
            return -1;
    }
    for (i = 0; i < 6; i++) {
        joint_lower[i] = result[i].lower;
        joint_upper[i] = result[i].upper;
    }
    return 0;
}

#endif
