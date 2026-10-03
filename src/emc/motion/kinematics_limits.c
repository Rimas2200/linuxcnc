/* Conservative preflight checks for coordinated joint limits.
 * License: GPL Version 2
 */
#include <limits.h>
#include "kinematics_limits.h"
#include "kinematics_bounds.h"
#include "blendmath.h"
#include "rtapi_math.h"

enum { PATH_MAX_JOINTS = sizeof(unsigned) * CHAR_BIT };

typedef struct {
    double start, end;
    int depth;
} PathInterval;

static int finite_pose(const EmcPose *pose)
{
    return isfinite(pose->tran.x) && isfinite(pose->tran.y)
        && isfinite(pose->tran.z) && isfinite(pose->a) && isfinite(pose->b)
        && isfinite(pose->c) && isfinite(pose->u) && isfinite(pose->v)
        && isfinite(pose->w);
}

static double interpolate(double start, double end, double t)
{
    if (t == 0) return start;
    if (t == 1) return end;
    return (1 - t) * start + t * end;
}

static int path_point(const EmcPose *start, const EmcPose *end,
                      const PmCircle *circle, const SpiralArcLengthFit *fit,
                      double t, EmcPose *point, double *angle)
{
    point->a = interpolate(start->a, end->a, t);
    point->b = interpolate(start->b, end->b, t);
    point->c = interpolate(start->c, end->c, t);
    point->u = interpolate(start->u, end->u, t);
    point->v = interpolate(start->v, end->v, t);
    point->w = interpolate(start->w, end->w, t);
    if (circle) {
        if (t == 0)
            *angle = 0;
        else if (t == 1)
            *angle = circle->angle;
        else if (pmCircleAngleFromParam(circle, fit, t, angle))
            return -1;
        if (!isfinite(*angle) || *angle < 0 || *angle > circle->angle
                || pmCirclePoint(circle, *angle, &point->tran))
            return -1;
    } else {
        point->tran.x = interpolate(start->tran.x, end->tran.x, t);
        point->tran.y = interpolate(start->tran.y, end->tran.y, t);
        point->tran.z = interpolate(start->tran.z, end->tran.z, t);
    }
    return finite_pose(point) ? 0 : -1;
}

static void linear_bounds(KINEMATICS_PATH *path)
{
#define COORDINATE_BOUNDS(member) \
    path->lower.member = fmin(path->start.member, path->end.member); \
    path->upper.member = fmax(path->start.member, path->end.member)
    COORDINATE_BOUNDS(tran.x);
    COORDINATE_BOUNDS(tran.y);
    COORDINATE_BOUNDS(tran.z);
    COORDINATE_BOUNDS(a);
    COORDINATE_BOUNDS(b);
    COORDINATE_BOUNDS(c);
    COORDINATE_BOUNDS(u);
    COORDINATE_BOUNDS(v);
    COORDINATE_BOUNDS(w);
#undef COORDINATE_BOUNDS
}

static KinematicsInterval interval_product(KinematicsInterval a,
                                           KinematicsInterval b)
{
    double aa = a.lower * b.lower, ab = a.lower * b.upper;
    double ba = a.upper * b.lower, bb = a.upper * b.upper;
    return kinematicsIntervalWiden(fmin(fmin(aa, ab), fmin(ba, bb)),
                                    fmax(fmax(aa, ab), fmax(ba, bb)));
}

/* Bound radius, phase and helix separately to avoid root searches at each
 * subdivision. Providers can use the subcircle for tighter bounds. */
static int circle_box(const PmCircle *circle, PmCartesian *lower,
                       PmCartesian *upper)
{
    KinematicsInterval angle = kinematicsInterval(0, circle->angle);
    KinematicsInterval cosine = kinematicsIntervalTrig(angle, 1);
    KinematicsInterval sine = kinematicsIntervalTrig(angle, 0);
    double end_scale = 1 + circle->spiral / circle->radius;
    KinematicsInterval radius = kinematicsIntervalWiden(fmin(1, end_scale),
                                                         fmax(1, end_scale));
    KinematicsInterval component;
    if (!kinematicsIntervalValid(radius))
        return -1;
#define CIRCLE_COORDINATE(member) \
    component = kinematicsIntervalAdd( \
        kinematicsIntervalScale(cosine, circle->rTan.member), \
        kinematicsIntervalScale(sine, circle->rPerp.member)); \
    if (!kinematicsIntervalValid(component)) return -1; \
    component = interval_product(component, radius); \
    component = kinematicsIntervalAdd(component, kinematicsInterval( \
        fmin(0, circle->rHelix.member), fmax(0, circle->rHelix.member))); \
    component = kinematicsIntervalAdd(component, \
        kinematicsInterval(circle->center.member, circle->center.member)); \
    if (!kinematicsIntervalValid(component)) return -1; \
    lower->member = component.lower; \
    upper->member = component.upper
    CIRCLE_COORDINATE(x);
    CIRCLE_COORDINATE(y);
    CIRCLE_COORDINATE(z);
#undef CIRCLE_COORDINATE
    return 0;
}

static int make_path(const EmcPose *start, const EmcPose *end,
                     const PmCircle *circle, const SpiralArcLengthFit *fit,
                     double t0, double t1, KINEMATICS_PATH *path,
                     PmCircle *subcircle)
{
    double a0 = 0, a1 = 0;
    if (path_point(start, end, circle, fit, t0, &path->start, &a0)
            || path_point(start, end, circle, fit, t1, &path->end, &a1))
        return -1;
    linear_bounds(path);
    path->circle = NULL;
    if (circle && t0 != t1) {
        if (pmCircleSubsegment(circle, a0, a1, subcircle)
                || circle_box(subcircle, &path->lower.tran, &path->upper.tran))
            return -1;
        path->circle = subcircle;
    }
    return finite_pose(&path->lower) && finite_pose(&path->upper) ? 0 : -1;
}

/* Pass a copy of the inverse flags and reject providers that change it. */
static int provider_bounds(KINEMATICS_INVERSE_BOUNDS bounds,
                           const KINEMATICS_PATH *path,
                           double *lower, double *upper, int num_joints,
                           unsigned active_mask, KINEMATICS_INVERSE_FLAGS flags)
{
    KINEMATICS_INVERSE_FLAGS copy = flags;
    int i, result;
    for (i = 0; i < num_joints; i++) {
        lower[i] = DBL_MAX;
        upper[i] = -DBL_MAX;
    }
    result = bounds(path, lower, upper, num_joints, &copy);
    if (copy != flags)
        return KINEMATICS_PATH_INVALID;
    if (result == KINEMATICS_BOUNDS_OK) {
        for (i = 0; i < num_joints; i++) {
            if ((active_mask & (1U << i)) && (!isfinite(lower[i])
                    || !isfinite(upper[i]) || lower[i] > upper[i]))
                return KINEMATICS_PATH_INVALID;
        }
    } else if (result != KINEMATICS_BOUNDS_UNSUPPORTED
            && result != KINEMATICS_BOUNDS_UNKNOWN) {
        return KINEMATICS_PATH_INVALID;
    }
    return result;
}

int kinematicsCheckPath(const EmcPose *start, const EmcPose *end,
                       const PmCircle *circle, KINEMATICS_INVERSE_BOUNDS bounds,
                       const double *joint_min, const double *joint_max,
                       int num_joints, unsigned active_mask,
                       const KINEMATICS_INVERSE_FLAGS *iflags,
                       int *failing_joint, int *direction)
{
    PathInterval stack[KINEMATICS_PATH_MAX_DEPTH + 1];
    double lower[PATH_MAX_JOINTS], upper[PATH_MAX_JOINTS];
    SpiralArcLengthFit fit;
    KINEMATICS_INVERSE_FLAGS flags;
    int i, count = 1, boxes = 0;
    const double tolerance = 1e-12;
    if (failing_joint) *failing_joint = -1;
    if (direction) *direction = 0;
    if (!bounds)
        return KINEMATICS_PATH_UNSUPPORTED;
    if (!start || !end || !joint_min || !joint_max || !iflags
            || num_joints <= 0 || num_joints > PATH_MAX_JOINTS
            || !finite_pose(start) || !finite_pose(end)
            || (num_joints < PATH_MAX_JOINTS && (active_mask >> num_joints)))
        return KINEMATICS_PATH_INVALID;
    for (i = 0; i < num_joints; i++) {
        if ((active_mask & (1U << i)) && (!isfinite(joint_min[i])
                || !isfinite(joint_max[i]) || joint_min[i] > joint_max[i]))
            return KINEMATICS_PATH_INVALID;
    }
    if (!active_mask)
        return KINEMATICS_PATH_OK;
    flags = *iflags;
    if (circle) {
        PmCircle validated;
        if (pmCircleSubsegment(circle, 0, circle->angle, &validated)
                || findSpiralArcLengthFit(circle, &fit)
                || !isfinite(fit.b0) || !isfinite(fit.b1)
                || !isfinite(fit.total_planar_length) || fit.total_planar_length <= 0)
            return KINEMATICS_PATH_INVALID;
    }
    stack[0].start = 0;
    stack[0].end = 1;
    stack[0].depth = 0;
    while (count) {
        KINEMATICS_PATH path;
        PmCircle subcircle;
        PathInterval interval = stack[--count];
        double middle = (interval.start + interval.end) / 2;
        int result, contained = 1;
        if (++boxes > KINEMATICS_PATH_MAX_BOXES)
            return KINEMATICS_PATH_UNCERTAIN;
        if (make_path(start, end, circle, &fit, interval.start, interval.end,
                      &path, &subcircle))
            return KINEMATICS_PATH_INVALID;
        result = provider_bounds(bounds, &path, lower, upper, num_joints,
                                 active_mask, flags);
        if (result == KINEMATICS_BOUNDS_UNSUPPORTED)
            return boxes == 1 ? KINEMATICS_PATH_UNSUPPORTED : KINEMATICS_PATH_INVALID;
        if (result == KINEMATICS_PATH_INVALID)
            return result;
        if (result == KINEMATICS_BOUNDS_OK) {
            for (i = 0; i < num_joints; i++) {
                if ((active_mask & (1U << i)) && (lower[i] < joint_min[i] - tolerance
                        || upper[i] > joint_max[i] + tolerance))
                    contained = 0;
            }
            if (contained)
                continue;
        }

        /* The box may be too wide. Check a path point before reporting a
         * limit violation; its whole enclosure must be outside the limit. */
        if (make_path(start, end, circle, &fit, middle, middle, &path, &subcircle))
            return KINEMATICS_PATH_INVALID;
        result = provider_bounds(bounds, &path, lower, upper, num_joints,
                                 active_mask, flags);
        if (result == KINEMATICS_BOUNDS_UNSUPPORTED || result == KINEMATICS_PATH_INVALID)
            return KINEMATICS_PATH_INVALID;
        if (result == KINEMATICS_BOUNDS_OK) {
            for (i = 0; i < num_joints; i++) {
                if (!(active_mask & (1U << i)))
                    continue;
                if (upper[i] < joint_min[i] - tolerance
                        || lower[i] > joint_max[i] + tolerance) {
                    if (failing_joint) *failing_joint = i;
                    if (direction) *direction = upper[i] < joint_min[i] - tolerance ? -1 : 1;
                    return KINEMATICS_PATH_LIMIT;
                }
            }
        }
        if (interval.depth >= KINEMATICS_PATH_MAX_DEPTH)
            return KINEMATICS_PATH_UNCERTAIN;
        stack[count].start = middle;
        stack[count].end = interval.end;
        stack[count++].depth = interval.depth + 1;
        stack[count].start = interval.start;
        stack[count].end = middle;
        stack[count++].depth = interval.depth + 1;
    }
    return KINEMATICS_PATH_OK;
}
