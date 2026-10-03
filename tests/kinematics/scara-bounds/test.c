#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "kinematics_bounds.h"

static void require(int condition, const char *message)
{
    if (!condition) {
        fprintf(stderr, "%s\n", message);
        exit(1);
    }
}

/* Match the point inverse in scarakins.c for comparison with the bounds. */
static void inverse(const EmcPose *world, const double d[6], unsigned long flags,
                    double joint[6])
{
    double angle = world->c * (PM_PI / 180);
    double x = world->tran.x - d[5] * cos(angle);
    double y = world->tran.y - d[5] * sin(angle);
    double cosine = (x * x + y * y - d[1] * d[1] - d[3] * d[3])
        / (2 * d[1] * d[3]);
    double q0, q1;
    if (cosine < -1) cosine = -1;
    if (cosine > 1) cosine = 1;
    q1 = acos(cosine);
    if (flags) q1 = -q1;
    q0 = atan2(y, x) - atan2(d[3] * sin(q1), d[1] + d[3] * cos(q1));
    q0 *= 180 / PM_PI;
    q1 *= 180 / PM_PI;
    joint[0] = q0;
    joint[1] = q1;
    joint[2] = d[0] + d[2] - d[4] - world->tran.z;
    joint[3] = world->c - (q0 + q1);
    joint[4] = world->a;
    joint[5] = world->b;
}

static double random_unit(void)
{
    static unsigned state = 3839;
    state = 1664525U * state + 1013904223U;
    return (state >> 8) / 16777216.0;
}

static double mix(double low, double high, double t)
{
    if (t == 0) return low;
    if (t == 1) return high;
    return (1 - t) * low + t * high;
}

static void check_box(EmcPose low, EmcPose high, const double d[6])
{
    double lower[6], upper[6], point[6];
    unsigned flag;
    int sample, joint;
    for (flag = 0; flag < 2; flag++) {
        require(scaraInverseBounds(&low, &high, d, flag, lower, upper) == 0,
                "valid box rejected");
        for (sample = 0; sample < 128; sample++) {
            EmcPose p = {0};
            double fraction[6];
            for (joint = 0; joint < 6; joint++)
                fraction[joint] = sample < 64 ? !!(sample & (1 << joint)) : random_unit();
            p.tran.x = mix(low.tran.x, high.tran.x, fraction[0]);
            p.tran.y = mix(low.tran.y, high.tran.y, fraction[1]);
            p.tran.z = mix(low.tran.z, high.tran.z, fraction[2]);
            p.a = mix(low.a, high.a, fraction[3]);
            p.b = mix(low.b, high.b, fraction[4]);
            p.c = mix(low.c, high.c, fraction[5]);
            inverse(&p, d, flag, point);
            for (joint = 0; joint < 6; joint++) {
                if (!isfinite(point[joint]) || !isfinite(lower[joint])
                        || !isfinite(upper[joint]) || point[joint] < lower[joint]
                        || point[joint] > upper[joint]) {
                    fprintf(stderr, "joint %d, elbow %u: %.17g outside [%.17g, %.17g]\n",
                            joint, flag, point[joint], lower[joint], upper[joint]);
                    exit(1);
                }
            }
        }
    }
}

static void analytic_cases(void)
{
    const double d[6] = {10, 6, 2, 4, 3, 0};
    double lower[6], upper[6];
    EmcPose p = {{6, 4, 7}, 12, -34, 0, 0, 0, 0};
    require(scaraInverseBounds(&p, &p, d, 0, lower, upper) == 0, "point failed");
    require(lower[0] <= 0 && upper[0] >= 0 && upper[0] - lower[0] < 1e-9,
            "right-angle shoulder bound");
    require(lower[1] <= 90 && upper[1] >= 90 && upper[1] - lower[1] < 1e-9,
            "right-angle elbow bound");
    require(lower[2] <= 2 && upper[2] >= 2, "vertical translation bound");
    require(lower[3] <= -90 && upper[3] >= -90, "wrist rotation bound");
    require(lower[4] == 12 && upper[4] == 12 && lower[5] == -34 && upper[5] == -34,
            "pass-through bounds");
    check_box(p, p, d);
    p.tran.x = -8;
    p.tran.y = 0;
    {
        EmcPose low = p, high = p;
        low.tran.y = -0.01;
        high.tran.y = 0.01;
        require(scaraInverseBounds(&low, &high, d, 0, lower, upper) == 0,
                "branch cut failed");
        require(upper[0] - lower[0] >= 360, "branch cut must include both sides");
        check_box(low, high, d);
    }
    p.tran.x = p.tran.y = 0;
    check_box(p, p, d);
    p.tran.x = 100; /* The production inverse clamps unreachable radii. */
    check_box(p, p, d);
    p.tran.x = 1;
    check_box(p, p, d);
    p.tran.x = 10;
    check_box(p, p, d);
    {
        const double reverse[6] = {10, 4, 2, 6, 3, 0};
        const double equal[6] = {10, 5, 2, 5, 3, 0};
        p.tran.x = 0;
        check_box(p, p, reverse);
        check_box(p, p, equal);
    }
    {
        /* The elbow angle is largest at the midpoint, nearest the shoulder. */
        EmcPose low = {{-6, 8, 0}, 0, 0, 0, 0, 0, 0};
        EmcPose high = {{6, 8, 0}, 0, 0, 0, 0, 0, 0};
        double maximum = acos(0.25) * 180 / PM_PI;
        require(scaraInverseBounds(&low, &high, d, 0, lower, upper) == 0,
                "interior elbow extremum failed");
        require(lower[1] <= 0 && upper[1] >= maximum
                && upper[1] - maximum < 1e-9, "interior elbow extremum missing");
        check_box(low, high, d);
    }
    {
        const double offset[6] = {490, 340, 50, 250, 50, 50};
        EmcPose low = {{240, 150, -10}, -20, -30, -1080, 0, 0, 0};
        EmcPose high = {{250, 160, 10}, 20, 30, 1080, 0, 0, 0};
        check_box(low, high, offset);
        low.c = high.c = 1e100;
        check_box(low, high, offset);
    }
}

static void invalid_cases(void)
{
    EmcPose low = {0}, high = {0};
    double d[6] = {490, 340, 50, 250, 50, 50};
    double lower[6], upper[6];
    require(scaraInverseBounds(NULL, &high, d, 0, lower, upper) < 0, "null accepted");
    low.tran.x = 1;
    require(scaraInverseBounds(&low, &high, d, 0, lower, upper) < 0, "reversed box accepted");
    low.tran.x = NAN;
    require(scaraInverseBounds(&low, &high, d, 0, lower, upper) < 0, "NaN accepted");
    low.tran.x = 0;
    high.w = INFINITY;
    require(scaraInverseBounds(&low, &high, d, 0, lower, upper) < 0, "infinity accepted");
    high.w = 0;
    d[1] = 0;
    require(scaraInverseBounds(&low, &high, d, 0, lower, upper) < 0, "zero arm accepted");
    d[1] = -1;
    require(scaraInverseBounds(&low, &high, d, 0, lower, upper) < 0, "negative arm accepted");
    d[1] = INFINITY;
    require(scaraInverseBounds(&low, &high, d, 0, lower, upper) < 0, "infinite arm accepted");
    d[1] = DBL_MAX;
    require(scaraInverseBounds(&low, &high, d, 0, lower, upper) < 0, "overflow accepted");
    d[1] = DBL_MIN;
    d[3] = DBL_MIN;
    require(scaraInverseBounds(&low, &high, d, 0, lower, upper) < 0, "underflow accepted");
}

int main(void)
{
    int trial;
    analytic_cases();
    invalid_cases();
    for (trial = 0; trial < 2000; trial++) {
        EmcPose low = {0}, high = {0};
        double d[6] = {490, 1 + 500 * random_unit(), 50,
                      1 + 500 * random_unit(), 50, 200 * random_unit() - 100};
        double span = trial % 3 ? 100 : 1e-6;
        low.tran.x = 1600 * random_unit() - 800;
        low.tran.y = 1600 * random_unit() - 800;
        low.tran.z = 200 * random_unit() - 100;
        low.a = 100 * random_unit() - 50;
        low.b = 100 * random_unit() - 50;
        low.c = 2000 * random_unit() - 1000;
        high.tran.x = low.tran.x + span * random_unit();
        high.tran.y = low.tran.y + span * random_unit();
        high.tran.z = low.tran.z + span * random_unit();
        high.a = low.a + span * random_unit();
        high.b = low.b + span * random_unit();
        high.c = low.c + span * random_unit();
        check_box(low, high, d);
    }
    puts("SCARA inverse bounds passed");
    return 0;
}
