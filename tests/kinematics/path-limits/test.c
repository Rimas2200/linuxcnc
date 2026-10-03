#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include "kinematics_limits.h"
#include "kinematics_bounds.h"
#include "blendmath.h"

/* Standalone replacement for the progress helper's RTAPI diagnostics. */
void rtapi_print_msg(msg_level_t level, const char *format, ...)
{
    va_list args;
    (void)level;
    va_start(args, format);
    vfprintf(stderr, format, args);
    va_end(args);
}

static unsigned calls;
static int behavior;
static PmCircle progress_circle;
static SpiralArcLengthFit progress_fit;

static void require(int condition, const char *message)
{
    if (!condition) {
        fprintf(stderr, "%s\n", message);
        exit(1);
    }
}

static void check_budget(void)
{
    require(calls <= 2 * KINEMATICS_PATH_MAX_BOXES, "callback budget exceeded");
}

static double bump(double x)
{
    double offset = (x - 0.3125) / 0.005;
    return 2 * exp(-offset * offset);
}

static int fake_bounds(const KINEMATICS_PATH *path, double *lower, double *upper,
                       int joints, const KINEMATICS_INVERSE_FLAGS *flags)
{
    int i;
    calls++;
    require(*flags == 42, "inverse flags changed between calls");
    for (i = 0; i < joints; i++)
        lower[i] = upper[i] = 0;
    switch (behavior) {
    case 0:
        lower[0] = path->lower.tran.x;
        upper[0] = path->upper.tran.x;
        return KINEMATICS_BOUNDS_OK;
    case 1: {
        double a = bump(path->lower.tran.x), b = bump(path->upper.tran.x);
        lower[0] = fmin(a, b);
        upper[0] = path->lower.tran.x <= 0.3125 && path->upper.tran.x >= 0.3125
            ? 2 : fmax(a, b);
        return KINEMATICS_BOUNDS_OK;
    }
    case 2:
        return KINEMATICS_BOUNDS_UNSUPPORTED;
    case 3:
        return calls == 1 ? KINEMATICS_BOUNDS_UNKNOWN : KINEMATICS_BOUNDS_UNSUPPORTED;
    case 4:
        return KINEMATICS_BOUNDS_UNKNOWN;
    case 5:
        lower[0] = NAN;
        return KINEMATICS_BOUNDS_OK;
    case 6:
        lower[0] = 1;
        upper[0] = -1;
        return KINEMATICS_BOUNDS_OK;
    case 7:
        lower[1] = upper[1] = NAN;
        return KINEMATICS_BOUNDS_OK;
    case 8:
        *(KINEMATICS_INVERSE_FLAGS *)flags = 0;
        return KINEMATICS_BOUNDS_OK;
    case 9:
        return path->upper.tran.x - path->lower.tran.x > 1.0 / 64
            ? KINEMATICS_BOUNDS_UNKNOWN : KINEMATICS_BOUNDS_OK;
    case 10:
        if (calls == 1) return KINEMATICS_BOUNDS_UNKNOWN;
        return calls == 2 ? KINEMATICS_BOUNDS_OK : KINEMATICS_BOUNDS_UNSUPPORTED;
    case 11:
        return 7;
    case 12:
        upper[0] = INFINITY;
        return KINEMATICS_BOUNDS_OK;
    }
    return KINEMATICS_BOUNDS_UNKNOWN;
}

static int scara_bounds(const KINEMATICS_PATH *path, double *lower, double *upper,
                        int joints, const KINEMATICS_INVERSE_FLAGS *flags)
{
    const double dimensions[6] = {10, 6, 2, 4, 3, 0};
    calls++;
    require(joints == 6, "SCARA joint count");
    return scaraInverseBounds(&path->lower, &path->upper, dimensions, *flags,
                              lower, upper);
}

/* Check XYZ bounds during subdivision. ABCUVW use planner progress, which
 * is not proportional to the angle of a spiral. */
static int progress_bounds(const KINEMATICS_PATH *path, double *lower, double *upper,
                           int joints, const KINEMATICS_INVERSE_FLAGS *flags)
{
    int i;
    double t0 = path->start.c, t1 = path->end.c;
    (void)flags;
    require(joints == 1, "progress joint count");
    calls++;
    for (i = 0; i <= 8; i++) {
        double t = (1 - i / 8.0) * t0 + i / 8.0 * t1;
        double angle;
        PmCartesian point;
        require(pmCircleAngleFromParam(&progress_circle, &progress_fit, t, &angle) == 0,
                "reference angle failed");
        require(pmCirclePoint(&progress_circle, angle, &point) == 0,
                "reference point failed");
        require(point.x >= path->lower.tran.x - 1e-10
                && point.x <= path->upper.tran.x + 1e-10
                && point.y >= path->lower.tran.y - 1e-10
                && point.y <= path->upper.tran.y + 1e-10
                && point.z >= path->lower.tran.z - 1e-10
                && point.z <= path->upper.tran.z + 1e-10,
                "spiral subpath does not contain actual path");
        if (i == 0)
            require(fabs(point.x - path->start.tran.x) < 1e-10
                    && fabs(point.y - path->start.tran.y) < 1e-10
                    && fabs(point.z - path->start.tran.z) < 1e-10,
                    "spiral start has wrong progress");
        if (i == 8)
            require(fabs(point.x - path->end.tran.x) < 1e-10
                    && fabs(point.y - path->end.tran.y) < 1e-10
                    && fabs(point.z - path->end.tran.z) < 1e-10,
                    "spiral end has wrong progress");
    }
    require((t0 == t1) == (path->circle == NULL), "point/circle geometry tag");
    lower[0] = upper[0] = 0;
    return t1 - t0 > 0.04 ? KINEMATICS_BOUNDS_UNKNOWN : KINEMATICS_BOUNDS_OK;
}

static void generic_cases(void)
{
    EmcPose start = {0}, end = {{1, 0, 0}, 0, 0, 0, 0, 0, 0};
    double low[2] = {-0.1, -1}, high[2] = {1, 1};
    KINEMATICS_INVERSE_FLAGS flags = 42;
    int joint, direction, result;
    const int expected[] = {
        KINEMATICS_PATH_OK, KINEMATICS_PATH_LIMIT, KINEMATICS_PATH_UNSUPPORTED,
        KINEMATICS_PATH_INVALID, KINEMATICS_PATH_UNCERTAIN,
        KINEMATICS_PATH_INVALID, KINEMATICS_PATH_INVALID, KINEMATICS_PATH_OK,
        KINEMATICS_PATH_INVALID, KINEMATICS_PATH_UNCERTAIN,
        KINEMATICS_PATH_INVALID, KINEMATICS_PATH_INVALID, KINEMATICS_PATH_INVALID
    };
    for (behavior = 0; behavior < (int)(sizeof(expected) / sizeof(expected[0])); behavior++) {
        calls = 0;
        result = kinematicsCheckPath(&start, &end, NULL, fake_bounds, low, high,
                                     2, 1, &flags, &joint, &direction);
        if (result != expected[behavior]) {
            fprintf(stderr, "case %d: result %d expected %d\n",
                    behavior, result, expected[behavior]);
            exit(1);
        }
        require(flags == 42, "provider altered caller flags");
        check_budget();
        if (behavior == 1)
            require(joint == 0 && direction == 1, "narrow interior violation missing");
        if (behavior == 9)
            require(calls >= KINEMATICS_PATH_MAX_BOXES, "box budget not exercised");
    }
    behavior = 0;
    require(kinematicsCheckPath(&start, &end, NULL, NULL, low, high,
                2, 1, &flags, &joint, &direction) == KINEMATICS_PATH_UNSUPPORTED,
            "missing provider accepted");
    require(kinematicsCheckPath(&start, &end, NULL, fake_bounds, low, high,
                2, 4, &flags, &joint, &direction) == KINEMATICS_PATH_INVALID,
            "invalid joint mask accepted");
    require(kinematicsCheckPath(&start, &end, NULL, fake_bounds, low, high,
                0, 0, &flags, &joint, &direction) == KINEMATICS_PATH_INVALID,
            "invalid joint count accepted");
    start.b = NAN;
    require(kinematicsCheckPath(&start, &end, NULL, fake_bounds, low, high,
                2, 1, &flags, &joint, &direction) == KINEMATICS_PATH_INVALID,
            "nonfinite pose accepted");
    start.b = 0;
    high[0] = -1;
    require(kinematicsCheckPath(&start, &end, NULL, fake_bounds, low, high,
                2, 1, &flags, &joint, &direction) == KINEMATICS_PATH_INVALID,
            "reversed joint limits accepted");
}

static void scara_cases(void)
{
    EmcPose start = {{-6, 8, 0}, 0, 0, 0, 0, 0, 0};
    EmcPose end = {{6, 8, 0}, 0, 0, 0, 0, 0, 0};
    double low[6] = {-360, -1, -100, -360, -100, -100};
    double high[6] = {360, 60, 100, 360, 100, 100};
    KINEMATICS_INVERSE_FLAGS flags = 0;
    int joint, direction;
    calls = 0;
    require(kinematicsCheckPath(&start, &end, NULL, scara_bounds, low, high,
                6, 63, &flags, &joint, &direction) == KINEMATICS_PATH_LIMIT,
            "SCARA elbow interior limit missing");
    require(joint == 1 && direction == 1, "wrong SCARA joint/direction");
    high[1] = 80;
    require(kinematicsCheckPath(&start, &end, NULL, scara_bounds, low, high,
                6, 63, &flags, &joint, &direction) == KINEMATICS_PATH_OK,
            "safe SCARA line rejected");
    flags = 1;
    low[1] = -60;
    high[1] = 1;
    require(kinematicsCheckPath(&start, &end, NULL, scara_bounds, low, high,
                6, 63, &flags, &joint, &direction) == KINEMATICS_PATH_LIMIT,
            "negative SCARA elbow interior limit missing");
    require(joint == 1 && direction == -1, "wrong negative SCARA joint/direction");
    check_budget();
}

static void progress_cases(void)
{
    EmcPose start = {{1, 0, 0}, 0, 0, 0, 0, 0, 0};
    EmcPose end = {{2, 0, 1}, 0, 0, 1, 0, 0, 0};
    PmCartesian center = {0, 0, 0}, normal = {0, 0, 1};
    KINEMATICS_INVERSE_FLAGS flags = 42;
    double low = -1, high = 1;
    int joint, direction, reverse;
    for (reverse = 0; reverse < 2; reverse++) {
        require(pmCircleInit(&progress_circle, &start.tran, &end.tran,
                    &center, &normal, 1) == 0, "spiral init failed");
        require(findSpiralArcLengthFit(&progress_circle, &progress_fit) == 0,
                "spiral fit failed");
        calls = 0;
        require(kinematicsCheckPath(&start, &end, &progress_circle, progress_bounds,
                    &low, &high, 1, 1, &flags, &joint, &direction) == KINEMATICS_PATH_OK,
                "spiral path failed");
        require(calls > 2, "spiral subdivision not exercised");
        check_budget();
        start.tran.x = 2;
        end.tran.x = 1;
    }
    {
        /* Cover tilted spirals, both directions, and multiple turns. */
        const int turns[] = {3, -4, 37, -37};
        const PmCartesian normals[] = {{1, 2, 3}, {-2, 1, 3},
                                      {0, 1, 0}, {1, 0, 1}};
        int trial;
        for (trial = 0; trial < 4; trial++) {
            start.tran.x = 2;
            start.tran.y = -1;
            start.tran.z = 0.5;
            end.tran.x = 1;
            end.tran.y = 3;
            end.tran.z = -2;
            normal = normals[trial];
            require(pmCircleInit(&progress_circle, &start.tran, &end.tran,
                        &center, &normal, turns[trial]) == 0, "tilted spiral init failed");
            require(findSpiralArcLengthFit(&progress_circle, &progress_fit) == 0,
                    "tilted spiral fit failed");
            calls = 0;
            require(kinematicsCheckPath(&start, &end, &progress_circle, progress_bounds,
                        &low, &high, 1, 1, &flags, &joint, &direction) == KINEMATICS_PATH_OK,
                    "tilted spiral enclosure failed");
            check_budget();
        }
    }
    {
        PmCartesian radial = {4, 0, 0}, finish = {0, 4, 3};
        normal.x = normal.y = 0;
        normal.z = 1;
        start.tran = radial;
        end.tran = finish;
        require(pmCircleInit(&progress_circle, &radial, &finish,
                    &center, &normal, 1) == 0, "constant-radius helix init failed");
        require(findSpiralArcLengthFit(&progress_circle, &progress_fit) == 0,
                "constant-radius helix fit failed");
        calls = 0;
        require(kinematicsCheckPath(&start, &end, &progress_circle, progress_bounds,
                    &low, &high, 1, 1, &flags, &joint, &direction) == KINEMATICS_PATH_OK,
                "constant-radius helix enclosure failed");
        check_budget();
    }
    progress_circle.angle = NAN;
    require(kinematicsCheckPath(&start, &end, &progress_circle, progress_bounds,
                &low, &high, 1, 1, &flags, &joint, &direction) == KINEMATICS_PATH_INVALID,
            "invalid circle accepted");
}

int main(void)
{
    generic_cases();
    scara_cases();
    progress_cases();
    puts("kinematics path limits passed");
    return 0;
}
