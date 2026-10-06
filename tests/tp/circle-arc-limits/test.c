/* SPDX-License-Identifier: GPL-2.0-only */
#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <string.h>

#include <rtapi.h>
#include <posemath.h>
#include "tp/blendmath.h"
#include "tp/tc.h"
#include "motion/motion.h"

/* TC reads this status during arc-limit calculations. */
static emcmot_status_t status;
emcmot_status_t *emcmotStatus = &status;

void rtapi_print_msg(msg_level_t level, const char *format, ...)
{
    (void)level;
    (void)format;
}

static int failures;

#define CHECK(condition) do { \
    if (!(condition)) { \
        fprintf(stderr, "%s:%d: %s\n", __FILE__, __LINE__, #condition); \
        failures++; \
    } \
} while (0)

static void near(double actual, double expected, const char *name)
{
    if (!isfinite(actual) || !isfinite(expected) ||
        fabs(actual - expected) > 1e-12 * fmax(1.0, fabs(expected))) {
        fprintf(stderr, "%s: actual %.17g, expected %.17g\n",
                name, actual, expected);
        failures++;
    }
}

static TC_STRUCT segment(double radius, double spiral, double angle, double rise)
{
    TC_STRUCT tc = {0};
    tc.motion_type = TC_CIRCULAR;
    tc.term_cond = TC_TERM_COND_STOP;
    tc.maxvel = 1000.0;
    tc.maxaccel = 100.0;
    tc.cycle_time = 0.001;
    tc.target = 100.0;
    tc.acc_ratio_tan = 0.123;
    tc.coords.circle.xyz.radius = radius;
    tc.coords.circle.xyz.spiral = spiral;
    tc.coords.circle.xyz.angle = angle;
    tc.coords.circle.xyz.rTan.x = radius;
    tc.coords.circle.xyz.rPerp.y = radius;
    tc.coords.circle.xyz.normal.z = 1.0;
    tc.coords.circle.xyz.rHelix.z = rise;
    tc.coords.circle.abc.tmag_zero = 1;
    tc.coords.circle.uvw.tmag_zero = 1;
    return tc;
}

static void test_helix_acceleration(void)
{
    status.planner_type = 0;
    TC_STRUCT full = segment(10.0, 0.0, PM_2_PI, 20.0);
    double rho = 10.0 + pow(20.0 / PM_2_PI, 2) / 10.0;
    CHECK(tcUpdateArcLimits(&full) == TP_ERR_OK);
    near(full.maxvel, sqrt(BLEND_ACC_RATIO_NORMAL * 100.0 * rho),
            "full-turn geometric acceleration cap");
    near(full.acc_ratio_tan, BLEND_ACC_RATIO_TANGENTIAL,
            "full-turn tangential budget");

    /* The old cap remains tighter for this short sweep. */
    TC_STRUCT short_arc = segment(10.0, 0.0, 0.5, 1.0);
    CHECK(tcUpdateArcLimits(&short_arc) == TP_ERR_OK);
    near(short_arc.maxvel, sqrt(BLEND_ACC_RATIO_NORMAL * 100.0 * 10.1),
            "short-turn legacy compatibility cap");

    TC_STRUCT reversed_rise = segment(10.0, 0.0, PM_2_PI, -20.0);
    CHECK(tcUpdateArcLimits(&reversed_rise) == TP_ERR_OK);
    near(reversed_rise.maxvel, full.maxvel, "rise sign does not change cap");

    TC_STRUCT slow = segment(10.0, 0.0, PM_2_PI, 20.0);
    slow.maxvel = 0.25;
    CHECK(tcUpdateArcLimits(&slow) == TP_ERR_OK);
    near(slow.maxvel, 0.25, "existing segment velocity cap");
}

static void test_spiral_compatibility(void)
{
    status.planner_type = 0;
    TC_STRUCT spiral = segment(1.0, 19.0, 19.0, 190.0);
    spiral.maxvel = 10000.0;
    double old_radius = pmCircleLegacyMinRadius(&spiral.coords.circle.xyz);
    double expected = sqrt(BLEND_ACC_RATIO_NORMAL * spiral.maxaccel * old_radius);
    CHECK(tcUpdateArcLimits(&spiral) == TP_ERR_OK);
    near(spiral.maxvel, expected, "general spiral keeps legacy motion cap");

    /* A small nonzero spiral still uses the old motion cap. */
    TC_STRUCT tiny = segment(10.0, 1e-10, PM_2_PI, 20.0);
    old_radius = pmCircleLegacyMinRadius(&tiny.coords.circle.xyz);
    CHECK(tcUpdateArcLimits(&tiny) == TP_ERR_OK);
    near(tiny.maxvel, sqrt(BLEND_ACC_RATIO_NORMAL * 100.0 * old_radius),
            "nonzero spiral compatibility");
}

static void test_individual_jerk_caps(void)
{
    status.planner_type = 1;
    /* Choose jerk and period to make each old cap active in turn. */
    const double jerks[] = {400.0, 10.0, 400.0};
    const double periods[] = {1.0, 1.0, 0.001};
    const double expected[] = {20.0, 2.0 / 3.0, 2.0};
    for (unsigned i = 0; i < 3; ++i) {
        TC_STRUCT tc = segment(10.0, 0.0, 0.5, 0.0);
        tc.cycle_time = periods[i];
        status.jerk = jerks[i];
        CHECK(tcUpdateArcLimits(&tc) == TP_ERR_OK);
        near(tc.maxvel, expected[i], "isolated legacy jerk cap");
    }

    /* The helix correction must not relax any old jerk cap. */
    const double sweeps[] = {0.5, PM_2_PI};
    for (unsigned i = 0; i < 2; ++i) {
        TC_STRUCT tc = segment(10.0, 0.0, sweeps[i], 20.0 * sweeps[i] / PM_2_PI);
        double old_radius = pmCircleLegacyMinRadius(&tc.coords.circle.xyz);
        status.jerk = 400.0;
        double old_steady = cbrt(old_radius * old_radius * status.jerk *
                sweeps[i] / (2.0 + sweeps[i]));
        double old_tangent = status.jerk * old_radius /
                (3.0 * BLEND_ACC_RATIO_TANGENTIAL * tc.maxaccel);
        double old_entry = sqrt(status.jerk * old_radius * tc.cycle_time);
        CHECK(tcUpdateArcLimits(&tc) == TP_ERR_OK);
        CHECK(isfinite(tc.maxvel));
        CHECK(tc.maxvel <= old_steady * (1.0 + 1e-12));
        CHECK(tc.maxvel <= old_tangent * (1.0 + 1e-12));
        CHECK(tc.maxvel <= old_entry * (1.0 + 1e-12));
    }
}

static void test_errors_are_atomic(void)
{
    status.planner_type = 0;
    for (unsigned i = 0; i < 7; ++i) {
        TC_STRUCT tc = segment(10.0, 0.0, PM_2_PI, 20.0);
        switch (i) {
            case 0: tc.coords.circle.xyz.angle = 0.0; break;
            case 1: tc.coords.circle.xyz.radius = -1.0; break;
            case 2: tc.coords.circle.xyz.spiral = -10.0; break;
            case 3: tc.coords.circle.xyz.rHelix.z = NAN; break;
            case 4: tc.maxvel = NAN; break;
            case 5: tc.maxaccel = 0.0; break;
            case 6: tc.cycle_time = 0.0; break;
        }
        TC_STRUCT before = tc;
        CHECK(tcUpdateArcLimits(&tc) == TP_ERR_FAIL);
        CHECK(memcmp(&tc, &before, sizeof(tc)) == 0);
    }
    CHECK(tcUpdateArcLimits(NULL) == TP_ERR_FAIL);

    TC_STRUCT invalid = segment(10.0, 0.0, 0.0, 20.0);
    invalid.target = 0.01; /* finalization would first clamp velocity to 10 */
    TC_STRUCT before = invalid;
    CHECK(tcFinalizeLength(&invalid) == TP_ERR_FAIL);
    CHECK(memcmp(&invalid, &before, sizeof(invalid)) == 0);
    CHECK(tcFinalizeLength(NULL) == TP_ERR_FAIL);
}

static void test_finalization_and_nonarcs(void)
{
    status.planner_type = 0;
    TC_STRUCT line = segment(10.0, 0.0, PM_2_PI, 20.0);
    line.motion_type = TC_LINEAR;
    line.target = 0.01;
    TC_STRUCT before = line;
    CHECK(tcUpdateArcLimits(&line) == TP_ERR_NO_ACTION);
    CHECK(memcmp(&line, &before, sizeof(line)) == 0);
    CHECK(tcFinalizeLength(&line) == TP_ERR_OK);
    CHECK(line.finalized);
    near(line.maxvel, 10.0, "length cap for line");
    before = line;
    CHECK(tcFinalizeLength(&line) == TP_ERR_NO_ACTION);
    CHECK(memcmp(&line, &before, sizeof(line)) == 0);

    TC_STRUCT circle = segment(10.0, 0.0, PM_2_PI, 20.0);
    circle.target = 0.001;
    CHECK(tcFinalizeLength(&circle) == TP_ERR_OK);
    CHECK(circle.finalized);
    near(circle.maxvel, 1.0, "length cap for helix");
}

static void test_circle_initialization(void)
{
    EmcPose start = {0}, end = {0};
    PmCartesian center = {0, 0, 0}, normal = {0, 0, 1};
    PmCircle9 circle = {0};
    start.tran.x = end.tran.x = 10.0;
    end.tran.z = 20.0;
    CHECK(pmCircle9Init(&circle, &start, &end, &center, &normal, 0) == TP_ERR_OK);
    near(circle.xyz.angle, PM_2_PI, "initialized helix sweep");
    near(pmCircle9Target(&circle), sqrt(pow(10.0 * PM_2_PI, 2) + 400.0),
            "initialized helix spatial target");

    /* Posemath accepts this axis endpoint; TP rejects it. */
    end.tran.x = 0.0;
    CHECK(pmCircle9Init(&circle, &start, &end, &center, &normal, 0) == TP_ERR_FAIL);

    normal.z = 0.0;
    CHECK(pmCircle9Init(&circle, &start, &end, &center, &normal, 0) == TP_ERR_FAIL);
}

static void test_circle_resize_is_atomic(void)
{
    status.planner_type = 0;
    TC_STRUCT tc = segment(10.0, 0.0, PM_2_PI, 0.0);
    tc.cycle_time = 1.0;
    PmCircle trimmed = tc.coords.circle.xyz;
    CHECK(pmCircleStretch(&trimmed, 0.05, 0) == PM_OK);
    CHECK(tcSetCircleXYZ(&tc, &trimmed) == TP_ERR_OK);
    near(tc.target, 0.5, "trimmed circle target");
    near(tc.maxvel, 0.5, "trimmed circle length cap");
    near(tc.coords.circle.xyz.angle, 0.05, "trimmed circle sweep");
    CHECK(!tc.finalized);

    for (unsigned i = 0; i < 3; ++i) {
        TC_STRUCT candidate = tc;
        PmCircle invalid = trimmed;
        if (i == 0) {
            invalid.spiral = -invalid.radius;
        } else if (i == 1) {
            /* Valid scalar curvature, unrepresentable fitted target. */
            invalid.radius = invalid.rTan.x = invalid.rPerp.y = 1e200;
            invalid.angle = 1e200;
        } else {
            candidate.maxaccel = 0.0;
        }
        TC_STRUCT before = candidate;
        CHECK(tcSetCircleXYZ(&candidate, &invalid) == TP_ERR_FAIL);
        CHECK(memcmp(&candidate, &before, sizeof(candidate)) == 0);
    }
    CHECK(tcSetCircleXYZ(NULL, &trimmed) == TP_ERR_FAIL);
    CHECK(tcSetCircleXYZ(&tc, NULL) == TP_ERR_FAIL);
}

static void test_prepared_blend_budget(void)
{
    for (int overlap = 0; overlap < 2; ++overlap) {
        TC_STRUCT previous = segment(10.0, 0.0, PM_2_PI, 0.0);
        TC_STRUCT next = previous;
        previous.term_cond = TC_TERM_COND_PARABOLIC;
        previous.blend_prev = overlap;
        previous.kink_accel_reduce = 0.7;
        previous.kink_accel_reduce_prev = 0.2;
        next.term_cond = overlap ? TC_TERM_COND_PARABOLIC : TC_TERM_COND_STOP;
        next.blend_prev = 1;
        next.kink_accel_reduce = 0.1;
        next.kink_accel_reduce_prev = 0.7;

        double previous_budget = tcGetArcBlendMaxAccel(&previous, 1);
        double next_budget = tcGetArcBlendMaxAccel(&next, 0);
        PmCircle9 geometry;
        double target, maxvel, ratio;
        CHECK(tcPrepareCircleXYZ(&previous, &previous.coords.circle.xyz,
                    previous_budget, &geometry, &target, &maxvel, &ratio) == TP_ERR_OK);

        /* The prepared budget must match the final blend flags. */
        CHECK(tcRemoveKinkProperties(&previous, &next) == 0);
        CHECK(tcSetTermCond(&previous, &next, TC_TERM_COND_TANGENT) == 0);
        near(previous_budget, tcGetOverallMaxAccel(&previous), "previous blend budget");
        near(next_budget, tcGetOverallMaxAccel(&next), "next blend budget");
        CHECK(tcSetCircleXYZ(&previous, &previous.coords.circle.xyz) == TP_ERR_OK);
        near(maxvel, previous.maxvel, "preflight has no temporary kink slowdown");
        near(ratio, previous.acc_ratio_tan, "preflight tangent budget");
    }
}

int main(void)
{
    test_helix_acceleration();
    test_spiral_compatibility();
    test_individual_jerk_caps();
    test_errors_are_atomic();
    test_finalization_and_nonarcs();
    test_circle_initialization();
    test_circle_resize_is_atomic();
    test_prepared_blend_budget();
    if (failures) {
        fprintf(stderr, "circle arc limits: %d failure(s)\n", failures);
        return 1;
    }
    puts("circle arc limits: OK");
    return 0;
}
