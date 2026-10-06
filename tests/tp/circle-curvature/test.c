/* SPDX-License-Identifier: GPL-2.0-only */
#include <float.h>
#include <math.h>
#include <stdint.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>

#include "circle_curvature.h"
#include "blendmath.h"

#define PI 3.141592653589793238462643383279502884L
#define SEED UINT32_C(20261005)
#define TOL 1e-12L

static uint32_t random_state = SEED;
static double input[4];
static unsigned int case_number;

/* Print unexpected diagnostics from the production functions. */
void rtapi_print_msg(msg_level_t level, const char *fmt, ...)
{
    va_list ap;
    (void)level;
    va_start(ap, fmt);
    vfprintf(stderr, fmt, ap);
    va_end(ap);
}

static void fail(const char *label, long double actual, long double expected)
{
    fprintf(stderr,
            "%s: actual=%.21Lg expected=%.21Lg; case=%u seed=%u "
            "r0=%.17g delta=%.17g sweep=%.17g H=%.17g\n",
            label, actual, expected, case_number, SEED,
            input[0], input[1], input[2], input[3]);
    exit(EXIT_FAILURE);
}

static void require(int condition, const char *label)
{
    if (!condition)
        fail(label, 0, 1);
}

static void close_to(const char *label, long double actual,
                     long double expected, long double tolerance)
{
    if (!isfinite(actual) || !isfinite(expected) ||
            fabsl(actual - expected) > tolerance * fabsl(expected))
        fail(label, actual, expected);
}

static void context(double r0, double delta, double sweep, double rise)
{
    input[0] = r0;
    input[1] = delta;
    input[2] = sweep;
    input[3] = rise;
    ++case_number;
}

static double kernel(double r0, double delta, double sweep, double rise)
{
    double result = -1234.5;
    context(r0, delta, sweep, rise);
    require(tpCircleMaxCurvature(r0, delta, sweep, rise, &result) == 0,
            "kernel success");
    require(isfinite(result) && result > 0, "positive finite curvature");
    return result;
}

/* Cartesian derivative oracle, independent of the production formula. */
static long double vector_curvature(long double radius, long double k,
                                    long double h, long double theta)
{
    long double c = cosl(theta), s = sinl(theta);
    long double d1[3] = { k*c-radius*s, k*s+radius*c, h };
    long double d2[3] = { -radius*c-2*k*s, -radius*s+2*k*c, 0 };
    long double cross[3] = {
        d1[1]*d2[2]-d1[2]*d2[1],
        d1[2]*d2[0]-d1[0]*d2[2],
        d1[0]*d2[1]-d1[1]*d2[0]
    };
    long double speed = hypotl(hypotl(d1[0], d1[1]), d1[2]);
    return hypotl(hypotl(cross[0], cross[1]), cross[2]) /
           speed / speed / speed;
}

/* Unimodality permits a numerical search without the production root formula. */
static long double reference_max(double r0, double delta,
                                 double sweep, double rise)
{
    long double start = r0, end = start + (long double)delta;
    long double lo = fminl(start, end), hi = fmaxl(start, end);
    long double k = (long double)delta / sweep;
    long double h = (long double)rise / sweep;
    long double maximum = fmaxl(vector_curvature(lo, k, h, 0.37L),
                                vector_curvature(hi, k, h, 0.37L));
    int i;
    for (i = 0; i < 180 && hi > lo; ++i) {
        long double left = lo + (hi-lo)/3;
        long double right = hi - (hi-lo)/3;
        if (vector_curvature(left, k, h, 0.37L) <
                vector_curvature(right, k, h, 0.37L))
            lo = left;
        else
            hi = right;
    }
    return fmaxl(maximum, vector_curvature((lo+hi)/2, k, h, 0.37L));
}

static void reference_case(double r0, double delta, double sweep, double rise)
{
    double actual = kernel(r0, delta, sweep, rise);
    long double expected = reference_max(r0, delta, sweep, rise);
    close_to("vector maximum", actual, expected, TOL);
}

static void test_helices(void)
{
    const int splits[] = { 1, 2, 4, 8, 100 };
    const int turns[] = { 2, 10, 1000 };
    const double sweep = (double)(2*PI);
    const double h = 20/sweep;
    const long double expected = 10/(100+(long double)h*h);
    size_t i;
    close_to("H01 circle", kernel(10, 0, sweep, 0), 0.1L, TOL);
    for (i = 0; i < sizeof(splits)/sizeof(splits[0]); ++i) {
        int count = splits[i];
        close_to("H02/H03 helix segments",
                 kernel(10, 0, sweep/count, 20.0/count), expected, TOL);
    }
    close_to("H04 partial turn", kernel(10, 0, sweep/4, 5), expected, TOL);
    for (i = 0; i < sizeof(turns)/sizeof(turns[0]); ++i)
        close_to("H05 multiple turns",
                 kernel(10, 0, sweep*turns[i], 20.0*turns[i]), expected, TOL);
    close_to("H06 negative rise", kernel(10, 0, sweep, -20), expected, TOL);
    close_to("H09 subradian helix", kernel(10, 0, 0.1, h*0.1), expected, TOL);
    close_to("H10 vanishing rise", kernel(10, 0, 1, 1e-100), 0.1L, TOL);
    close_to("H10 large pitch", kernel(10, 0, 1, 1e100), 1e-199L, TOL);
    close_to("N03 tiny positive sweep", kernel(10, 0, 1e-200, h*1e-200),
             expected, TOL);
    close_to("stable norm with extreme length ratio", kernel(1e-300, 0, 1, 1e-100),
             1e-100L, TOL);
}

static void test_spirals(void)
{
    const double h_threshold = sqrt(8.0);
    const double thresholds[] = { nextafter(h_threshold, 0), h_threshold,
                                 nextafter(h_threshold, INFINITY) };
    double whole, left, right;
    long double root = sqrtl(-3+sqrtl(9301)); /* k=1, h=10 */
    size_t i;
    close_to("S01/S10 planar outward", kernel(1, 19, 19, 0),
             3/(2*sqrtl(2)), TOL);
    close_to("S02 planar inward", kernel(20, -19, 19, 0),
             3/(2*sqrtl(2)), TOL);
    reference_case(1, 19, 19, 19*2);
    for (i = 0; i < sizeof(thresholds)/sizeof(thresholds[0]); ++i) {
        /* Exercise a root near zero. */
        reference_case(1e-9, 1, 1, thresholds[i]);
    }
    whole = kernel(1, 19, 19, 190);
    close_to("S06 internal maximum", whole,
             vector_curvature(root, 1, 10, 0.83L), TOL);
    require(whole > 1.5*vector_curvature(1, 1, 10, 0),
            "S06 greater than starting curvature");
    require(whole > vector_curvature(20, 1, 10, 0),
            "S06 greater than ending curvature");
    reference_case(1, 2, 2, 20);     /* root above interval */
    reference_case(12, 8, 8, 80);   /* root below interval */
    reference_case(1, (double)root-1, (double)root-1, 10*((double)root-1));
    reference_case((double)root, 20-(double)root,
                   20-(double)root, 10*(20-(double)root));
    left = kernel(1, 19*0.375, 19*0.375, 190*0.375);
    right = kernel(1+19*0.375, 19*0.625, 19*0.625, 190*0.625);
    close_to("S08 nonuniform partition", fmax(left, right), whole, TOL);
    close_to("S09 almost constant outward", kernel(10, 1e-14, 2, 8),
             10.0L/116, TOL);
    close_to("S09 almost constant inward", kernel(10, -1e-14, 2, 8),
             10.0L/116, TOL);
    /* radius + delta rounds to radius, but delta/sweep is 1000. */
    reference_case(1, 1e-17, 1e-20, 0);
    reference_case(1, -1e-17, 1e-20, 0);
}

static double random_unit(void)
{
    random_state = UINT32_C(1664525)*random_state + UINT32_C(1013904223);
    return (random_state+0.5)/4294967296.0;
}

static void test_random_and_scale(void)
{
    const double scales[] = { 1e-120, 1e-12, 1e-6, 1, 1e6, 1e12, 1e120 };
    const long double expected = reference_max(1, 19, 19, 190);
    size_t i;
    unsigned int j;
    for (i = 0; i < sizeof(scales)/sizeof(scales[0]); ++i) {
        double scale = scales[i];
        close_to("N01/N02 scale", kernel(scale, 19*scale, 19, 190*scale),
                 expected/scale, TOL);
    }
    for (j = 0; j < 2000; ++j) {
        double scale = pow(10.0, -12+24*random_unit());
        double r0 = scale*(0.01+20*random_unit());
        double r1 = scale*(0.01+20*random_unit());
        double delta = r1-r0;
        double sweep = pow(10.0, -3+6*random_unit());
        double rise = scale*sweep*(-30+60*random_unit());
        double whole = kernel(r0, delta, sweep, rise);
        double reversed = kernel(r0+delta, -delta, sweep, -rise);
        double fraction = 0.1+0.8*random_unit();
        double left = kernel(r0, delta*fraction, sweep*fraction, rise*fraction);
        double right = kernel(r0+delta*fraction, delta*(1-fraction),
                              sweep*(1-fraction), rise*(1-fraction));
        long double ref = reference_max(r0, delta, sweep, rise);
        int sample;
        context(r0, delta, sweep, rise);
        close_to("random reference maximum", whole, ref, TOL);
        close_to("random reversal", reversed, whole, TOL);
        close_to("random partition", fmax(left, right), whole, TOL);
        for (sample = 0; sample < 65; ++sample) {
            long double t = sample/64.0L;
            long double local = vector_curvature(r0+t*delta,
                    (long double)delta/sweep, (long double)rise/sweep, t*sweep);
            if (!isfinite(local) || local > whole*(1+TOL))
                fail("random sampled curvature exceeds maximum", local, whole);
        }
    }
}

static void invalid_kernel(double r0, double delta, double sweep, double rise)
{
    double sentinel = -1234.5;
    context(r0, delta, sweep, rise);
    require(tpCircleMaxCurvature(r0, delta, sweep, rise, &sentinel) != 0,
            "N04/N05/N06 invalid input rejected");
    require(sentinel == -1234.5, "invalid input preserves output");
}

static void test_invalid_kernel(void)
{
    const double bad[] = { NAN, INFINITY, -INFINITY };
    size_t i;
    int field;
    for (i = 0; i < sizeof(bad)/sizeof(bad[0]); ++i) {
        for (field = 0; field < 4; ++field) {
            double values[] = { 1, 2, 3, 4 };
            values[field] = bad[i];
            invalid_kernel(values[0], values[1], values[2], values[3]);
        }
    }
    invalid_kernel(0, 1, 1, 0);
    invalid_kernel(-1, 2, 1, 0);
    invalid_kernel(1, -1, 1, 0);
    invalid_kernel(1, -2, 1, 0);
    invalid_kernel(1, 1, 0, 0);
    invalid_kernel(1, 1, -1, 0);
    invalid_kernel(DBL_MAX, DBL_MAX, 1, 0);
    invalid_kernel(1, 1, DBL_MIN, DBL_MAX);
    invalid_kernel(1, DBL_MAX, DBL_MIN, 0);
    invalid_kernel(DBL_MIN/16, 0, 1, 0); /* curvature is not representable */
    invalid_kernel(DBL_MIN, 0, 1, DBL_MAX); /* curvature rounds to zero */
    invalid_kernel(1e-120, 0, 1, 1e120); /* underflow within finite inputs */
    require(tpCircleMaxCurvature(1, 0, 1, 0, NULL) != 0,
            "N07 null kernel output rejected");
}

static PmCartesian point(PmCartesian center, PmCartesian u, PmCartesian v,
                         PmCartesian normal, double radius, double theta,
                         double rise)
{
    PmCartesian result;
    result.x = center.x+radius*(u.x*cos(theta)+v.x*sin(theta))+rise*normal.x;
    result.y = center.y+radius*(u.y*cos(theta)+v.y*sin(theta))+rise*normal.y;
    result.z = center.z+radius*(u.z*cos(theta)+v.z*sin(theta))+rise*normal.z;
    return result;
}

static double adapter(const PmCircle *circle)
{
    double result = -1234.5;
    context(circle->radius, circle->spiral, circle->angle,
            hypot(hypot(circle->rHelix.x, circle->rHelix.y), circle->rHelix.z));
    require(pmCircleMaxCurvature(circle, &result) == 0, "adapter success");
    require(isfinite(result) && result > 0, "adapter finite curvature");
    close_to("radius compatibility wrapper", pmCircleEffectiveMinRadius(circle),
             1/(long double)result, TOL);
    return result;
}

static void close_point(const char *label, PmCartesian actual, PmCartesian expected)
{
    double distance = hypot(hypot(actual.x-expected.x, actual.y-expected.y),
                            actual.z-expected.z);
    double scale = 1+hypot(hypot(expected.x, expected.y), expected.z);
    if (!isfinite(distance) || distance > 2e-12*scale)
        fail(label, distance, 0);
}

static void test_posemath(void)
{
    const double a = 1/sqrt(2.0), b = 1/sqrt(6.0), c = 1/sqrt(3.0);
    const PmCartesian bases[][3] = {
        { {1,0,0}, {0,1,0}, {0,0,1} },
        { {0,0,1}, {1,0,0}, {0,1,0} },
        { {0,1,0}, {0,0,1}, {1,0,0} },
        { {a,a,0}, {-b,b,2*b}, {c,-c,c} }
    };
    const PmCartesian center = { 5, -7, 11 };
    size_t basis;
    int spiral, sign;
    for (basis = 0; basis < sizeof(bases)/sizeof(bases[0]); ++basis) {
        PmCartesian u = bases[basis][0], v = bases[basis][1], n = bases[basis][2];
        for (spiral = 0; spiral <= 1; ++spiral) {
            for (sign = -1; sign <= 1; sign += 2) {
                const double theta = 4.2, r0 = 10, delta = spiral*3;
                double rise = sign*6;
                PmCartesian start = point(center, u, v, n, r0, 0, 0);
                PmCartesian end = point(center, u, v, n, r0+delta, theta, rise);
                PmCircle circle, reverse;
                double actual;
                int sample;
                require(pmCircleInit(&circle, &start, &end, &center, &n, 0) == 0,
                        "pmCircleInit forward");
                actual = adapter(&circle);
                close_to("H07/H08 transformed circle", actual,
                         reference_max(r0, delta, theta, rise), TOL);
                require(pmCircleInit(&reverse, &end, &start, &center, &n, -1) == 0,
                        "pmCircleInit reverse");
                close_to("H07 reverse constructor", adapter(&reverse), actual, TOL);
                for (sample = 0; sample <= 8; ++sample) {
                    double t = sample/8.0;
                    PmCartesian obtained, reversed;
                    PmCartesian expected = point(center, u, v, n,
                            r0+t*delta, t*theta, t*rise);
                    require(pmCirclePoint(&circle, t*circle.angle, &obtained) == 0,
                            "pmCirclePoint forward");
                    require(pmCirclePoint(&reverse, (1-t)*reverse.angle,
                                          &reversed) == 0, "pmCirclePoint reverse");
                    close_point("pmCirclePoint geometric model", obtained, expected);
                    close_point("pmCirclePoint reversal", reversed, expected);
                }
                {
                    PmCartesian middle;
                    PmCircle first, second;
                    require(pmCirclePoint(&circle, 0.375*circle.angle, &middle) == 0,
                            "nonuniform partition point");
                    require(pmCircleInit(&first, &start, &middle, &center, &n, 0) == 0,
                            "nonuniform first constructor");
                    require(pmCircleInit(&second, &middle, &end, &center, &n, 0) == 0,
                            "nonuniform second constructor");
                    close_to("S08 posemath nonuniform partition",
                             fmax(adapter(&first), adapter(&second)), actual, TOL);
                }
            }
        }
    }
    {
        const PmCartesian u = {1,0,0}, v = {0,1,0}, n = {0,0,1};
        const int splits[] = { 2, 4, 8, 100 };
        PmCartesian start = point(center, u, v, n, 10, 0, 0);
        PmCartesian end = point(center, u, v, n, 10, 0, 20);
        PmCircle whole;
        double maximum;
        size_t split;
        require(pmCircleInit(&whole, &start, &end, &center, &n, 0) == 0,
                "whole helix constructor");
        maximum = adapter(&whole);
        for (split = 0; split < sizeof(splits)/sizeof(splits[0]); ++split) {
            int count = splits[split], part;
            for (part = 0; part < count; ++part) {
                PmCartesian p0, p1;
                PmCircle piece;
                require(pmCirclePoint(&whole, whole.angle*part/count, &p0) == 0,
                        "helix partition start");
                require(pmCirclePoint(&whole, whole.angle*(part+1)/count, &p1) == 0,
                        "helix partition end");
                require(pmCircleInit(&piece, &p0, &p1, &center, &n, 0) == 0,
                        "helix partition constructor");
                close_to("H03 posemath segmentation", adapter(&piece), maximum, TOL);
            }
        }
        {
            const int turns[] = { 2, 10, 1000 };
            size_t turn;
            for (turn = 0; turn < sizeof(turns)/sizeof(turns[0]); ++turn) {
                PmCircle multiple, reverse;
                end = point(center, u, v, n, 10, 0, 20.0*turns[turn]);
                require(pmCircleInit(&multiple, &start, &end, &center, &n,
                                    turns[turn]-1) == 0, "multiturn constructor");
                require(pmCircleInit(&reverse, &end, &start, &center, &n,
                                    -turns[turn]) == 0, "reverse multiturn constructor");
                close_to("H05 posemath multiple turns", adapter(&multiple), maximum, TOL);
                close_to("H07 posemath reverse multiple turns", adapter(&reverse), maximum, TOL);
            }
        }
        /* Trimming changes the minimum radius of this inward spiral. */
        end = point(center, u, v, n, 2, 2, 0);
        require(pmCircleInit(&whole, &start, &end, &center, &n, 0) == 0,
                "planar stretch constructor");
        require(pmCircleStretch(&whole, 1.25, 0) == 0, "planar stretch");
        close_to("I06 curvature after stretch", adapter(&whole),
                 reference_max(10, -5, 1.25, 0), TOL);
    }
}

static void invalid_adapter(PmCircle circle, const char *label)
{
    double sentinel = -1234.5;
    require(pmCircleMaxCurvature(&circle, &sentinel) != 0, label);
    require(sentinel == -1234.5, "invalid adapter preserves output");
    require(pmCircleEffectiveMinRadius(&circle) == 0,
            "invalid wrapper returns zero radius");
}

static void test_invalid_adapter(void)
{
    const PmCartesian center = {0,0,0}, start = {10,0,0}, end = {0,10,2};
    const PmCartesian normal = {0,0,1};
    PmCircle circle, bad;
    double sentinel = -1234.5;
    require(pmCircleInit(&circle, &start, &end, &center, &normal, 0) == 0,
            "invalid adapter base constructor");
    require(pmCircleMaxCurvature(NULL, &sentinel) != 0,
            "N07 null circle rejected");
    require(sentinel == -1234.5, "null circle preserves output");
    require(pmCircleMaxCurvature(&circle, NULL) != 0,
            "N07 null adapter output rejected");
#define INVALID_FIELD(field, value) do { \
        bad = circle; bad.field = (value); \
        invalid_adapter(bad, "invalid adapter " #field); \
    } while (0)
    INVALID_FIELD(radius, 0);
    INVALID_FIELD(spiral, -10);
    INVALID_FIELD(angle, 0);
    INVALID_FIELD(angle, -1);
    INVALID_FIELD(radius, NAN);
    INVALID_FIELD(spiral, INFINITY);
    INVALID_FIELD(angle, INFINITY);
    INVALID_FIELD(rHelix.x, NAN);
    INVALID_FIELD(rHelix.y, INFINITY);
    INVALID_FIELD(rHelix.z, INFINITY);
#undef INVALID_FIELD
}

/* A spiral's progress need not equal its spatial arc length. */
static void test_progress_mapping(void)
{
    const PmCartesian center = {0,0,0}, normal = {0,0,1}, start = {1,0,0};
    PmCartesian end = {20*cos(19), 20*sin(19), 190};
    PmCircle circle;
    SpiralArcLengthFit fit;
    double total, theta, step, measured, previous_error = INFINITY;
    long double expected;
    int iteration;
    require(pmCircleInit(&circle, &start, &end, &center, &normal, 3) == 0,
            "progress spiral constructor");
    close_to("progress spiral sweep", circle.angle, 19, TOL);
    require(findSpiralArcLengthFit(&circle, &fit) == 0, "production spiral fit");
    total = hypot(fit.total_planar_length, 190);
    expected = sqrtl(102)*fit.total_planar_length/(total*sqrtl(2));
    require(expected > 5.2L && expected < 5.3L, "progress counterexample magnitude");
    for (iteration = 0, step = 1e-4; iteration < 4; ++iteration, step /= 10) {
        PmCartesian p;
        double error;
        require(pmCircleAngleFromProgress(&circle, &fit, step, &theta) == 0,
                "production progress inversion");
        require(pmCirclePoint(&circle, theta, &p) == 0,
                "production progress coordinate");
        measured = hypot(hypot(p.x-start.x, p.y-start.y), p.z-start.z)/step;
        error = fabs(measured-(double)expected);
        require(error < previous_error, "progress derivative converges with step");
        previous_error = error;
    }
    close_to("I07 measured progress derivative", measured, expected, 2e-7L);
    end.x = cos(1.2); end.y = sin(1.2); end.z = 3;
    require(pmCircleInit(&circle, &start, &end, &center, &normal, 0) == 0,
            "progress helix constructor");
    require(findSpiralArcLengthFit(&circle, &fit) == 0, "helix production fit");
    total = hypot(fit.total_planar_length, 3);
    require(pmCircleAngleFromProgress(&circle, &fit, total*0.375, &theta) == 0,
            "helix progress inversion");
    close_to("helix progress is arc length", theta, circle.angle*0.375, TOL);
}

int main(void)
{
    test_helices();
    test_spirals();
    test_random_and_scale();
    test_invalid_kernel();
    test_posemath();
    test_invalid_adapter();
    test_progress_mapping();
    puts("circle curvature: OK");
    return EXIT_SUCCESS;
}
