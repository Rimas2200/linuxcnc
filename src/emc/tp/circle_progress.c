/********************************************************************
 * Description: circle_progress.c
 *
 * Circular and spiral arc progress parameterization, shared by motion
 * checks and the trajectory planner.
 *
 * Author: Robert W. Ellenberg
 * License: GPL Version 2
 *
 * Copyright (c) 2014 All rights reserved.
 ********************************************************************/

#include "posemath.h"
#include "tc_types.h"
#include "tp_types.h"
#include "rtapi_math.h"
#include "blendmath.h"
#include "tp_debug.h"

/** @section spiralfuncs Functions to approximate spiral arc length */

/**
 * Intermediate function to find the angle for a parameter from 0..1 along the
 * spiral arc.
 */
int pmCircleAngleFromParam(PmCircle const * const circle,
        SpiralArcLengthFit const * const fit,
        double t,
        double * const angle)
{
    if (fit->spiral_in) {
        t = 1.0 - t;
    }
    //TODO error or cleanup input to prevent param outside 0..1
    double s_in = t * fit->total_planar_length;

    // Quadratic formula to invert arc length -> angle

    double A = fit->b0;
    double B = fit->b1;
    double C = -s_in;

    double disc = pmSq(B) - 4.0 * A * C ;
    if (disc < 0) {
        rtapi_print_msg(RTAPI_MSG_ERR, "discriminant %f is negative in angle calculation\n",disc);
        return TP_ERR_FAIL;
    }

    /*
     * Stability of inverting the arc-length relationship.
     * Since the b1 coefficient is analogous to arc radius, we can be
     * reasonably assured that it will be large enough not to cause numerical
     * errors. If this is not the case, then the arc itself is degenerate (very
     * small radius), and this condition should be caught well before here.
     *
     * Since an arc with a very small spiral coefficient will have a small b0
     * coefficient in the fit, we use the Citardauq Formula to ensure that the
     * positive root does not lose precision due to subtracting near-similar values.
     *
     * For more information, see:
     * http://people.csail.mit.edu/bkph/articles/Quadratics.pdf
     */

    double angle_out = (2.0 * C) / ( -B - pmSqrt(disc));

    if (fit->spiral_in) {
        // Spiral fit assumes that we're spiraling out, so
        // parameterize from opposite end
        angle_out = circle->angle - angle_out;
    }

    *angle = angle_out;
    return TP_ERR_OK;
}


static void printSpiralArcLengthFit(SpiralArcLengthFit const * const fit)
{
    (void)fit; /* Debug printing can be compiled out. */
    tp_debug_print("Spiral fit: b0 = %.12f, b1 = %.12f, length = %.12f, spiral_in = %d\n",
            fit->b0,
            fit->b1,
            fit->total_planar_length,
            fit->spiral_in);
}

/**
 * Approximate the arc length function of a general spiral.
 *
 * The closed-form arc length of a general archimedean spiral is rather
 * computationally messy to work with.
 * See http://mathworld.wolfram.com/ArchimedesSpiral.html for the actual form.
 *
 * The simplification here is made possible by a few assumptions:
 *  1) That the spiral starts with a nonzero radius
 *  2) The spiral coefficient (i.e. change in radius / angle) is not too large
 *  3) The spiral coefficient has some minimum magnitude ("perfect" circles are handled as a special case)
 *
 * The 2nd-order fit below works by matching slope at the start and end of the
 * arc length vs. angle curve. This completely specifies the 2nd order fit.
 * Also, this fit predicts a total arc length >= the true arc length, which
 * means the true speed along the curve will be the same or slower than the
 * nominal speed.
 */
int findSpiralArcLengthFit(PmCircle const * const circle,
        SpiralArcLengthFit * const fit)
{
    // Additional data for arc length approximation
    double spiral_coef = circle->spiral / circle->angle;
    double min_radius = circle->radius;

    if (circle->spiral < 0.0) {
        // Treat as positive spiral, parameterized in opposite
        // direction
        spiral_coef*=-1.0;
        // Treat final radius as starting radius for fit, so we add the
        // negative spiral term to get the minimum radius
        //
        min_radius+=circle->spiral;
        fit->spiral_in = true;
    } else {
        fit->spiral_in = false;
    }
    tp_debug_print("radius = %.12f, angle = %.12f\n", min_radius, circle->angle);
    tp_debug_print("spiral_coef = %.12f\n", spiral_coef);


    //Compute the slope of the arc length vs. angle curve at the start and end of the segment
    double slope_start = pmSqrt(pmSq(min_radius) + pmSq(spiral_coef));
    double slope_end = pmSqrt(pmSq(min_radius + spiral_coef * circle->angle) + pmSq(spiral_coef));

    fit->b0 = (slope_end - slope_start) / (2.0 * circle->angle);
    fit->b1 = slope_start;

    fit->total_planar_length = fit->b0 * pmSq(circle->angle) + fit->b1 * circle->angle;
    printSpiralArcLengthFit(fit);

    // Check against start and end angle
    double angle_end_chk = 0.0;
    int res_angle = pmCircleAngleFromParam(circle, fit, 1.0, &angle_end_chk);
    if (res_angle != TP_ERR_OK) {
        //TODO better error message
        rtapi_print_msg(RTAPI_MSG_ERR,
                "Spiral fit failed\n");
        return TP_ERR_FAIL;
    }

    // Check fit against angle
    double fit_err = angle_end_chk - circle->angle;
    if (fabs(fit_err) > TP_ANGLE_EPSILON) {
        rtapi_print_msg(RTAPI_MSG_ERR,
                "Spiral fit angle difference is %e, maximum allowed is %e\n",
                fit_err,
                TP_ANGLE_EPSILON);
        return TP_ERR_FAIL;
    }

    return TP_ERR_OK;
}


/**
 * Compute the angle around a circular segment from the total progress along
 * the curve.
 */
int pmCircleAngleFromProgress(PmCircle const * const circle,
        SpiralArcLengthFit const * const fit,
        double progress,
        double * const angle)
{
    double h2;
    pmCartMagSq(&circle->rHelix, &h2);
    double s_end = pmSqrt(pmSq(fit->total_planar_length) + h2);
    // Parameterize by total progress along helix
    double t = progress / s_end;
    return pmCircleAngleFromParam(circle, fit, t, angle);
}
