/* SPDX-License-Identifier: GPL-2.0-only */
#ifndef TP_CIRCLE_CURVATURE_H
#define TP_CIRCLE_CURVATURE_H

/* Maximum curvature for R(theta)=radius0+radial_delta*theta/sweep,
 * with axial_displacement over the full sweep. radial_delta is not the
 * final radius. Both endpoint radii and sweep must be positive.
 * Returns -1 on invalid or unrepresentable input/result and leaves
 * *kappa_max unchanged. The result is not a directed-rounding bound.
 */
int tpCircleMaxCurvature(double radius0, double radial_delta, double sweep,
                        double axial_displacement, double *kappa_max);

#endif
