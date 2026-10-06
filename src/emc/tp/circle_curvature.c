/* SPDX-License-Identifier: GPL-2.0-only */
#include <rtapi_math.h>
#include "circle_curvature.h"

/* rtapi_math lacks hypot on some realtime targets. Inputs are nonnegative. */
static double norm2(double x, double y)
{
    double hi = fmax(x, y);
    double lo = fmin(x, y);
    if (hi == 0.0) {
        return 0.0;
    }
    double ratio = lo / hi;
    return hi * sqrt(1.0 + ratio * ratio);
}

static int curvatureAtRadius(double radius, double k, double h, double *value)
{
    double scale = fmax(radius, fmax(k, h));
    double r = radius / scale;
    double c = k / scale;
    double z = h / scale;
    double speed = sqrt(r * r + c * c + z * z);
    double cross = norm2(r * r + 2.0 * c * c, z * norm2(r, 2.0 * c));
    double curvature = (cross / (speed * speed * speed)) / scale;
    if (!isfinite(curvature) || curvature <= 0.0) {
        return -1;
    }
    *value = curvature;
    return 0;
}

int tpCircleMaxCurvature(double radius0, double radial_delta, double sweep,
                        double axial_displacement, double *kappa_max)
{
    if (!kappa_max || !isfinite(radius0) || !isfinite(radial_delta) ||
            !isfinite(sweep) || !isfinite(axial_displacement) ||
            radius0 <= 0.0 || sweep <= 0.0) {
        return -1;
    }
    double radius1 = radius0 + radial_delta;
    double k = fabs(radial_delta / sweep);
    double h = fabs(axial_displacement / sweep);
    if (!isfinite(radius1) || radius1 <= 0.0 || !isfinite(k) || !isfinite(h) ||
            (radial_delta != 0.0 && k == 0.0) ||
            (axial_displacement != 0.0 && h == 0.0)) {
        return -1;
    }

    double maximum, end;
    if (curvatureAtRadius(radius0, k, h, &maximum) ||
            curvatureAtRadius(radius1, k, h, &end)) {
        return -1;
    }
    maximum = fmax(maximum, end);

    /* With x=R^2, a=k^2, b=h^2, the derivative has the sign of
     * (b-8a)(b+a)-6ax-x^2. Any positive root is a maximum.
     * Rationalize the root to avoid cancellation near b=8a. */
    if (radial_delta != 0.0 && k < h) {
        double q = k / h;
        double a = q * q;
        if (a < 0.125) {
            double root = h * sqrt(((1.0 - 8.0 * a) * (1.0 + a)) /
                                  (sqrt(a * a - 7.0 * a + 1.0) + 3.0 * a));
            if (root > fmin(radius0, radius1) && root < fmax(radius0, radius1)) {
                double interior;
                if (curvatureAtRadius(root, k, h, &interior)) {
                    return -1;
                }
                maximum = fmax(maximum, interior);
            }
        }
    }
    *kappa_max = maximum;
    return 0;
}
