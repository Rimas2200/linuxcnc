/* Bounded preflight checks for coordinated paths.  License: GPL Version 2. */
#ifndef KINEMATICS_LIMITS_H
#define KINEMATICS_LIMITS_H

#include "kinematics.h"

enum {
    KINEMATICS_PATH_OK = 0,
    KINEMATICS_PATH_UNSUPPORTED = 1,
    KINEMATICS_PATH_LIMIT = -1,
    KINEMATICS_PATH_UNCERTAIN = -2,
    KINEMATICS_PATH_INVALID = -3
};

enum {
    KINEMATICS_PATH_MAX_DEPTH = 20,
    KINEMATICS_PATH_MAX_BOXES = 64
};

/* Check the whole path against active joint limits. UNCERTAIN means the
 * subdivision budget was exhausted and the move must be rejected.
 * LIMIT sets failing_joint and direction (-1 lower, +1 upper).
 * Does not change flags or kinematics state. */
int kinematicsCheckPath(const EmcPose *start, const EmcPose *end,
                       const PmCircle *circle, KINEMATICS_INVERSE_BOUNDS bounds,
                       const double *joint_min, const double *joint_max,
                       int num_joints, unsigned active_mask,
                       const KINEMATICS_INVERSE_FLAGS *iflags,
                       int *failing_joint, int *direction);

#endif
