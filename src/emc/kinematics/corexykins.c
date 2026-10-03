/********************************************************************
* Description: kinematics for corexy
* Adapted from trivkins.c
* ref: http://corexy.com/theory.html
********************************************************************/

#include "motion.h"
#include "hal.h"
#include "rtapi.h"
#include "rtapi.h"
#include "rtapi_app.h"
#include "rtapi_math.h"
#include "rtapi_string.h"
#include "kinematics.h"

static struct data {
    hal_s32_t joints[EMCMOT_MAX_JOINTS];
} *data;

int kinematicsForward(const double *joints
                     ,EmcPose *pos
                     ,const KINEMATICS_FORWARD_FLAGS *fflags
                     ,KINEMATICS_INVERSE_FLAGS *iflags
                     ) {
    pos->tran.x = 0.5 * (joints[0] + joints[1]);
    pos->tran.y = 0.5 * (joints[0] - joints[1]);
    pos->tran.z = joints[2];
    pos->a      = joints[3];
    pos->b      = joints[4];
    pos->c      = joints[5];
    pos->u      = joints[6];
    pos->v      = joints[7];
    pos->w      = joints[8];

    return 0;
}

int kinematicsInverse(const EmcPose *pos
                     ,double *joints
                     ,const KINEMATICS_INVERSE_FLAGS *iflags
                     ,KINEMATICS_FORWARD_FLAGS *fflags
                     ) {
    joints[0] = pos->tran.x + pos->tran.y;
    joints[1] = pos->tran.x - pos->tran.y;
    joints[2] = pos->tran.z;
    joints[3] = pos->a;
    joints[4] = pos->b;
    joints[5] = pos->c;
    joints[6] = pos->u;
    joints[7] = pos->v;
    joints[8] = pos->w;

    return 0;
}

int kinematicsHome(EmcPose *world
                  ,double *joint
                  ,KINEMATICS_FORWARD_FLAGS *fflags
                  ,KINEMATICS_INVERSE_FLAGS *iflags
                  ) {
    *fflags = 0;
    *iflags = 0;
    return kinematicsForward(joint, world, fflags, iflags);
}

KINEMATICS_TYPE kinematicsType() { return KINEMATICS_BOTH; }

/* Bound X+Y and X-Y directly. Summing separate X/Y bounds can reject
   valid arcs tangent to a joint limit. */
int kinematicsInverseBounds(const KINEMATICS_PATH *path,
                            double *joint_lower, double *joint_upper,
                            int num_joints,
                            const KINEMATICS_INVERSE_FLAGS *iflags)
{
    double lower[9], upper[9];
    int i;
    if (!path || !joint_lower || !joint_upper || num_joints < 1 || num_joints > 9)
        return KINEMATICS_BOUNDS_UNKNOWN;

    if (path->circle) {
        const PmCartesian sum = {1, 1, 0}, difference = {1, -1, 0}, z = {0, 0, 1};
        if (pmCircleProjectionBounds(path->circle, &sum, &lower[0], &upper[0])
            || pmCircleProjectionBounds(path->circle, &difference, &lower[1], &upper[1])
            || pmCircleProjectionBounds(path->circle, &z, &lower[2], &upper[2]))
            return KINEMATICS_BOUNDS_UNKNOWN;
    } else {
        double start = path->start.tran.x + path->start.tran.y;
        double end = path->end.tran.x + path->end.tran.y;
        lower[0] = fmin(start, end);
        upper[0] = fmax(start, end);
        start = path->start.tran.x - path->start.tran.y;
        end = path->end.tran.x - path->end.tran.y;
        lower[1] = fmin(start, end);
        upper[1] = fmax(start, end);
        lower[2] = path->lower.tran.z; upper[2] = path->upper.tran.z;
    }
    lower[3] = path->lower.a; upper[3] = path->upper.a;
    lower[4] = path->lower.b; upper[4] = path->upper.b;
    lower[5] = path->lower.c; upper[5] = path->upper.c;
    lower[6] = path->lower.u; upper[6] = path->upper.u;
    lower[7] = path->lower.v; upper[7] = path->upper.v;
    lower[8] = path->lower.w; upper[8] = path->upper.w;
    for (i = 0; i < num_joints; i++) {
        if (!isfinite(lower[i]) || !isfinite(upper[i]) || lower[i] > upper[i])
            return KINEMATICS_BOUNDS_UNKNOWN;
    }
    for (i = 0; i < num_joints; i++) {
        joint_lower[i] = lower[i];
        joint_upper[i] = upper[i];
    }
    return KINEMATICS_BOUNDS_OK;
}

KINS_NOT_SWITCHABLE
EXPORT_SYMBOL(kinematicsType);
EXPORT_SYMBOL(kinematicsForward);
EXPORT_SYMBOL(kinematicsInverse);
EXPORT_SYMBOL(kinematicsInverseBounds);
MODULE_LICENSE("GPL");

static int comp_id;
int rtapi_app_main(void) {
    comp_id = hal_init("corexykins");
    if(comp_id < 0) return comp_id;

    data = hal_malloc(sizeof(struct data));

    hal_ready(comp_id);
    return 0;
}

void rtapi_app_exit(void) { hal_exit(comp_id); }
