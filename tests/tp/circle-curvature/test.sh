#!/bin/bash
set -eu

here=$(cd -- "$(dirname -- "$0")" && pwd)
root=$(cd -- "$here/../../.." && pwd)
build=$(mktemp -d)
trap 'rm -rf -- "$build"' EXIT

# Drop unused controller functions at link time.
set -f
# shellcheck disable=SC2086
"${CC:-cc}" ${CFLAGS:--O2 -g} -std=gnu11 -Wall -Wextra -fno-fast-math \
    -ffunction-sections -fdata-sections -DULAPI \
    -I"$root/src/emc/tp" -I"$root/src/libposemath" -I"$root/src/rtapi" \
    -I"$root/src/emc/nml_intf" -I"$root/src/emc/motion" \
    -I"$root/src/emc/kinematics" -I"$root/src/hal" -I"$root/include" \
    "$here/test.c" "$root/src/emc/tp/circle_curvature.c" \
    "$root/src/emc/tp/blendmath.c" "$root/src/libposemath/_posemath.c" \
    ${LDFLAGS:-} -Wl,--gc-sections -lm -o "$build/test"
"$build/test"
