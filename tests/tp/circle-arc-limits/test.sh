#!/bin/sh
set -eu

# Drop unused TP functions at link time.
TOPDIR=$(dirname "$HEADERS")
trap 'rm -f test-circle-arc-limits' EXIT HUP INT TERM
set -f
# shellcheck disable=SC2086
"${CC:-cc}" ${CFLAGS:--O2} -std=gnu11 -Wall -Wextra -fno-fast-math \
    -ffunction-sections -fdata-sections -DULAPI \
    -I"$HEADERS" -I"$TOPDIR/src" -I"$TOPDIR/src/emc" -I"$TOPDIR/src/emc/tp" \
    test.c "$TOPDIR/src/emc/tp/tc.c" "$TOPDIR/src/emc/tp/blendmath.c" \
    "$TOPDIR/src/emc/tp/circle_curvature.c" "$TOPDIR/src/emc/tp/cruckig/roots.c" \
    -L"$LIBDIR" -Wl,-rpath,"$LIBDIR" ${LDFLAGS:-} -Wl,--gc-sections -lposemath -lm \
    -o test-circle-arc-limits
./test-circle-arc-limits
