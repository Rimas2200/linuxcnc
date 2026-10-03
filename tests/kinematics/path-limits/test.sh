#!/bin/sh
set -eu
trap 'rm -f test' EXIT
gcc -O2 -Wall -Wextra -DULAPI -I"$HEADERS" \
    -I../../../src/emc/kinematics -I../../../src/emc/motion \
    -I../../../src/emc/tp test.c \
    ../../../src/emc/motion/kinematics_limits.c \
    ../../../src/emc/tp/circle_progress.c \
    -L"$LIBDIR" -Wl,-rpath,"$LIBDIR" -lposemath -lm -o test
./test
