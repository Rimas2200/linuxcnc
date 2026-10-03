#!/bin/sh
set -eu
trap 'rm -f test' EXIT
gcc -O2 -Wall -Wextra -DULAPI -I"$HEADERS" \
    -I../../../src/emc/kinematics test.c -lm -o test
./test
