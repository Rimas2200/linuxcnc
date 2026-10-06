#!/usr/bin/env python3
"""Headless motion smoke test"""
import csv
import math
from pathlib import Path
import time

import linuxcnc


command = linuxcnc.command()
status = linuxcnc.stat()
errors = linuxcnc.error_channel()
Path("result.log").unlink(missing_ok=True)


def complete():
    result = command.wait_complete(10)
    if result not in (linuxcnc.RCS_DONE,):
        raise RuntimeError(f"command failed or timed out: {result}")


def check_errors():
    while True:
        error = errors.poll()
        if error is None:
            return
        if error[0] in (linuxcnc.NML_ERROR, linuxcnc.OPERATOR_ERROR):
            raise RuntimeError(error[1])


def move(gcode, endpoint, samples, writer):
    command.mdi(gcode)
    deadline = time.monotonic() + 20
    while True:
        result = command.wait_complete(0.01)
        status.poll()
        check_errors()
        xyz = status.position[:3]
        if not all(math.isfinite(value) for value in xyz):
            raise RuntimeError(f"nonfinite position: {xyz}")
        writer.writerow((samples[0], time.monotonic(), *xyz))
        samples[0] += 1
        if result == linuxcnc.RCS_DONE:
            break
        if result == linuxcnc.RCS_ERROR or time.monotonic() > deadline:
            raise RuntimeError(f"move failed or timed out: {gcode}")
    if max(abs(a - b) for a, b in zip(xyz, endpoint)) > 1e-4:
        raise RuntimeError(f"wrong endpoint for {gcode}: {xyz}, expected {endpoint}")


try:
    command.state(linuxcnc.STATE_ESTOP_RESET)
    complete()
    command.state(linuxcnc.STATE_ON)
    complete()
    command.mode(linuxcnc.MODE_MDI)
    complete()
    with open("positions.csv", "w", newline="") as output:
        writer = csv.writer(output)
        writer.writerow(("sample", "monotonic_seconds", "x", "y", "z"))
        samples = [0]
        move("G21 G90 G17 G61 F1200 G1 X10 Y0 Z0", (10, 0, 0), samples, writer)
        move("G3 X10 Y0 I-10 J0", (10, 0, 0), samples, writer)
        move("G3 X10 Y0 Z20 I-10 J0", (10, 0, 20), samples, writer)
        move("G2 X10 Y0 Z0 I-10 J0", (10, 0, 0), samples, writer)
        # Equivalent one-turn helix as four commands.
        for block, endpoint in (
            ("G3 X0 Y10 Z5 I-10 J0", (0, 10, 5)),
            ("G3 X-10 Y0 Z10 I0 J-10", (-10, 0, 10)),
            ("G3 X0 Y-10 Z15 I10 J0", (0, -10, 15)),
            ("G3 X10 Y0 Z20 I0 J10", (10, 0, 20)),
        ):
            move(block, endpoint, samples, writer)
        # A partial helix with a sweep below one radian.
        x, y = 10 * math.cos(0.5), 10 * math.sin(0.5)
        move(f"G3 X{x:.12f} Y{y:.12f} Z21 I-10 J0", (x, y, 21), samples, writer)
        if samples[0] < 20:
            raise RuntimeError("insufficient motion samples")
    Path("result.log").write_text("helix motion: OK\n")
finally:
    command.abort()
    command.wait_complete(5)
    command.state(linuxcnc.STATE_OFF)
    command.wait_complete(5)
