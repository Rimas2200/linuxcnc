#!/usr/bin/env python3
"""Reject paths whose joint limits are exceeded between legal endpoints."""

import math
import os
import time

import linuxcnc
import linuxcnc_util


ini = linuxcnc.ini(os.environ["INI_FILE_NAME"])
kinematics = ini.find("KINS", "KINEMATICS").split()[0]
joint_count = int(ini.find("KINS", "JOINTS"))
c = linuxcnc.command()
s = linuxcnc.stat()
e = linuxcnc.error_channel()
machine = linuxcnc_util.LinuxCNC(command=c, status=s, error=e)


def errors():
    messages = []
    while True:
        error = e.poll()
        if error is None:
            return messages
        messages.append(error[1])


def drain_after_abort():
    quiet_until = time.monotonic() + 0.1
    deadline = time.monotonic() + 2
    while time.monotonic() < deadline:
        if errors():
            quiet_until = time.monotonic() + 0.1
        elif time.monotonic() >= quiet_until:
            return
        time.sleep(0.01)
    raise AssertionError("error channel did not settle after abort")


def mdi(command):
    c.mdi(command)
    result = c.wait_complete(10)
    messages = errors()
    assert result == linuxcnc.RCS_DONE, (command, result, messages)
    assert not messages, (command, messages)


def rejected(command, move_type, joint, direction):
    expected = "would exceed joint {}'s {} limit".format(joint, direction)
    s.poll()
    start = tuple(s.joint_position[:joint_count])
    assert not errors()
    c.mdi(command)
    result = c.wait_complete(10)
    messages = []
    deadline = time.monotonic() + 1
    while time.monotonic() < deadline:
        messages.extend(errors())
        if any(move_type + " move" in message and expected in message
               for message in messages):
            break
        time.sleep(0.01)
    s.poll()
    position = tuple(s.joint_position[:joint_count])
    diagnostic = (command, result, start, position, messages)
    assert all(abs(a - b) < 1e-8 for a, b in zip(start, position)), diagnostic
    assert any(move_type + " move" in message and expected in message
               for message in messages), diagnostic
    print("Rejected before joint motion: " + command, flush=True)
    c.abort()
    c.wait_complete(5)
    drain_after_abort()


machine.wait_for_linuxcnc_startup()
c.state(linuxcnc.STATE_ESTOP_RESET)
c.state(linuxcnc.STATE_ON)
c.mode(linuxcnc.MODE_MANUAL)
c.wait_complete(5)
c.home(-1)
machine.wait_for_home([int(i < joint_count) for i in range(9)])
c.mode(linuxcnc.MODE_MDI)
c.wait_complete(5)
mdi("G21 G90 G17 G40 G49 G54 G61 G94 F600")

if kinematics == "corexykins":
    # All four arcs stay within the XYZ limits.  Each has identical start
    # and end points, but an interior joint extremum exceeds +/-10.
    for command, joint, direction in [
            ("G2 I4 J4", 0, "positive"),
            ("G3 I-4 J-4", 0, "negative"),
            ("G2 I4 J-4", 1, "positive"),
            ("G3 I-4 J4", 1, "negative")]:
        rejected(command, "Circular", joint, direction)

    mdi("G2 I2 J2")
    mdi("G3 I2 J2 P3 Z1")
    mdi("G1 X4 Y4 Z0")
    # The XYZ box overestimates this circle's joint range of [-8, 8].
    mdi("G2 I-4 J-4")
    # Accept tangency and a line with constant joint 0 position.
    mdi("G1 X5 Y5")
    mdi("G2 I-5 J-5")
    mdi("G1 X4 Y6")
    mdi("G1 X0 Y0")

elif kinematics == "scarakins":
    y = 5 * math.sqrt(3)
    s.poll()
    assert abs(s.position[0] - 5) < 1e-8, s.position
    assert abs(s.position[1] - y) < 1e-8, s.position
    # Home and each destination have elbow angle 120 degrees.  The line
    # passes nearer the shoulder, requiring 151 degrees in its interior.
    rejected("G1 X5 Y{:.12f}".format(-y), "Linear", 1, "positive")
    # The left semicircle also goes inside the permitted elbow envelope.
    rejected("G3 X5 Y{:.12f} J{:.12f}".format(-y, -y),
             "Circular", 1, "positive")
    # The opposite semicircle stays between elbow angles 93 and 120.
    mdi("G2 X5 Y{:.12f} J{:.12f}".format(-y, -y))
    mdi("G3 X5 Y{:.12f} J{:.12f}".format(y, y))
    mdi("G1 X6 Y{:.12f}".format(y))
    mdi("G1 X5 Y{:.12f}".format(y))
else:
    raise AssertionError("unexpected kinematics: " + kinematics)

print(kinematics + " joint soft-limit checks passed", flush=True)
