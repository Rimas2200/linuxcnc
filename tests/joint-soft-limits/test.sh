#!/bin/bash
set -e
for configuration in corexy scara; do
    rm -f sim.var sim.var.bak
    touch sim.var
    linuxcnc -r "$configuration.ini"
done
