#!/bin/bash
# Script to generate cuboid meshes as replacement for tendon and muscle, to test the program and coupling schemes.

# usage: generate_cuboid_mesh.py <output filename> <n points x> <n points y> <n points z> <mesh width x> <mesh width y> <mesh width z> [<offset x> <offset y> <offset z>]

# generate mesh for tendon
generate_cuboid_mesh.py tendon_box.bin 9 9 21     0.5 0.5 0.3     -2.0 -2.0 -6
#generate_cuboid_mesh.py tendon_box.bin 5 5 13     0.2 0.2 0.6     -2.0 -2.0 -6

# generate mesh for muscle
generate_cuboid_mesh.py muscle_box.bin 5 5 1481   1.0 1.0 0.01    -2.0 -2.0 0.0

# examine generated files
echo ""
echo ""
echo "tendon_box.bin"
echo "--------------"
examine_bin_fibers.py tendon_box.bin

echo ""
echo ""
echo "muscle_box.bin"
echo "--------------"
examine_bin_fibers.py muscle_box.bin
