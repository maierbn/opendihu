#!/bin/bash
# Script to generate cuboid meshes as replacement for tendon and muscle, to test the program and coupling schemes.

# usage: generate_cuboid_mesh.py <output filename> <n points x> <n points y> <n points z> <mesh width x> <mesh width y> <mesh width z> [<offset x> <offset y> <offset z>]

# generate mesh for tendon
generate_cuboid_mesh.py tendon_box.bin   11 11 22   0.3 0.3 0.6    -2.0 -2.0 -0.3

# generate mesh for muscle 1 and 2
generate_cuboid_mesh.py muscle1_box.bin  11 11 33   0.3 0.3 0.9    -2.0 -2.0 -1.2
generate_cuboid_mesh.py muscle2_box.bin  11 11 33   0.3 0.3 0.9    -2.0 -2.0  0.3

# examine generated files
echo ""
echo ""
echo "tendon_box.bin"
echo "--------------"
examine_bin_fibers.py tendon_box.bin

echo ""
echo ""
echo "muscle1_box.bin"
echo "--------------"
examine_bin_fibers.py muscle1_box.bin

echo ""
echo ""
echo "muscle2_box.bin"
echo "--------------"
examine_bin_fibers.py muscle2_box.bin