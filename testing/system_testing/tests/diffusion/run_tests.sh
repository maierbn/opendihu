#!/bin/bash

echo "running example $(pwd)"

workdir=$(pwd)
variant="debug"
variant="release"

mkdir -p build_${variant}
cd build_${variant}

# remove old output data
rm -rf out

# command arguments: <name> <number elements>

./diffusion_1d_structured_regular_fixed_linear    ../settings_1d.py 1d_structured_regular_fixed_linear    10
./diffusion_1d_structured_regular_fixed_quadratic ../settings_1d.py 1d_structured_regular_fixed_quadratic 10
./diffusion_1d_structured_regular_fixed_hermite   ../settings_1d.py 1d_structured_regular_fixed_hermite   6   # has lower number of elements
./diffusion_1d_structured_deformable_linear       ../settings_1d.py 1d_structured_deformable_linear       10
./diffusion_1d_structured_deformable_quadratic    ../settings_1d.py 1d_structured_deformable_quadratic    6
./diffusion_1d_structured_deformable_hermite      ../settings_1d.py 1d_structured_deformable_hermite      6
./diffusion_1d_unstructured_deformable_linear     ../settings_1d.py 1d_unstructured_deformable_linear     10
./diffusion_1d_unstructured_deformable_quadratic  ../settings_1d.py 1d_unstructured_deformable_quadratic  6
./diffusion_1d_unstructured_deformable_hermite    ../settings_1d.py 1d_unstructured_deformable_hermite    6


./diffusion_2d_structured_regular_fixed_linear    ../settings_2d.py 2d_structured_regular_fixed_linear    10
./diffusion_2d_structured_regular_fixed_quadratic ../settings_2d.py 2d_structured_regular_fixed_quadratic 10
./diffusion_2d_structured_regular_fixed_hermite   ../settings_2d.py 2d_structured_regular_fixed_hermite   6   # has lower number of elements
./diffusion_2d_structured_deformable_linear       ../settings_2d.py 2d_structured_deformable_linear       10
./diffusion_2d_structured_deformable_quadratic    ../settings_2d.py 2d_structured_deformable_quadratic    6
./diffusion_2d_structured_deformable_hermite      ../settings_2d.py 2d_structured_deformable_hermite      6
./diffusion_2d_unstructured_deformable_linear     ../settings_2d.py 2d_unstructured_deformable_linear     10
./diffusion_2d_unstructured_deformable_quadratic  ../settings_2d.py 2d_unstructured_deformable_quadratic  6
./diffusion_2d_unstructured_deformable_hermite    ../settings_2d.py 2d_unstructured_deformable_hermite    6

./diffusion_3d_structured_regular_fixed_linear    ../settings_3d.py 3d_structured_regular_fixed_linear    10
./diffusion_3d_structured_regular_fixed_quadratic ../settings_3d.py 3d_structured_regular_fixed_quadratic 10
./diffusion_3d_structured_regular_fixed_hermite   ../settings_3d.py 3d_structured_regular_fixed_hermite   6   # has lower number of elements
./diffusion_3d_structured_deformable_linear       ../settings_3d.py 3d_structured_deformable_linear       10
./diffusion_3d_structured_deformable_quadratic    ../settings_3d.py 3d_structured_deformable_quadratic    6
./diffusion_3d_structured_deformable_hermite      ../settings_3d.py 3d_structured_deformable_hermite      6
./diffusion_3d_unstructured_deformable_linear     ../settings_3d.py 3d_unstructured_deformable_linear     10
./diffusion_3d_unstructured_deformable_quadratic  ../settings_3d.py 3d_unstructured_deformable_quadratic  6
./diffusion_3d_unstructured_deformable_hermite    ../settings_3d.py 3d_unstructured_deformable_hermite    6

cd $workdir
