#!/bin/bash

echo "running example $(pwd)"

export OMP_NUM_THREADS=1

workdir=$(pwd)
variant="debug"
#variant="release_with_debug_info"
variant="release"

# generate mesh from biceps stl file
cd $workdir/scripts
. run.sh
cd $workdir

mkdir -p build_${variant}
cd build_${variant}

# remove old output data
rm -rf out



# command arguments: <name>
./laplace_structured_linear ../settings.py laplace3d_structured_linear
../scripts/csv2mesh.py laplace3d_structured_linear.csv laplace3d_structured_linear.stl
../scripts/read_csv_fibres.py laplace3d_structured_linear.csv laplace3d_structured_linear
../scripts/read_csv_fibres.py laplace3d_structured_linear_raw.csv laplace3d_structured_linear_raw


./laplace_structured_quadratic ../settings.py laplace3d_structured_quadratic
../scripts/csv2mesh.py laplace3d_structured_quadratic.csv laplace3d_structured_quadratic.stl
../scripts/read_csv_fibres.py laplace3d_structured_quadratic.csv laplace3d_structured_quadratic

#./laplace_structured_hermite ../settings.py laplace3d_structured_hermite
#../scripts/csv2mesh.py laplace3d_structured_hermite.csv laplace3d_structured_hermite.stl

./laplace_unstructured_linear ../settings.py laplace3d_unstructured_linear
../scripts/csv2mesh.py laplace3d_unstructured_linear.csv laplace3d_unstructured_linear.stl
../scripts/read_csv_fibres.py laplace3d_unstructured_linear.csv laplace3d_unstructured_linear

./laplace_unstructured_quadratic ../settings.py laplace3d_unstructured_quadratic
../scripts/csv2mesh.py laplace3d_unstructured_quadratic.csv laplace3d_unstructured_quadratic.stl
../scripts/read_csv_fibres.py laplace3d_unstructured_quadratic.csv laplace3d_unstructured_quadratic

#./laplace_unstructured_hermite ../settings.py laplace3d_unstructured_hermite
#../scripts/csv2mesh.py laplace3d_unstructured_hermite.csv laplace3d_unstructured_hermite.stl


cd $workdir
