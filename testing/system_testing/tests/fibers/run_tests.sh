#!/bin/bash

echo "running example $(pwd)"

workdir=$(pwd)
variant="debug"
#variant="release"

mkdir -p build_${variant}
cd build_${variant}

# remove old output data
rm -rf out

# command arguments: <name>
./laplace_structured_linear ../settings.py laplace3d_structured_linear
mv streamlines.csv laplace3d_structured_linear.csv
../scripts/csv2mesh.py laplace3d_structured_linear.csv laplace3d_structured_linear.stl
../scripts/read_csv_fibres.py laplace3d_structured_linear.csv laplace3d_structured_linear

./laplace_structured_quadratic ../settings.py laplace3d_structured_quadratic
mv streamlines.csv laplace3d_structured_quadratic.csv
../scripts/csv2mesh.py laplace3d_structured_quadratic.csv laplace3d_structured_quadratic.stl

#./laplace_structured_hermite ../settings.py laplace3d_structured_hermite
#mv streamlines.csv laplace3d_structured_hermite.csv
#../scripts/csv2mesh.py laplace3d_structured_hermite.csv laplace3d_structured_hermite.stl

./laplace_unstructured_linear ../settings.py laplace3d_unstructured_linear
mv streamlines.csv laplace3d_unstructured_linear.csv
../scripts/csv2mesh.py laplace3d_unstructured_linear.csv laplace3d_unstructured_linear.stl

./laplace_unstructured_quadratic ../settings.py laplace3d_unstructured_quadratic
mv streamlines.csv laplace3d_unstructured_quadratic.csv
../scripts/csv2mesh.py laplace3d_unstructured_quadratic.csv laplace3d_unstructured_quadratic.stl

#./laplace_unstructured_hermite ../settings.py laplace3d_unstructured_hermite
#mv streamlines.csv laplace3d_unstructured_hermite.csv
#../scripts/csv2mesh.py laplace3d_unstructured_hermite.csv laplace3d_unstructured_hermite.stl


cd $workdir
