#!/bin/bash

echo "running example $(pwd)"

export OMP_NUM_THREADS=1

workdir=$(pwd)
variant="debug"
#variant="release_with_debug_info"
variant="release"

# generate mesh from biceps stl file
#cd $workdir/scripts
#. run.sh
cd $workdir

mkdir -p build_${variant}
cd build_${variant}

# remove old output data
rm -rf out



# command arguments: <name>
./laplace_structured_linear ../settings.py laplace3d_structured_linear
python ../scripts/csv2mesh.py laplace3d_structured_linear.csv laplace3d_structured_linear.stl
../scripts/read_csv_fibres.py laplace3d_structured_linear.csv laplace3d_structured_linear
../scripts/read_csv_fibres.py laplace3d_structured_linear_raw.csv laplace3d_structured_linear_raw
meshlab laplace3d_structured_linear.stl &

cd $workdir
