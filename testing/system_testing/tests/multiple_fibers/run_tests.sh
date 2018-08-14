#!/bin/bash

echo "running example $(pwd)"

workdir=$(pwd)
variant="debug"
variant="release"

. input/copy.sh

mkdir -p build_${variant}
cd build_${variant}

# remove old output data
rm -rf out

export OMP_NUM_THREADS=1

# arguments: <fibre no> <cellml_file> <end_time>
./single_fibre ../single_fibre_settings.py 1 "../input/hodgkin_huxley_1952.c" 10

cd $workdir
