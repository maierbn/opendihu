#!/bin/bash

echo "running example $(pwd)"

workdir=$(pwd)
#variant="debug"
variant="release"

mkdir -p build_${variant}
cd build_${variant}

# remove old output data
rm -rf out
ln -s /data/scratch/maierbn/out out

export OMP_NUM_THREADS=1

# loop over fibres
for fibre_no in `seq 1 50`; do

  # command arguments: <fibre_no>
  ./single_fibre ../single_fibre_settings.py $fibre_no &
  
done

./single_fibre ../single_fibre_settings.py 0

cd $workdir
