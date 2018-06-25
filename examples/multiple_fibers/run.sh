#!/bin/bash

echo "running example $(pwd)"

workdir=$(pwd)
#variant="debug"
variant="release"

mkdir -p build_${variant}
cd build_${variant}

# remove old output data
rm -rf out

# loop over fibres
for fibre_no in `seq 1 100`; do

  # command arguments: <fibre_no>
  ./single_fibre ../single_fibre_settings.py $fibre_no &
  
done

cd $workdir
