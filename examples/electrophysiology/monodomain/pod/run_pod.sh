#!/bin/bash

echo "running example $(pwd)"

workdir=$(pwd)
variant="debug"
#variant="release"

scons BUILD_TYPE=${variant}
cd build_${variant}

# remove old output data
rm -rf out
rm  -rf out_snapshots

# creating snapshots
./hodgkin_huxley_godunov ../settings_hodgkin_huxley_godunov.py

# creating the snapshots matrix
#mkdir out_snapshots
#cp ./out/*.py ./out_snapshots
#python ../../scripts/snapshots.py ./out_snapshots

# run the reduced-order scheme
#./hodgkin_huxley_pod ../settings_hodgkin_huxley_pod.py
  
cd ${workdir}

