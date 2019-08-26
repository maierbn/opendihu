#!/bin/bash

echo "running example $(pwd)"

workdir=$(pwd)
variant="debug"
#variant="release"

scons BUILD_TYPE=${variant}
cd build_${variant}

# remove old output data
rm -rf out

./hodgkin_huxley_godunov ../settings_hodgkin_huxley_godunov.py
python ../scripts/snapshots.py --path ./out/
./svdexample 

cd $workdir

