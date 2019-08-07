#!/bin/bash

echo "running example $(pwd)"

workdir=$(pwd)
variant="debug"
#variant="release"

#mkdir -p build_${variant}
scons BUILD_TYPE=${variant}
cd build_${variant}

# remove old output data
rm -rf out
rm -rf out_snapshots

./hodgkin_huxley_godunov ../settings_hodgkin_huxley_godunov.py

mkdir -p out_snapshots
cd out
cp godunov* ../out_snapshots
cd ..

python ../scripts/snapshots.py --path ./out_snapshots/

./dmdexample  

cd ..

