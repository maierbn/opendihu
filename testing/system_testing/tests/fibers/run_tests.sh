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
./laplace3d    ../settings.py laplace3d

cd $workdir
