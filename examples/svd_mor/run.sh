#!/bin/bash

echo "running example $(pwd)"

workdir=$(pwd)
variant="debug"
variant="release"

mkdir -p build_${variant}
cd build_${variant}

# remove old output data
rm -rf out

./shorten_implicit ../settings.py
python ../check_results.py  

cd $workdir

