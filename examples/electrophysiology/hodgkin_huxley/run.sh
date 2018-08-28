#!/bin/bash

echo "running example $(pwd)"

workdir=$(pwd)
variant="debug"
variant="release"

mkdir -p build_${variant}
cd build_${variant}

# remove old output data
rm -rf out

# loop over fibres

./hodgkin_huxley_explicit ../settings.py
  

cd $workdir

