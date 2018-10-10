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
python ../write_csv.py
./shorten_pod ../settings_pod.py

  

cd $workdir

