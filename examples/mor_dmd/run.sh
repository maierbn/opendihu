#!/bin/bash

echo "running example $(pwd)"

workdir=$(pwd)
#variant="debug"
variant="release"

#mkdir -p build_${variant}
scons BUILD_TYPE=${variant}
cd build_${variant}

# remove old output data
rm -rf out

./hodgkin_huxley_godunov ../settings_hodgkin_huxley_godunov.py
python ../scripts/check_results.py ../build_${variant}/out/
./dmdexample  

cd ..

