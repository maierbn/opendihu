#!/bin/bash

basedir=$(pwd)

# for all tests
for name in \
"diffusion" \
"solid_mechanics"
do

echo ""
echo $name
echo "=================="

# change directory to test directory
cd $basedir/$name

# compile
../../../../dependencies/scons/scons.py BUILD_TYPE=DEBUG

# run tests and postprocessing
. run_tests.sh
. postprocess.sh

# recompile documents
cd $basedir/../document
make

done
