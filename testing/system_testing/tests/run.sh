#!/bin/bash

basedir=$(pwd)

# for all tests
for name in \
"laplace" \
"diffusion" \
"monodomain" \
"fibers" \
"multiple_fibers"
do

echo ""
echo $name
echo "=================="

# change directory to test directory
cd $basedir/$name

# compile
../../../../dependencies/scons/scons.py BUILD_TYPE=DEBUG
../../../../dependencies/scons/scons.py BUILD_TYPE=RELEASE

# run tests and postprocessing
. run_tests.sh
. postprocess.sh

# recompile documents
cd $basedir/../document
make

done
