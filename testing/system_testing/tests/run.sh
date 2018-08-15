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
echo "$name, starting at $(date +%T)"
echo "=============================="

# change directory to test directory
cd $basedir/$name


# compile
START=$(date +%s.%N)

../../../../dependencies/scons/scons.py BUILD_TYPE=DEBUG
../../../../dependencies/scons/scons.py BUILD_TYPE=RELEASE

END=$(date +%s.%N)
DIFF=$(echo "$END - $START" | bc)
echo "compilation took $(date -u -d @$DIFF +%T)"


# run tests
START=$(date +%s.%N)

. run_tests.sh

END=$(date +%s.%N)
DIFF=$(echo "$END - $START" | bc)
echo "running tests took $(date -u -d @$DIFF +%T)"


# run postprocessing
START=$(date +%s.%N)

. postprocess.sh

END=$(date +%s.%N)
DIFF=$(echo "$END - $START" | bc)
echo "postprocessing took $(date -u -d @$DIFF +%T)"


# recompile documents
cd $basedir/../document
make

done
