#!/bin/bash

basedir=$(pwd)

START_ALL=$(date +%s.%N)

# for all tests
for name in \
"laplace" \
"diffusion" \
"monodomain" \
"fibers" \
"multiple_fibers" \
"monodomain_timestep_widths"
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
DIFF=$(python -c "print $END - $START")
echo ""
echo "compilation took $(date -u -d @$DIFF +%T)"
echo ""


# run tests
START=$(date +%s.%N)

. run_tests.sh 1   # the last argument is the number of cores to use

END=$(date +%s.%N)
DIFF=$(python -c "print $END - $START")
echo ""
echo "running tests took $(date -u -d @$DIFF +%T)"
echo ""


# run postprocessing
START=$(date +%s.%N)

. postprocess.sh

END=$(date +%s.%N)
DIFF=$(python -c "print $END - $START")
echo ""
echo "postprocessing took $(date -u -d @$DIFF +%T)"
echo ""


# recompile documents
cd $basedir/../document
make

done

# output total duration
END_ALL=$(date +%s.%N)
DIFF_ALL=$(python -c "print $END_ALL - $START_ALL")
echo ""
echo "total duration of system tests: $(date -u -d @$DIFF_ALL +%T)"
echo ""
