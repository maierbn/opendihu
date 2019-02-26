#!/bin/bash

basedir=$(pwd)

START_ALL=$(date +%s.%N)
rm log.txt

# for all tests
for name in \
"laplace" \
"diffusion" \
"monodomain" \
"fibers" \
"multiple_fibers" \
"monodomain_timestep_widths"
do

echo "" | tee -a $basedir/log.txt
echo "$name, starting at $(date +%T)" | tee -a $basedir/log.txt
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
echo "compilation took $(date -u -d @$DIFF +%T)" | tee -a $basedir/log.txt
echo ""


# run tests
START=$(date +%s.%N)

. run_tests.sh 1   # the last argument is the number of cores to use

END=$(date +%s.%N)
DIFF=$(python -c "print $END - $START")
echo ""
echo "running tests took $(date -u -d @$DIFF +%T)" | tee -a $basedir/log.txt
echo ""


# run postprocessing
START=$(date +%s.%N)

. postprocess.sh

END=$(date +%s.%N)
DIFF=$(python -c "print $END - $START")
echo ""
echo "postprocessing took $(date -u -d @$DIFF +%T)" | tee -a $basedir/log.txt
echo ""

done

# recompile documents
cd $basedir/../document
make

# output total duration
END_ALL=$(date +%s.%N)
DIFF_ALL=$(python -c "print $END_ALL - $START_ALL")
echo "" | tee -a log.txt
echo "total duration of system tests: $(date -u -d @$DIFF_ALL +%T)" | tee -a log.txt
echo ""
