#!/bin/bash

workdir=$(pwd)

# for all tests
for name in \
"diffusion"
do

echo ""
echo $name
echo "=================="

# change directory to test directory
cd $workdir/$name

# run tests and postprocessing
. run_tests.sh
. postprocess.sh

# recomple documents
cd $workdir/../document
make all

done
