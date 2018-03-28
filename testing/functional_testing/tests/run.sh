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

. run_tests.sh
. postprocess.sh

done
