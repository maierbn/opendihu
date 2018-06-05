#!/bin/bash

workdir=$(pwd)

# set path and pythonpath to include scripts
export PYTHONPATH=$PYTHONPATH:$(pwd)/../../scripts/
export PATH=$PATH:$(pwd)/../../scripts/

# run all tests, this includes compiling the documents
cd $workdir/tests
./run.sh
