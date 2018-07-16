#!/bin/bash

workdir=$(pwd)

# set path and pythonpath to include scripts
export PYTHONPATH=$PYTHONPATH:$(pwd)/../../scripts/
export PATH=$PATH:$(pwd)/../../scripts/

if [[ $(hostname) == "neon" ]]; then # on neon
module load tex/2012 
fi

# run all tests, this includes compiling the documents
cd $workdir/tests
./run.sh
