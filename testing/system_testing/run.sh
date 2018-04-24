#!/bin/bash

workdir=$(pwd)

# run all tests, this includes compiling the documents
cd $workdir/tests
./run.sh
