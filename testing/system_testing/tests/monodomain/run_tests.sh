#!/bin/bash

echo "running example $(pwd)"

workdir=$(pwd)
variant="debug"
variant="release"

mkdir -p build_${variant}
cd build_${variant}

# remove old output data
rm -rf out

# arguments: <name> <cellml_file> <end_time>
./shorten_explicit ../settings.py shorten_opencor ../input/shorten_ocallaghan_davidson_soboleva_2007.c 1
./shorten_explicit ../settings.py shorten_opencmiss ../input/shorten_opencmiss.cpp 1

./hodgkin_huxley_explicit ../settings.py  hodgkin_huxley ../input/hodgkin_huxley_1952.c 1

cd $workdir

