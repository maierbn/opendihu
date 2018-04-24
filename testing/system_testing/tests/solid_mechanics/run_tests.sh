#!/bin/bash

echo "running example $(pwd)"

workdir=$(pwd)
variant="debug"

mkdir -p build_${variant}
cd build_${variant}

# remove old output data
rm -rf out

# command arguments: <name> <number elements>

./mooney_rivlin_incompressible_mixed_condensation ../settings.py  mooney_rivlin_incompressible_mixed_condensation  10
./mooney_rivlin_incompressible_penalty            ../settings.py  mooney_rivlin_incompressible_penalty  10
./mooney_rivlin_incompressible_mixed              ../settings.py  mooney_rivlin_incompressible_mixed  10

cd $workdir
