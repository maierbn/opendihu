#!/bin/bash

echo "running example $(pwd)"

workdir=$(pwd)
variant="debug"

mkdir -p build_${variant}
cd build_${variant}

# remove old output data
rm -rf out

# command arguments: <name> <number elements>
./mooney_rivlin_incompressible_penalty2d            ../settings_2d.py  0 mooney_rivlin_incompressible_penalty2d_numerical_jacobian
./mooney_rivlin_incompressible_penalty2d            ../settings_2d.py  1 mooney_rivlin_incompressible_penalty2d_analytical_jacobian
#./mooney_rivlin_incompressible_mixed_condensation ../settings.py  mooney_rivlin_incompressible_mixed_condensation  10
#./mooney_rivlin_incompressible_mixed              ../settings.py  mooney_rivlin_incompressible_mixed  10

cd $workdir
