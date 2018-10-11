#!/bin/bash

#../run_tests.py  # this generates the following calls

# hardcoded calls for hazel hen

while true; do

# 1 node
if [ $nnodes = "1" ]; then
aprun -N 24 -n 24 ./cuboid ../cuboid_settings.py 24 1 2400 Strong_scaling >> out${nnodes}_Strong_scaling.txt 2>&1   # 1 node
fi

# 2 nodes
if [ $nnodes = "2" ]; then
aprun -N 24 -n 48 ./cuboid ../cuboid_settings.py 48 1 4800 Strong_scaling >> out${nnodes}_Strong_scaling.txt 2>&1   # 2 nodes
fi

# 13 nodes
if [ $nnodes = "13" ]; then
aprun -N 24 -n 144 ./cuboid ../cuboid_settings.py 144 1 14400 Strong_scaling >> out${nnodes}_Strong_scaling.txt 2>&1   # 6 nodes
aprun -N 24 -n 312 ./cuboid ../cuboid_settings.py 312 1 31200 Strong_scaling >> out${nnodes}_Strong_scaling.txt 2>&1   # 13 nodes
fi

# 32 nodes
if [ $nnodes = "32" ]; then
aprun -N 24 -n 768 ./cuboid ../cuboid_settings.py 768 1 76800 Strong_scaling >> out${nnodes}_Strong_scaling.txt 2>&1   # 32 nodes
fi

# 75 nodes
if [ $nnodes = "75" ]; then
aprun -N 24 -n 1800 ./cuboid ../cuboid_settings.py 1800 1 180000 Strong_scaling >> out${nnodes}_Strong_scaling.txt 2>&1   # 75 nodes
fi

# 178 nodes
if [ $nnodes = "178" ]; then
aprun -N 24 -n 4272 ./cuboid ../cuboid_settings.py 4272 1 427200 Strong_scaling >> out${nnodes}_Strong_scaling.txt 2>&1   # 178 nodes
fi

# 422 nodes
if [ $nnodes = "422" ]; then
aprun -N 24 -n 10128 ./cuboid ../cuboid_settings.py 10128 1 1012800 Strong_scaling >> out${nnodes}_Strong_scaling.txt 2>&1   # 422 nodes
fi

# 1000 nodes
if [ $nnodes = "1000" ]; then
aprun -N 24 -n 24000 ./cuboid ../cuboid_settings.py 24000 1 2400000 Strong_scaling >> out${nnodes}_Strong_scaling.txt 2>&1   # 1000 nodes
fi

done

cd $workdir
