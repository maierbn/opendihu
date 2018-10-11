#!/bin/bash

#../run_tests.py  # this generates the following calls

# hardcoded calls for hazel hen

while true; do

# 1 node
if [ $nnodes = "1" ]; then
aprun -N 22  -n 22 ./cuboid ../cuboid_settings.py -5 100 1000 strong_scaling >> out${nnodes}_strong_scaling.txt 2>&1
aprun -N 5  -n 5 ./cuboid ../cuboid_settings.py -20 100 1000 strong_scaling >> out${nnodes}_strong_scaling.txt 2>&1
aprun -N 1  -n 1 ./cuboid ../cuboid_settings.py -100 100 1000 strong_scaling >> out${nnodes}_strong_scaling.txt 2>&1
aprun -N 3  -n 3 ./cuboid ../cuboid_settings.py 10 0 1000 weak_scaling >> out${nnodes}_weak_scaling.txt 2>&1
aprun -N 1  -n 1 ./cuboid ../cuboid_settings.py 10 0 1000 weak_scaling >> out${nnodes}_weak_scaling.txt 2>&1
aprun -N 10  -n 10 ./cuboid ../cuboid_settings.py 10 1 1000 weak_scaling >> out${nnodes}_weak_scaling.txt 2>&1
fi

# 14 nodes
if [ $nnodes = "14" ]; then
aprun -N 24  -n 32 ./cuboid ../cuboid_settings.py 10 3 1000 weak_scaling >> out${nnodes}_weak_scaling.txt 2>&1
aprun -N 24  -n 316 ./cuboid ../cuboid_settings.py 3 100 1000 strong_scaling >> out${nnodes}_strong_scaling.txt 2>&1
aprun -N 24  -n 178 ./cuboid ../cuboid_settings.py 2 100 1000 strong_scaling >> out${nnodes}_strong_scaling.txt 2>&1
aprun -N 24  -n 100 ./cuboid ../cuboid_settings.py 1 100 1000 strong_scaling >> out${nnodes}_strong_scaling.txt 2>&1
aprun -N 24  -n 100 ./cuboid ../cuboid_settings.py 10 10 1000 weak_scaling >> out${nnodes}_weak_scaling.txt 2>&1
aprun -N 24  -n 316 ./cuboid ../cuboid_settings.py 10 32 1000 weak_scaling >> out${nnodes}_weak_scaling.txt 2>&1
fi

# 24 nodes
if [ $nnodes = "24" ]; then
aprun -N 24  -n 562 ./cuboid ../cuboid_settings.py 6 100 1000 strong_scaling >> out${nnodes}_strong_scaling.txt 2>&1
fi

# 42 nodes
if [ $nnodes = "42" ]; then
aprun -N 24  -n 1000 ./cuboid ../cuboid_settings.py 10 100 1000 weak_scaling >> out${nnodes}_weak_scaling.txt 2>&1
aprun -N 24  -n 1000 ./cuboid ../cuboid_settings.py 10 100 1000 strong_scaling >> out${nnodes}_strong_scaling.txt 2>&1
fi

# 75 nodes
if [ $nnodes = "75" ]; then
aprun -N 24  -n 1778 ./cuboid ../cuboid_settings.py 18 100 1000 strong_scaling >> out${nnodes}_strong_scaling.txt 2>&1
fi

# 321 nodes
if [ $nnodes = "321" ]; then
aprun -N 24  -n 3162 ./cuboid ../cuboid_settings.py 10 316 1000 weak_scaling >> out${nnodes}_weak_scaling.txt 2>&1
aprun -N 24  -n 3162 ./cuboid ../cuboid_settings.py 32 100 1000 strong_scaling >> out${nnodes}_strong_scaling.txt 2>&1
fi

# 235 nodes
if [ $nnodes = "235" ]; then
aprun -N 24  -n 5623 ./cuboid ../cuboid_settings.py 56 100 1000 strong_scaling >> out${nnodes}_strong_scaling.txt 2>&1
fi

# 417 nodes
if [ $nnodes = "417" ]; then
aprun -N 24  -n 10000 ./cuboid ../cuboid_settings.py 100 100 1000 strong_scaling >> out${nnodes}_strong_scaling.txt 2>&1
aprun -N 24  -n 10000 ./cuboid ../cuboid_settings.py 10 1000 1000 weak_scaling >> out${nnodes}_weak_scaling.txt 2>&1
fi

done

cd $workdir
