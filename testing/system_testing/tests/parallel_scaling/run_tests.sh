#!/bin/bash

echo "running example $(pwd) with $1 processes"

workdir=$(pwd)
variant="debug"
variant="release"

cd $workdir

mkdir -p build_${variant}
cd build_${variant}

if [ -d /data/scratch/maierbn ]; then  # if the output path exists (on neon)
  export output_path=/data/scratch/maierbn/cuboid/out
else
  export output_path=out2
fi
# remove old output data
rm -rf out
mkdir -p $output_path
ln -s $output_path out

export OMP_NUM_THREADS=1

#../run_tests.py  # this generates the following calls

# hardcoded calls for hazel hen
aprun -N 24  -n 10000 ./cuboid ../cuboid_settings.py 100 100 1000 strong_scaling >> out_strong_scaling.txt 2>&1
aprun -N 24  -n 5623 ./cuboid ../cuboid_settings.py 56 100 1000 strong_scaling >> out_strong_scaling.txt 2>&1
aprun -N 24  -n 3162 ./cuboid ../cuboid_settings.py 32 100 1000 strong_scaling >> out_strong_scaling.txt 2>&1
aprun -N 24  -n 1778 ./cuboid ../cuboid_settings.py 18 100 1000 strong_scaling >> out_strong_scaling.txt 2>&1
aprun -N 24  -n 1000 ./cuboid ../cuboid_settings.py 10 100 1000 strong_scaling >> out_strong_scaling.txt 2>&1
aprun -N 24  -n 562 ./cuboid ../cuboid_settings.py 6 100 1000 strong_scaling >> out_strong_scaling.txt 2>&1
aprun -N 24  -n 316 ./cuboid ../cuboid_settings.py 3 100 1000 strong_scaling >> out_strong_scaling.txt 2>&1
aprun -N 24  -n 178 ./cuboid ../cuboid_settings.py 2 100 1000 strong_scaling >> out_strong_scaling.txt 2>&1
aprun -N 24  -n 100 ./cuboid ../cuboid_settings.py 1 100 1000 strong_scaling >> out_strong_scaling.txt 2>&1
echo strong scaling, multiple fibers per process
aprun -N 24  -n 22 ./cuboid ../cuboid_settings.py -5 100 1000 strong_scaling >> out_strong_scaling.txt 2>&1
aprun -N 24  -n 5 ./cuboid ../cuboid_settings.py -20 100 1000 strong_scaling >> out_strong_scaling.txt 2>&1
aprun -N 24  -n 1 ./cuboid ../cuboid_settings.py -100 100 1000 strong_scaling >> out_strong_scaling.txt 2>&1
echo weak scaling
aprun -N 24  -n 1 ./cuboid ../cuboid_settings.py 10 0 1000 weak_scaling >> out_weak_scaling.txt 2>&1
aprun -N 24  -n 3 ./cuboid ../cuboid_settings.py 10 0 1000 weak_scaling >> out_weak_scaling.txt 2>&1
aprun -N 24  -n 10 ./cuboid ../cuboid_settings.py 10 1 1000 weak_scaling >> out_weak_scaling.txt 2>&1
aprun -N 24  -n 32 ./cuboid ../cuboid_settings.py 10 3 1000 weak_scaling >> out_weak_scaling.txt 2>&1
aprun -N 24  -n 100 ./cuboid ../cuboid_settings.py 10 10 1000 weak_scaling >> out_weak_scaling.txt 2>&1
aprun -N 24  -n 316 ./cuboid ../cuboid_settings.py 10 32 1000 weak_scaling >> out_weak_scaling.txt 2>&1
aprun -N 24  -n 1000 ./cuboid ../cuboid_settings.py 10 100 1000 weak_scaling >> out_weak_scaling.txt 2>&1
aprun -N 24  -n 3162 ./cuboid ../cuboid_settings.py 10 316 1000 weak_scaling >> out_weak_scaling.txt 2>&1
aprun -N 24  -n 10000 ./cuboid ../cuboid_settings.py 10 1000 1000 weak_scaling >> out_weak_scaling.txt 2>&1

cd $workdir
