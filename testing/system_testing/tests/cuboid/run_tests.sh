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

# arguments: <n_processes_per_fiber>
#mpirun -n $1 ./cuboid ../cuboid_settings.py 1
./cuboid ../cuboid_settings.py 1

cd $workdir
