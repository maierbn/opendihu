#!/bin/bash
if [ "$#" -ne 1 ]; then
    echo "usage: run_tests.sh <n processes>"
    return
fi
echo "running example $(pwd) with $1 processes"

workdir=$(pwd)
variant="debug"
variant="release"

# copy fibre results from fibers system test to this folder
cd input
. copy.sh
cd $workdir

mkdir -p build_${variant}
cd build_${variant}

if [ -d /data/scratch/maierbn ]; then  # if the output path exists (on neon)
  export output_path=/data/scratch/maierbn/multiple_fibres/out2
else
  export output_path=out2
fi
# remove old output data
rm -rf out

mkdir -p $output_path
ln -s $output_path out

export OMP_NUM_THREADS=1

# arguments: <n_processes_per_fiber>
mpirun -n $1 ./multiple_fibers ../multiple_fibers_settings.py 1

cd $workdir
