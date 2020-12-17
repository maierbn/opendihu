# this file is created from the jobscript on hawk

#!/bin/bash
#PBS -N with_electrophysiology
#PBS -l select=2:node_type=rome:mpiprocs=64
#PBS -l walltime=01:00:00
#PBS -j oe
#PBS -o stdout

# load modules  
#. ~/load_modules_hawk.sh
#unset PYTHONPATH

# Change to the directory that the job was submitted from
export PBS_O_WORKDIR=$(pwd)    # comment this in for local execution without job
cd $PBS_O_WORKDIR

export OMP_NUM_THREADS=1
export PRGENV=gnu
#export OPENDIHU_HOME=/lustre/cray/ws9/2/ws/icbbnmai-opendihu1/opendihu-hawk-${PRGENV}

export EXAMPLE_HOME=${OPENDIHU_HOME}/examples/electrophysiology/fibers/fibers_contraction/with_tendons_precice/multiple_tendons_with_electrophysiology
export INPUT_DIR=${OPENDIHU_HOME}/examples/electrophysiology/input

# specify partitioning for muscle simulation [x,y,z,f], f is number of fibers in one coordinate direction
partitioning="[2,2,15,31]"
# values from weak scaling:
#"[2,2,1,7]" \
#"[3,3,2,13]" \
#"[4,4,4,25]" \
#"[6,6,4,37]" \
#"[7,8,8,67]" \
#"[12,12,8,109]" \
#"[15,15,16,187]" \
#"[22,22,16,277]" \
#"[24,24,32,427]" \
#"[29,29,32,523]" \

export X=`python -c "print($partitioning[0])"`
export Y=`python -c "print($partitioning[1])"`
export Z=`python -c "print($partitioning[2])"`
export F=`python -c "print($partitioning[3])"`
export XYZ=`python -c "print($X*$Y*$Z)"`

export OFFSET_TENDON_BOTTOM=`python -c "print($X*$Y*$Z)"`
export OFFSET_TENDON_TOP_A=`python -c "print($X*$Y*$Z+4)"`
export OFFSET_TENDON_TOP_B=`python -c "print($X*$Y*$Z+8)"`

echo "muscle: (X,Y,Z)=($X,$Y,$Z), n processes: $XYZ, n fibers: ${F}x${F}, n tasks: $PBS_TASKNUM, n nodes: $PBS_NODENUM"

# The ${PBS_NODEFILE} contains the list of hosts for the MPI ranks. It can be altered prior to HPE mpirun

# -------------------
# muscle

mpirun -np $XYZ \
${EXAMPLE_HOME}/muscle_electrophysiology_precice \
  ${EXAMPLE_HOME}/settings_muscle.py hawk_medium.py \
    --n_subdomains $X $Y $Z \
    --fiber_file   ${INPUT_DIR}/left_biceps_brachii_${F}x${F}fibers.bin \
  | tee ${PBS_O_WORKDIR}/stdout/${PBS_JOBID}.muscle.txt  &


# domain decomposition for tendons
export X=2
export Y=1
export Z=2
export XYZ=`python -c "print($X*$Y*$Z)"`

# -------------------
# tendon_bottom
mpirun -np $XYZ \
${EXAMPLE_HOME}/tendon_linear_precice_dynamic \
  ${EXAMPLE_HOME}/settings_tendon_bottom.py \
    --n_subdomains $X $Y $Z \
  | tee ${PBS_O_WORKDIR}/stdout/${PBS_JOBID}.tendon_bottom.txt  &

# -------------------
# tendon_top_a
mpirun -np $XYZ \
${EXAMPLE_HOME}/tendon_linear_precice_dynamic \
  ${EXAMPLE_HOME}/settings_tendon_top_a.py \
    --n_subdomains $X $Y $Z \
  | tee ${PBS_O_WORKDIR}/stdout/${PBS_JOBID}.tendon_top_a.txt  &

# -------------------
# tendon_top_b
mpirun -np $XYZ \
${EXAMPLE_HOME}/tendon_linear_precice_dynamic \
  ${EXAMPLE_HOME}/settings_tendon_top_b.py \
    --n_subdomains $X $Y $Z \
  | tee ${PBS_O_WORKDIR}/stdout/${PBS_JOBID}.tendon_top_b.txt 


echo "end"
