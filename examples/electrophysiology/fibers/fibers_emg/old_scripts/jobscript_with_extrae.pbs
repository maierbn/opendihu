#!/bin/sh
#PBS -l walltime=00:20:00
#PBS -l select=1:node_type=rome:mpiprocs=16
#PBS -j oe

# submit from build_release

# Change to directory the job was submitted from
cd $PBS_O_WORKDIR

# Load required modules
module load openmpi
module load petsc/3.12.2-int32-shared
module load extrae

module list

export TRACE_NAME=fast_fibers_emg.8.small.prv
export TMPDIR=$PBS_O_WORKDIR/tmp # besser einen ablsotuen Pfad in einen Workspace zu nehmen

# copy extrae script and config file to build_release
cp ../old_scripts/extrae_trace.sh .
chmod +x extrae_trace.sh
cp ../old_scripts/extrae.xml .

echo cwd: `pwd`
ls
echo now start program with 8 ranks

# Run
mpirun -np 8 ./extrae_trace.sh ./fast_fibers_emg ../settings_fibers_emg.py ramp_emg.py --dt_0D 1e-3 --dt_1D 2e-3 --dt_splitting 2e-3 --dt_3D 2e-3 --end_time 8e-3 --diffusion_solver_reltol 0 --diffusion_solver_maxit 5 --emg_solver_reltol 0 --emg_solver_maxit 5 --potential_flow_solver_reltol 0 --potential_flow_solver_maxit 5 

# open result with paraver
