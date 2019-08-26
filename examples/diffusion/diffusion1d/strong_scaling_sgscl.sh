#!/bin/bash
#SBATCH --job-name="strong scaling diffusion"
#SBATCH --nodes=6
#SBATCH --ntasks-per-node=8
#SBATCH --time=2:00:00

# strong scaling
# args: <n_nodes_per_fiber><scenario_name><solver_type> <preconditioner_type><petsc options>

for j in 100 200 500 1000 1480 2000 5000
do
  for i in $(seq 5); do
    for N in $(seq 6); do 
        for n in $(seq 8); do
           n_tot=$(python -c "print(($N-1)*8+$n)")
      	   srun mpirun.openmpi -N $N -n $n_tot ./build_release/diffusion_1d settings_diffusion_CN_strong_scaling.py $j strong_scaling_diffusion_lu lu none
	done
    done
  done
done
###########################################################
#for j in 100 200 500 1000 1480 2000 5000
#do
#  for i in $(seq 5); do
#    for N in $(seq 6); do 
#        for n in $(seq 8); do
#           n_tot=$(python -c "print(($N-1)*8+$n)")
#           srun mpirun.openmpi -N $N -n $n_tot ./build_release/diffusion_1d settings_diffusion_CN_strong_scaling.py $j strong_scaling_diffusion_cg cg none
#        done
#    done
#  done
#done
###########################################################
#for j in 100 200 500 1000 1480 2000 5000
#do
#  for i in $(seq 5); do
#    for N in $(seq 6); do 
#        for n in $(seq 8); do
#           n_tot=$(python -c "print(($N-1)*8+$n)")      
#           srun mpirun.openmpi -N $N -n $n_tot ./build_release/diffusion_1d settings_diffusion_CN_strong_scaling.py $j strong_scaling_diffusion_gmres gmres none
#	done
#    done
#  done
#done
###########################################################
#for j in 100 200 500 1000 1480 2000 5000
#do
#  for i in $(seq 5); do
#     for N in $(seq 6); do 
#        for n in $(seq 8); do
#           n_tot=$(python -c "print(($N-1)*8+$n)")      
#      	   srun mpirun.openmpi -N $N -n $n_tot ./build_release/diffusion_1d settings_diffusion_CN_strong_scaling.py $j strong_scaling_diffusion_gamg preonly gamg
#        done
#     done
#  done
#done

