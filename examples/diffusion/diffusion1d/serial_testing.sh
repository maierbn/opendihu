#!/bin/bash

# this generates the following calls
# args: <n_nodes_per_fiber><scenario_name><solver_type> <preconditioner_type><petsc options>

#for i in $(seq 5); do
#  for j in 100 200 500 1000 1480 2000 5000 
#  do   
#    mpiexec -n 1 ./build_release/diffusion_1d settings_diffusion_IE_strong_scaling.py $j serial_diffusion_lu lu none 
#    mpiexec -n 1 ./build_release/diffusion_1d settings_diffusion_IE_strong_scaling.py $j serial_diffusion_cg cg none 
#    mpiexec -n 1 ./build_release/diffusion_1d settings_diffusion_IE_strong_scaling.py $j serial_diffusion_gmres gmres none 
#    mpiexec -n 1 ./build_release/diffusion_1d settings_diffusion_IE_strong_scaling.py $j serial_diffusion_gamg preonly gamg
    # boomeramg cannot be use as the only solver (not to be used with the preonly as solver) because the solution is wrong
#    mpiexec -n 1 ./build_release/diffusion_1d settings_diffusion_IE_strong_scaling.py $j serial_diffusion_boomeramg cg pchypre boomeramg #-ksp_type cg -pc_type hypre -pc_hypre_type boomeramg
#  done
#done

for i in $(seq 5); do
  for j in 100 200 500 1000 1480 2000 5000 
  do   
    mpiexec -n 1 ./build_release/diffusion_1d settings_diffusion_CN_strong_scaling.py $j serial_diffusion_lu lu none 
    mpiexec -n 1 ./build_release/diffusion_1d settings_diffusion_CN_strong_scaling.py $j serial_diffusion_cg cg none 
    mpiexec -n 1 ./build_release/diffusion_1d settings_diffusion_CN_strong_scaling.py $j serial_diffusion_gmres gmres none 
    mpiexec -n 1 ./build_release/diffusion_1d settings_diffusion_CN_strong_scaling.py $j serial_diffusion_gamg preonly gamg
    # boomeramg cannot be use as the only solver (not to be used with the preonly as solver) because the solution is wrong
    mpiexec -n 1 ./build_release/diffusion_1d settings_diffusion_CN_strong_scaling.py $j serial_diffusion_boomeramg cg pchypre boomeramg #-ksp_type cg -pc_type hypre -pc_hypre_type boomeramg
  done
done
