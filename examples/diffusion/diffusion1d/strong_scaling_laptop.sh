# strong scaling
# args: <n_nodes_per_fiber><scenario_name><solver_type> <preconditioner_type><petsc options>

for j in 100 200 500 1000 1480
do
  for i in $(seq 5); do
    for n in $(seq 4); do      
      mpiexec -n $n ./build_release/diffusion_1d settings_diffusion_implicit_strong_scaling.py $j strong_scaling_diffusion_lu lu none
    done
  done
done
###########################################################
for j in 100 200 500 1000 1480
do
  for i in $(seq 5); do
    for n in $(seq 4); do      
      mpiexec -n $n ./build_release/diffusion_1d settings_diffusion_implicit_strong_scaling.py $j strong_scaling_diffusion_cg cg none
    done
  done
done
###########################################################
#for j in 100 200 500 1000 1480
#do
#  for i in $(seq 5); do
#    for n in $(seq 4); do      
#      mpiexec -n $n ./build_release/diffusion_1d settings_diffusion_implicit_strong_scaling.py $j strong_scaling_diffusion_gmres gmres none
#    done
#  done
#done
###########################################################
for j in 100 200 500 1000 1480
do
  for i in $(seq 5); do
    for n in $(seq 4); do      
      mpiexec -n $n ./build_release/diffusion_1d settings_diffusion_implicit_strong_scaling.py $j strong_scaling_diffusion_gamg preonly gamg
    done
  done
done
###########################################################
#for j in 100 200 500 1000 1480
#do
#  for i in $(seq 5); do
#    for n in $(seq 4); do      
#      mpiexec -n $n ./build_release/diffusion_1d settings_diffusion_implicit_strong_scaling.py $j strong_scaling_diffusion_boomeramg cg pchypre boomeramg
#    done
#  done
#done
###########################################################
