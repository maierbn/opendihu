#!/bin/bash

#../run_tests.py  # this generates the following calls
#<n_processes_per_fiber> <n_fibers> <n_nodes_per_fiber> <scenario_name> <solver_type> <preconditioner_type>

mpiexec -n 1 ./diffusion_1d_mg settings_scaling.py 100 cg pchypre cg_boomeramg -ksp_type cg -pc_type hypre -pc_hypre_type boomeramg
mpiexec -n 1 ./diffusion_1d_mg settings_scaling.py 1000 cg pchypre cg_boomeramg -ksp_type cg -pc_type hypre -pc_hypre_type boomeramg
mpiexec -n 1 ./diffusion_1d_mg settings_scaling.py 10000 cg pchypre cg_boomeramg -ksp_type cg -pc_type hypre -pc_hypre_type boomeramg
mpiexec -n 1 ./diffusion_1d_mg settings_scaling.py 100000 cg pchypre cg_boomeramg -ksp_type cg -pc_type hypre -pc_hypre_type boomeramg

###########################################################

mpiexec -n 1 ./diffusion_1d_mg settings_scaling.py 100 richardson pchypre boomeramg -ksp_type richardson -pc_type hypre -pc_hypre_type boomeramg
mpiexec -n 1 ./diffusion_1d_mg settings_scaling.py 1000 richardson pchypre boomeramg -ksp_type richardson -pc_type hypre -pc_hypre_type boomeramg
mpiexec -n 1 ./diffusion_1d_mg settings_scaling.py 10000 richardson pchypre boomeramg -ksp_type richardson -pc_type hypre -pc_hypre_type boomeramg
mpiexec -n 1 ./diffusion_1d_mg settings_scaling.py 100000 richardson pchypre boomeramg -ksp_type richardson -pc_type hypre -pc_hypre_type boomeramg

############################################################

mpiexec -n 1 ./diffusion_1d_mg settings_scaling.py 100 cg gamg cg_gamg -ksp_type cg -pc_type gamg
mpiexec -n 1 ./diffusion_1d_mg settings_scaling.py 1000 cg gamg cg_gamg -ksp_type cg -pc_type gamg
mpiexec -n 1 ./diffusion_1d_mg settings_scaling.py 10000 cg gamg cg_gamg -ksp_type cg -pc_type gamg
mpiexec -n 1 ./diffusion_1d_mg settings_scaling.py 100000 cg gamg cg_gamg -ksp_type cg -pc_type gamg

#############################################################

mpiexec -n 1 ./diffusion_1d_mg settings_scaling.py 100 preonly gamg gamg -ksp_type preonly -pc_type gamg
mpiexec -n 1 ./diffusion_1d_mg settings_scaling.py 1000 preonly gamg gamg -ksp_type preonly -pc_type gamg
mpiexec -n 1 ./diffusion_1d_mg settings_scaling.py 10000 preonly gamg gamg -ksp_type preonly -pc_type gamg
mpiexec -n 1 ./diffusion_1d_mg settings_scaling.py 100000 preonly gamg gamg -ksp_type preonly -pc_type gamg

#############################################################

mpiexec -n 1 ./diffusion_1d_mg settings_scaling.py 100 gmres none gmres -ksp_type gmres -pc_type none
mpiexec -n 1 ./diffusion_1d_mg settings_scaling.py 1000 gmres none gmres -ksp_type gmres -pc_type none
mpiexec -n 1 ./diffusion_1d_mg settings_scaling.py 10000 gmres none gmres -ksp_type gmres -pc_type none
mpiexec -n 1 ./diffusion_1d_mg settings_scaling.py 100000 gmres none gmres -ksp_type gmres -pc_type none
