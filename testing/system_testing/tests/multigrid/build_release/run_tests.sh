#!/bin/bash

#../run_tests.py  # this generates the following calls
#<n_processes_per_fiber> <n_fibers> <n_nodes_per_fiber> <scenario_name> <solver_type> <preconditioner_type>

mpiexec -n 1 ./cuboid ../cuboid_settings.py 1 1 100 solver_scaling gmres none      | tee -a out_solver_scaling_1.txt 2>&1 # 1 nodes
mpiexec -n 1 ./cuboid ../cuboid_settings.py 1 1 100 solver_scaling preonly gamg   | tee -a out_solver_scaling_1.txt 2>&1 # 1 nodes
mpiexec -n 1 ./cuboid ../cuboid_settings.py 1 1 100 solver_scaling richardson pchypre      | tee -a out_solver_scaling_1.txt 2>&1 # 1 nodes

mpiexec -n 4 ./cuboid ../cuboid_settings.py 4 1 400 solver_scaling gmres none      | tee -a out_solver_scaling_1.txt 2>&1 # 1 nodes
mpiexec -n 4 ./cuboid ../cuboid_settings.py 4 1 400 solver_scaling preonly gamg    | tee -a out_solver_scaling_1.txt 2>&1 # 1 nodes
mpiexec -n 4 ./cuboid ../cuboid_settings.py 4 1 400 solver_scaling richardson pchypre      | tee -a out_solver_scaling_1.txt 2>&1 # 1 nodes

mpiexec -n 8 ./cuboid ../cuboid_settings.py 8 1 800 solver_scaling gmres none      | tee -a out_solver_scaling_1.txt 2>&1 # 1 nodes
mpiexec -n 8 ./cuboid ../cuboid_settings.py 8 1 800 solver_scaling preonly gamg    | tee -a out_solver_scaling_1.txt 2>&1 # 1 nodes
mpiexec -n 8 ./cuboid ../cuboid_settings.py 8 1 800 solver_scaling richardson pchypre      | tee -a out_solver_scaling_1.txt 2>&1 # 1 nodes

#mpiexec -n 16 ./cuboid ../cuboid_settings.py 16 1 1600 solver_scaling gmres none      | tee -a out_solver_scaling_2.txt 2>&1 # 2 nodes
#mpiexec -n 16 ./cuboid ../cuboid_settings.py 16 1 1600 solver_scaling preonly gamg    | tee -a out_solver_scaling_2.txt 2>&1 # 2 nodes
#mpiexec -n 16 ./cuboid ../cuboid_settings.py 16 1 1600 solver_scaling richardson pchypre      | tee -a out_solver_scaling_2.txt 2>&1 # 2 nodes

#mpiexec -n 32 ./cuboid ../cuboid_settings.py 32 1 3200 solver_scaling gmres none      | tee -a out_solver_scaling_2.txt 2>&1 # 2 nodes
#mpiexec -n 32 ./cuboid ../cuboid_settings.py 32 1 3200 solver_scaling preonly gamg    | tee -a out_solver_scaling_2.txt 2>&1 # 2 nodes
#mpiexec -n 32 ./cuboid ../cuboid_settings.py 32 1 3200 solver_scaling richardson pchypre      | tee -a out_solver_scaling_2.txt 2>&1 # 2 nodes

#mpiexec -n 64 ./cuboid ../cuboid_settings.py 64 1 6400 solver_scaling gmres none      | tee -a out_solver_scaling_2.txt 2>&1 # 2 nodes
#mpiexec -n 64 ./cuboid ../cuboid_settings.py 64 1 6400 solver_scaling preonly gamg    | tee -a out_solver_scaling_2.txt 2>&1 # 2 nodes
#mpiexec -n 64 ./cuboid ../cuboid_settings.py 64 1 6400 solver_scaling richardson pchypre      | tee -a out_solver_scaling_2.txt 2>&1 # 2 nodes