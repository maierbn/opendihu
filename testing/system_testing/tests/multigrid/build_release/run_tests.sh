#!/bin/bash

#../run_tests.py  # this generates the following calls
#<n_processes_per_fiber> <n_fibers> <n_nodes_per_fiber> <scenario_name> <solver_type> <preconditioner_type>

mpiexec -n 1 ./cuboid ../cuboid_settings.py 1 1 100 solver_scaling gmres none      | tee -a out_solver_scaling_1.txt 2>&1 # 1 nodes
mpiexec -n 1 ./cuboid ../cuboid_settings.py 1 1 100 solver_scaling preonly gamg   | tee -a out_solver_scaling_1.txt 2>&1 # 1 nodes
mpiexec -n 1 ./cuboid ../cuboid_settings.py 1 1 100 solver_scaling richardson pchypre      | tee -a out_solver_scaling_1.txt 2>&1 # 1 nodes

mpiexec -n 4 ./cuboid ../cuboid_settings.py 4 1 400 solver_scaling gmres none      | tee -a out_solver_scaling_1.txt 2>&1 # 1 nodes
mpiexec -n 4 ./cuboid ../cuboid_settings.py 4 1 400 solver_scaling preonly gamg    | tee -a out_solver_scaling_1.txt 2>&1 # 1 nodes
mpiexec -n 4 ./cuboid ../cuboid_settings.py 4 1 400 solver_scaling richardson pchypre      | tee -a out_solver_scaling_1.txt 2>&1 # 1 nodes

mpiexec -n 6 ./cuboid ../cuboid_settings.py 6 1 600 solver_scaling gmres none      | tee -a out_solver_scaling_1.txt 2>&1 # 1 nodes
mpiexec -n 6 ./cuboid ../cuboid_settings.py 6 1 600 solver_scaling preonly gamg    | tee -a out_solver_scaling_1.txt 2>&1 # 1 nodes
mpiexec -n 6 ./cuboid ../cuboid_settings.py 6 1 600 solver_scaling richardson pchypre      | tee -a out_solver_scaling_1.txt 2>&1 # 1 nodes

mpiexec -n 12 ./cuboid ../cuboid_settings.py 12 1 1200 solver_scaling gmres none      | tee -a out_solver_scaling_1.txt 2>&1 # 1 nodes
mpiexec -n 12 ./cuboid ../cuboid_settings.py 12 1 1200 solver_scaling preonly gamg    | tee -a out_solver_scaling_1.txt 2>&1 # 1 nodes
mpiexec -n 12 ./cuboid ../cuboid_settings.py 12 1 1200 solver_scaling richardson pchypre      | tee -a out_solver_scaling_1.txt 2>&1 # 1 nodes

#mpiexec -n 18 ./cuboid ../cuboid_settings.py 18 1 1800 solver_scaling gmres none      | tee -a out_solver_scaling_2.txt 2>&1 # 2 nodes
#mpiexec -n 18 ./cuboid ../cuboid_settings.py 18 1 1800 solver_scaling preonly gamg    | tee -a out_solver_scaling_2.txt 2>&1 # 2 nodes
#mpiexec -n 18 ./cuboid ../cuboid_settings.py 18 1 1800 solver_scaling richardson pchypre      | tee -a out_solver_scaling_2.txt 2>&1 # 2 nodes

#mpiexec -n 30 ./cuboid ../cuboid_settings.py 30 1 3000 solver_scaling gmres none      | tee -a out_solver_scaling_2.txt 2>&1 # 2 nodes
#mpiexec -n 30 ./cuboid ../cuboid_settings.py 30 1 3000 solver_scaling preonly gamg    | tee -a out_solver_scaling_2.txt 2>&1 # 2 nodes
#mpiexec -n 30 ./cuboid ../cuboid_settings.py 30 1 3000 solver_scaling richardson pchypre      | tee -a out_solver_scaling_2.txt 2>&1 # 2 nodes

#mpiexec -n 36 ./cuboid ../cuboid_settings.py 64 1 6400 solver_scaling gmres none      | tee -a out_solver_scaling_2.txt 2>&1 # 2 nodes
#mpiexec -n 36 ./cuboid ../cuboid_settings.py 64 1 6400 solver_scaling preonly gamg    | tee -a out_solver_scaling_2.txt 2>&1 # 2 nodes
#mpiexec -n 36 ./cuboid ../cuboid_settings.py 64 1 6400 solver_scaling richardson pchypre      | tee -a out_solver_scaling_2.txt 2>&1 # 2 nodes

#mpiexec -n 46 ./cuboid ../cuboid_settings.py 46 1 4600 solver_scaling gmres none      | tee -a out_solver_scaling_2.txt 2>&1 # 2 nodes
#mpiexec -n 46 ./cuboid ../cuboid_settings.py 46 1 4600 solver_scaling preonly gamg    | tee -a out_solver_scaling_2.txt 2>&1 # 2 nodes
#mpiexec -n 46 ./cuboid ../cuboid_settings.py 46 1 4600 solver_scaling richardson pchypre      | tee -a out_solver_scaling_2.txt 2>&1 # 2 nodes

#mpiexec -n 54 ./cuboid ../cuboid_settings.py 54 1 5400 solver_scaling gmres none      | tee -a out_solver_scaling_2.txt 2>&1 # 2 nodes
#mpiexec -n 54 ./cuboid ../cuboid_settings.py 54 1 5400 solver_scaling preonly gamg    | tee -a out_solver_scaling_2.txt 2>&1 # 2 nodes
#mpiexec -n 54 ./cuboid ../cuboid_settings.py 54 1 5400 solver_scaling richardson pchypre      | tee -a out_solver_scaling_2.txt 2>&1 # 2 nodes

#mpiexec -n 64 ./cuboid ../cuboid_settings.py 64 1 6400 solver_scaling gmres none      | tee -a out_solver_scaling_2.txt 2>&1 # 2 nodes
#mpiexec -n 64 ./cuboid ../cuboid_settings.py 64 1 6400 solver_scaling preonly gamg    | tee -a out_solver_scaling_2.txt 2>&1 # 2 nodes
#mpiexec -n 64 ./cuboid ../cuboid_settings.py 64 1 6400 solver_scaling richardson pchypre      | tee -a out_solver_scaling_2.txt 2>&1 # 2 nodes

#mpiexec -n 72 ./cuboid ../cuboid_settings.py 72 1 7200 solver_scaling gmres none      | tee -a out_solver_scaling_2.txt 2>&1 # 2 nodes
#mpiexec -n 72 ./cuboid ../cuboid_settings.py 72 1 7200 solver_scaling preonly gamg    | tee -a out_solver_scaling_2.txt 2>&1 # 2 nodes
#mpiexec -n 72 ./cuboid ../cuboid_settings.py 72 1 7200 solver_scaling richardson pchypre      | tee -a out_solver_scaling_2.txt 2>&1 # 2 nodes