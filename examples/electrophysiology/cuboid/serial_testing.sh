#!/bin/bash

#../run_tests.py  # this generates the following calls
#<n_processes_per_fiber> <n_fibers> <n_nodes_per_fiber> <scenario_name> <solver_type> <preconditioner_type>

mpiexec -n 1 ./cuboid ../cuboid_settings.py 1 1 100 solver_scaling gmres none
mpiexec -n 1 ./cuboid ../cuboid_settings.py 1 1 1000 solver_scaling gmres none
mpiexec -n 1 ./cuboid ../cuboid_settings.py 1 1 10000 solver_scaling gmres none
mpiexec -n 1 ./cuboid ../cuboid_settings.py 1 1 100000 solver_scaling gmres none

mpiexec -n 1 ./cuboid ../cuboid_settings.py 1 1 100 solver_scaling cg pchypre -pc_type hypre -pc_hypre_type boomeramg
mpiexec -n 1 ./cuboid ../cuboid_settings.py 1 1 1000 solver_scaling cg pchypre -pc_type hypre -pc_hypre_type boomeramg
mpiexec -n 1 ./cuboid ../cuboid_settings.py 1 1 10000 solver_scaling cg pchypre -pc_type hypre -pc_hypre_type boomeramg
mpiexec -n 1 ./cuboid ../cuboid_settings.py 1 1 100000 solver_scaling cg pchypre -pc_type hypre -pc_hypre_type boomeramg

mpiexec -n 1 ./cuboid ../cuboid_settings.py 1 1 100 solver_scaling richardson pchypre -pc_type hypre -pc_hypre_type boomeramg
mpiexec -n 1 ./cuboid ../cuboid_settings.py 1 1 1000 solver_scaling richardson pchypre -pc_type hypre -pc_hypre_type boomeramg
mpiexec -n 1 ./cuboid ../cuboid_settings.py 1 1 10000 solver_scaling richardson pchypre -pc_type hypre -pc_hypre_type boomeramg
mpiexec -n 1 ./cuboid ../cuboid_settings.py 1 1 100000 solver_scaling richardson pchypre -pc_type hypre -pc_hypre_type boomeramg

mpiexec -n 1 ./cuboid ../cuboid_settings.py 1 1 100 solver_scaling cg gamg
mpiexec -n 1 ./cuboid ../cuboid_settings.py 1 1 1000 solver_scaling cg gamg
mpiexec -n 1 ./cuboid ../cuboid_settings.py 1 1 10000 solver_scaling cg gamg
mpiexec -n 1 ./cuboid ../cuboid_settings.py 1 1 100000 solver_scaling cg gamg

mpiexec -n 1 ./cuboid ../cuboid_settings.py 1 1 100 solver_scaling preonly gamg
mpiexec -n 1 ./cuboid ../cuboid_settings.py 1 1 1000 solver_scaling preonly gamg
mpiexec -n 1 ./cuboid ../cuboid_settings.py 1 1 10000 solver_scaling preonly gamg
mpiexec -n 1 ./cuboid ../cuboid_settings.py 1 1 100000 solver_scaling preonly gamg
