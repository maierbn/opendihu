#!/bin/bash

#../run_tests.py  # this generates the following calls
#<n_processes_per_fiber> <n_fibers> <n_nodes_per_fiber> <scenario_name> <solver_type> <preconditioner_type>

mpiexec -n 1 ./cuboid ../cuboid_settings.py 1 1 10000 solver_scaling gmres none
mpiexec -n 4 ./cuboid ../cuboid_settings.py 4 1 10000 solver_scaling gmres none
mpiexec -n 8 ./cuboid ../cuboid_settings.py 8 1 10000 solver_scaling gmres none
mpiexec -n 18 ./cuboid ../cuboid_settings.py 18 1 10000 solver_scaling gmres none
mpiexec -n 30 ./cuboid ../cuboid_settings.py 30 1 10000 solver_scaling gmres none   
mpiexec -n 36 ./cuboid ../cuboid_settings.py 36 1 10000 solver_scaling gmres none     
mpiexec -n 42 ./cuboid ../cuboid_settings.py 42 1 10000 solver_scaling gmres none  
mpiexec -n 48 ./cuboid ../cuboid_settings.py 48 1 10000 solver_scaling gmres none   
mpiexec -n 54 ./cuboid ../cuboid_settings.py 54 1 10000 solver_scaling gmres none    
mpiexec -n 60 ./cuboid ../cuboid_settings.py 60 1 10000 solver_scaling gmres none 
mpiexec -n 66 ./cuboid ../cuboid_settings.py 66 1 10000 solver_scaling gmres none    
mpiexec -n 72 ./cuboid ../cuboid_settings.py 72 1 10000 solver_scaling gmres none 

mpiexec -n 1 ./cuboid ../cuboid_settings.py 1 1 10000 solver_scaling cg pchypre -pc_type hypre -pc_hypre_type boomeramg
mpiexec -n 4 ./cuboid ../cuboid_settings.py 4 1 10000 solver_scaling cg pchypre -pc_type hypre -pc_hypre_type boomeramg
mpiexec -n 8 ./cuboid ../cuboid_settings.py 8 1 10000 solver_scaling cg pchypre -pc_type hypre -pc_hypre_type boomeramg
mpiexec -n 18 ./cuboid ../cuboid_settings.py 18 1 10000 solver_scaling cg pchypre -pc_type hypre -pc_hypre_type boomeramg
mpiexec -n 30 ./cuboid ../cuboid_settings.py 30 1 10000 solver_scaling cg pchypre -pc_type hypre -pc_hypre_type boomeramg
mpiexec -n 36 ./cuboid ../cuboid_settings.py 36 1 10000 solver_scaling cg pchypre -pc_type hypre -pc_hypre_type boomeramg    
mpiexec -n 42 ./cuboid ../cuboid_settings.py 42 1 10000 solver_scaling cg pchypre -pc_type hypre -pc_hypre_type boomeramg
mpiexec -n 48 ./cuboid ../cuboid_settings.py 48 1 10000 solver_scaling cg pchypre -pc_type hypre -pc_hypre_type boomeramg
mpiexec -n 54 ./cuboid ../cuboid_settings.py 54 1 10000 solver_scaling cg pchypre -pc_type hypre -pc_hypre_type boomeramg
mpiexec -n 60 ./cuboid ../cuboid_settings.py 60 1 10000 solver_scaling cg pchypre -pc_type hypre -pc_hypre_type boomeramg
mpiexec -n 66 ./cuboid ../cuboid_settings.py 66 1 10000 solver_scaling cg pchypre -pc_type hypre -pc_hypre_type boomeramg
mpiexec -n 72 ./cuboid ../cuboid_settings.py 72 1 10000 solver_scaling cg pchypre -pc_type hypre -pc_hypre_type boomeramg

mpiexec -n 1 ./cuboid ../cuboid_settings.py 1 1 10000 solver_scaling richardson pchypre -pc_type hypre -pc_hypre_type boomeramg
mpiexec -n 4 ./cuboid ../cuboid_settings.py 4 1 10000 solver_scaling richardson pchypre -pc_type hypre -pc_hypre_type boomeramg
mpiexec -n 8 ./cuboid ../cuboid_settings.py 8 1 10000 solver_scaling richardson pchypre -pc_type hypre -pc_hypre_type boomeramg
mpiexec -n 18 ./cuboid ../cuboid_settings.py 18 1 10000 solver_scaling richardson pchypre -pc_type hypre -pc_hypre_type boomeramg
mpiexec -n 30 ./cuboid ../cuboid_settings.py 30 1 10000 solver_scaling richardson pchypre -pc_type hypre -pc_hypre_type boomeramg
mpiexec -n 36 ./cuboid ../cuboid_settings.py 36 1 10000 solver_scaling richardson pchypre -pc_type hypre -pc_hypre_type boomeramg
mpiexec -n 42 ./cuboid ../cuboid_settings.py 42 1 10000 solver_scaling richardson pchypre -pc_type hypre -pc_hypre_type boomeramg
mpiexec -n 48 ./cuboid ../cuboid_settings.py 48 1 10000 solver_scaling richardson pchypre -pc_type hypre -pc_hypre_type boomeramg
mpiexec -n 54 ./cuboid ../cuboid_settings.py 54 1 10000 solver_scaling richardson pchypre -pc_type hypre -pc_hypre_type boomeramg
mpiexec -n 60 ./cuboid ../cuboid_settings.py 60 1 10000 solver_scaling richardson pchypre -pc_type hypre -pc_hypre_type boomeramg
mpiexec -n 66 ./cuboid ../cuboid_settings.py 66 1 10000 solver_scaling richardson pchypre -pc_type hypre -pc_hypre_type boomeramg
mpiexec -n 72 ./cuboid ../cuboid_settings.py 72 1 10000 solver_scaling richardson pchypre -pc_type hypre -pc_hypre_type boomeramg

mpiexec -n 1 ./cuboid ../cuboid_settings.py 1 1 10000 solver_scaling cg gamg
mpiexec -n 4 ./cuboid ../cuboid_settings.py 4 1 10000 solver_scaling cg gamg
mpiexec -n 8 ./cuboid ../cuboid_settings.py 8 1 10000 solver_scaling cg gamg
mpiexec -n 18 ./cuboid ../cuboid_settings.py 18 1 10000 solver_scaling cg gamg
mpiexec -n 30 ./cuboid ../cuboid_settings.py 30 1 10000 solver_scaling cg gamg   
mpiexec -n 36 ./cuboid ../cuboid_settings.py 36 1 10000 solver_scaling cg gamg     
mpiexec -n 42 ./cuboid ../cuboid_settings.py 42 1 10000 solver_scaling cg gamg  
mpiexec -n 48 ./cuboid ../cuboid_settings.py 48 1 10000 solver_scaling cg gamg   
mpiexec -n 54 ./cuboid ../cuboid_settings.py 54 1 10000 solver_scaling cg gamg    
mpiexec -n 60 ./cuboid ../cuboid_settings.py 60 1 10000 solver_scaling cg gamg 
mpiexec -n 66 ./cuboid ../cuboid_settings.py 66 1 10000 solver_scaling cg gamg    
mpiexec -n 72 ./cuboid ../cuboid_settings.py 72 1 10000 solver_scaling cg gamg 

mpiexec -n 1 ./cuboid ../cuboid_settings.py 1 1 10000 solver_scaling preonly gamg
mpiexec -n 4 ./cuboid ../cuboid_settings.py 4 1 10000 solver_scaling preonly gamg
mpiexec -n 8 ./cuboid ../cuboid_settings.py 8 1 10000 solver_scaling preonly gamg
mpiexec -n 18 ./cuboid ../cuboid_settings.py 18 1 10000 solver_scaling preonly gamg
mpiexec -n 30 ./cuboid ../cuboid_settings.py 30 1 10000 solver_scaling preonly gamg   
mpiexec -n 36 ./cuboid ../cuboid_settings.py 36 1 10000 solver_scaling preonly gamg     
mpiexec -n 42 ./cuboid ../cuboid_settings.py 42 1 10000 solver_scaling preonly gamg  
mpiexec -n 48 ./cuboid ../cuboid_settings.py 48 1 10000 solver_scaling preonly gamg   
mpiexec -n 54 ./cuboid ../cuboid_settings.py 54 1 10000 solver_scaling preonly gamg    
mpiexec -n 60 ./cuboid ../cuboid_settings.py 60 1 10000 solver_scaling preonly gamg 
mpiexec -n 66 ./cuboid ../cuboid_settings.py 66 1 10000 solver_scaling preonly gamg    
mpiexec -n 72 ./cuboid ../cuboid_settings.py 72 1 10000 solver_scaling preonly gamg 