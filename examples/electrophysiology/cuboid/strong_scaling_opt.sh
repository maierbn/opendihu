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

mpiexec -n 1 ./cuboid ../cuboid_settings.py 1 1 10000 solver_scaling cg pchypre -pc_type hypre -pc_hypre_type boomeramg -pc_hypre_boomeramg_coarsen_type HMIS -pc_hypre_boomeramg_agg_nl 2 -pc_hypre_boomeramg_strong_threshold 0.3
mpiexec -n 4 ./cuboid ../cuboid_settings.py 4 1 10000 solver_scaling cg pchypre -pc_type hypre -pc_hypre_type boomeramg -pc_hypre_boomeramg_coarsen_type HMIS -pc_hypre_boomeramg_agg_nl 2 -pc_hypre_boomeramg_strong_threshold 0.3
mpiexec -n 8 ./cuboid ../cuboid_settings.py 8 1 10000 solver_scaling cg pchypre -pc_type hypre -pc_hypre_type boomeramg -pc_hypre_boomeramg_coarsen_type HMIS -pc_hypre_boomeramg_agg_nl 2 -pc_hypre_boomeramg_strong_threshold 0.3
mpiexec -n 18 ./cuboid ../cuboid_settings.py 18 1 10000 solver_scaling cg pchypre -pc_type hypre -pc_hypre_type boomeramg -pc_hypre_boomeramg_coarsen_type HMIS -pc_hypre_boomeramg_agg_nl 2 -pc_hypre_boomeramg_strong_threshold 0.3
mpiexec -n 30 ./cuboid ../cuboid_settings.py 30 1 10000 solver_scaling cg pchypre -pc_type hypre -pc_hypre_type boomeramg -pc_hypre_boomeramg_coarsen_type HMIS -pc_hypre_boomeramg_agg_nl 2 -pc_hypre_boomeramg_strong_threshold 0.3
mpiexec -n 36 ./cuboid ../cuboid_settings.py 36 1 10000 solver_scaling cg pchypre -pc_type hypre -pc_hypre_type boomeramg -pc_hypre_boomeramg_coarsen_type HMIS -pc_hypre_boomeramg_agg_nl 2 -pc_hypre_boomeramg_strong_threshold 0.3
mpiexec -n 42 ./cuboid ../cuboid_settings.py 42 1 10000 solver_scaling cg pchypre -pc_type hypre -pc_hypre_type boomeramg -pc_hypre_boomeramg_coarsen_type HMIS -pc_hypre_boomeramg_agg_nl 2 -pc_hypre_boomeramg_strong_threshold 0.3
mpiexec -n 48 ./cuboid ../cuboid_settings.py 48 1 10000 solver_scaling cg pchypre -pc_type hypre -pc_hypre_type boomeramg -pc_hypre_boomeramg_coarsen_type HMIS -pc_hypre_boomeramg_agg_nl 2 -pc_hypre_boomeramg_strong_threshold 0.3
mpiexec -n 54 ./cuboid ../cuboid_settings.py 54 1 10000 solver_scaling cg pchypre -pc_type hypre -pc_hypre_type boomeramg -pc_hypre_boomeramg_coarsen_type HMIS -pc_hypre_boomeramg_agg_nl 2 -pc_hypre_boomeramg_strong_threshold 0.3
mpiexec -n 60 ./cuboid ../cuboid_settings.py 60 1 10000 solver_scaling cg pchypre -pc_type hypre -pc_hypre_type boomeramg -pc_hypre_boomeramg_coarsen_type HMIS -pc_hypre_boomeramg_agg_nl 2 -pc_hypre_boomeramg_strong_threshold 0.3
mpiexec -n 66 ./cuboid ../cuboid_settings.py 66 1 10000 solver_scaling cg pchypre -pc_type hypre -pc_hypre_type boomeramg -pc_hypre_boomeramg_coarsen_type HMIS -pc_hypre_boomeramg_agg_nl 2 -pc_hypre_boomeramg_strong_threshold 0.3
mpiexec -n 72 ./cuboid ../cuboid_settings.py 72 1 10000 solver_scaling cg pchypre -pc_type hypre -pc_hypre_type boomeramg -pc_hypre_boomeramg_coarsen_type HMIS -pc_hypre_boomeramg_agg_nl 2 -pc_hypre_boomeramg_strong_threshold 0.3

mpiexec -n 1 ./cuboid ../cuboid_settings.py 1 1 10000 solver_scaling richardson pchypre -pc_type hypre -pc_hypre_type boomeramg -pc_hypre_boomeramg_coarsen_type Falgout -pc_hypre_boomeramg_agg_nl 0 -pc_hypre_boomeramg_strong_threshold 0.3
mpiexec -n 4 ./cuboid ../cuboid_settings.py 4 1 10000 solver_scaling richardson pchypre -pc_type hypre -pc_hypre_type boomeramg -pc_hypre_boomeramg_coarsen_type Falgout -pc_hypre_boomeramg_agg_nl 0 -pc_hypre_boomeramg_strong_threshold 0.3
mpiexec -n 8 ./cuboid ../cuboid_settings.py 8 1 10000 solver_scaling richardson pchypre -pc_type hypre -pc_hypre_type boomeramg -pc_hypre_boomeramg_coarsen_type Falgout -pc_hypre_boomeramg_agg_nl 0 -pc_hypre_boomeramg_strong_threshold 0.3
mpiexec -n 18 ./cuboid ../cuboid_settings.py 18 1 10000 solver_scaling richardson pchypre -pc_type hypre -pc_hypre_type boomeramg -pc_hypre_boomeramg_coarsen_type Falgout -pc_hypre_boomeramg_agg_nl 0 -pc_hypre_boomeramg_strong_threshold 0.3
mpiexec -n 30 ./cuboid ../cuboid_settings.py 30 1 10000 solver_scaling richardson pchypre -pc_type hypre -pc_hypre_type boomeramg -pc_hypre_boomeramg_coarsen_type Falgout -pc_hypre_boomeramg_agg_nl 0 -pc_hypre_boomeramg_strong_threshold 0.3
mpiexec -n 36 ./cuboid ../cuboid_settings.py 36 1 10000 solver_scaling richardson pchypre -pc_type hypre -pc_hypre_type boomeramg -pc_hypre_boomeramg_coarsen_type Falgout -pc_hypre_boomeramg_agg_nl 0 -pc_hypre_boomeramg_strong_threshold 0.3
mpiexec -n 42 ./cuboid ../cuboid_settings.py 42 1 10000 solver_scaling richardson pchypre -pc_type hypre -pc_hypre_type boomeramg -pc_hypre_boomeramg_coarsen_type Falgout -pc_hypre_boomeramg_agg_nl 0 -pc_hypre_boomeramg_strong_threshold 0.3
mpiexec -n 48 ./cuboid ../cuboid_settings.py 48 1 10000 solver_scaling richardson pchypre -pc_type hypre -pc_hypre_type boomeramg -pc_hypre_boomeramg_coarsen_type Falgout -pc_hypre_boomeramg_agg_nl 0 -pc_hypre_boomeramg_strong_threshold 0.3
mpiexec -n 54 ./cuboid ../cuboid_settings.py 54 1 10000 solver_scaling richardson pchypre -pc_type hypre -pc_hypre_type boomeramg -pc_hypre_boomeramg_coarsen_type Falgout -pc_hypre_boomeramg_agg_nl 0 -pc_hypre_boomeramg_strong_threshold 0.3
mpiexec -n 60 ./cuboid ../cuboid_settings.py 60 1 10000 solver_scaling richardson pchypre -pc_type hypre -pc_hypre_type boomeramg -pc_hypre_boomeramg_coarsen_type Falgout -pc_hypre_boomeramg_agg_nl 0 -pc_hypre_boomeramg_strong_threshold 0.3
mpiexec -n 66 ./cuboid ../cuboid_settings.py 66 1 10000 solver_scaling richardson pchypre -pc_type hypre -pc_hypre_type boomeramg -pc_hypre_boomeramg_coarsen_type Falgout -pc_hypre_boomeramg_agg_nl 0 -pc_hypre_boomeramg_strong_threshold 0.3
mpiexec -n 72 ./cuboid ../cuboid_settings.py 72 1 10000 solver_scaling richardson pchypre -pc_type hypre -pc_hypre_type boomeramg -pc_hypre_boomeramg_coarsen_type Falgout -pc_hypre_boomeramg_agg_nl 0 -pc_hypre_boomeramg_strong_threshold 0.3

mpiexec -n 1 ./cuboid ../cuboid_settings.py 1 1 10000 solver_scaling cg gamg -pc_gamg_square_graph 4 -pc_gamg_threshold 0.06
mpiexec -n 4 ./cuboid ../cuboid_settings.py 4 1 10000 solver_scaling cg gamg -pc_gamg_square_graph 4 -pc_gamg_threshold 0.06
mpiexec -n 8 ./cuboid ../cuboid_settings.py 8 1 10000 solver_scaling cg gamg -pc_gamg_square_graph 4 -pc_gamg_threshold 0.06
mpiexec -n 18 ./cuboid ../cuboid_settings.py 18 1 10000 solver_scaling cg gamg -pc_gamg_square_graph 4 -pc_gamg_threshold 0.06
mpiexec -n 30 ./cuboid ../cuboid_settings.py 30 1 10000 solver_scaling cg gamg -pc_gamg_square_graph 4 -pc_gamg_threshold 0.06
mpiexec -n 36 ./cuboid ../cuboid_settings.py 36 1 10000 solver_scaling cg gamg -pc_gamg_square_graph 4 -pc_gamg_threshold 0.06
mpiexec -n 42 ./cuboid ../cuboid_settings.py 42 1 10000 solver_scaling cg gamg -pc_gamg_square_graph 4 -pc_gamg_threshold 0.06
mpiexec -n 48 ./cuboid ../cuboid_settings.py 48 1 10000 solver_scaling cg gamg -pc_gamg_square_graph 4 -pc_gamg_threshold 0.06
mpiexec -n 54 ./cuboid ../cuboid_settings.py 54 1 10000 solver_scaling cg gamg -pc_gamg_square_graph 4 -pc_gamg_threshold 0.06
mpiexec -n 60 ./cuboid ../cuboid_settings.py 60 1 10000 solver_scaling cg gamg -pc_gamg_square_graph 4 -pc_gamg_threshold 0.06
mpiexec -n 66 ./cuboid ../cuboid_settings.py 66 1 10000 solver_scaling cg gamg -pc_gamg_square_graph 4 -pc_gamg_threshold 0.06
mpiexec -n 72 ./cuboid ../cuboid_settings.py 72 1 10000 solver_scaling cg gamg -pc_gamg_square_graph 4 -pc_gamg_threshold 0.06

mpiexec -n 1 ./cuboid ../cuboid_settings.py 1 1 10000 solver_scaling preonly gamg -pc_gamg_square_graph 4 -pc_gamg_threshold 0.06
mpiexec -n 4 ./cuboid ../cuboid_settings.py 4 1 10000 solver_scaling preonly gamg -pc_gamg_square_graph 4 -pc_gamg_threshold 0.06
mpiexec -n 8 ./cuboid ../cuboid_settings.py 8 1 10000 solver_scaling preonly gamg -pc_gamg_square_graph 4 -pc_gamg_threshold 0.06
mpiexec -n 18 ./cuboid ../cuboid_settings.py 18 1 10000 solver_scaling preonly gamg -pc_gamg_square_graph 4 -pc_gamg_threshold 0.06
mpiexec -n 30 ./cuboid ../cuboid_settings.py 30 1 10000 solver_scaling preonly gamg -pc_gamg_square_graph 4 -pc_gamg_threshold 0.06
mpiexec -n 36 ./cuboid ../cuboid_settings.py 36 1 10000 solver_scaling preonly gamg -pc_gamg_square_graph 4 -pc_gamg_threshold 0.06
mpiexec -n 42 ./cuboid ../cuboid_settings.py 42 1 10000 solver_scaling preonly gamg -pc_gamg_square_graph 4 -pc_gamg_threshold 0.06
mpiexec -n 48 ./cuboid ../cuboid_settings.py 48 1 10000 solver_scaling preonly gamg -pc_gamg_square_graph 4 -pc_gamg_threshold 0.06
mpiexec -n 54 ./cuboid ../cuboid_settings.py 54 1 10000 solver_scaling preonly gamg -pc_gamg_square_graph 4 -pc_gamg_threshold 0.06
mpiexec -n 60 ./cuboid ../cuboid_settings.py 60 1 10000 solver_scaling preonly gamg -pc_gamg_square_graph 4 -pc_gamg_threshold 0.06
mpiexec -n 66 ./cuboid ../cuboid_settings.py 66 1 10000 solver_scaling preonly gamg -pc_gamg_square_graph 4 -pc_gamg_threshold 0.06  
mpiexec -n 72 ./cuboid ../cuboid_settings.py 72 1 10000 solver_scaling preonly gamg -pc_gamg_square_graph 4 -pc_gamg_threshold 0.06 

