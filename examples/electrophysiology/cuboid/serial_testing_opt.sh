#!/bin/bash

#../run_tests.py  # this generates the following calls
#<n_processes_per_fiber> <n_fibers> <n_nodes_per_fiber> <scenario_name> <solver_type> <preconditioner_type>

mpiexec -n 1 ./cuboid ../cuboid_settings.py 1 1 100 solver_scaling gmres none
mpiexec -n 1 ./cuboid ../cuboid_settings.py 1 1 1000 solver_scaling gmres none
mpiexec -n 1 ./cuboid ../cuboid_settings.py 1 1 10000 solver_scaling gmres none
mpiexec -n 1 ./cuboid ../cuboid_settings.py 1 1 100000 solver_scaling gmres none

mpiexec -n 1 ./cuboid ../cuboid_settings.py 1 1 100 solver_scaling cg pchypre -pc_type hypre -pc_hypre_type boomeramg -pc_hypre_boomeramg_coarsen_type HMIS -pc_hypre_boomeramg_agg_nl 2 -pc_hypre_boomeramg_strong_threshold 0.5 -pc_hypre_boomeramg_interp_type ext+i -pc_hypre_boomeramg_agg_num_paths 3 -pc_hypre_boom eramg_max_levels 8
mpiexec -n 1 ./cuboid ../cuboid_settings.py 1 1 1000 solver_scaling cg pchypre -pc_type hypre -pc_hypre_type boomeramg -pc_hypre_boomeramg_coarsen_type HMIS -pc_hypre_boomeramg_agg_nl 2 -pc_hypre_boomeramg_strong_threshold 0.5 -pc_hypre_boomeramg_interp_type ext+i -pc_hypre_boomeramg_agg_num_paths 3 -pc_hypre_boom eramg_max_levels 8
mpiexec -n 1 ./cuboid ../cuboid_settings.py 1 1 10000 solver_scaling cg pchypre -pc_type hypre -pc_hypre_type boomeramg -pc_hypre_boomeramg_coarsen_type HMIS -pc_hypre_boomeramg_agg_nl 2 -pc_hypre_boomeramg_strong_threshold 0.5 -pc_hypre_boomeramg_interp_type ext+i -pc_hypre_boomeramg_agg_num_paths 3 -pc_hypre_boom eramg_max_levels 8
mpiexec -n 1 ./cuboid ../cuboid_settings.py 1 1 100000 solver_scaling cg pchypre -pc_type hypre -pc_hypre_type boomeramg -pc_hypre_boomeramg_coarsen_type HMIS -pc_hypre_boomeramg_agg_nl 2 -pc_hypre_boomeramg_strong_threshold 0.5 -pc_hypre_boomeramg_interp_type ext+i -pc_hypre_boomeramg_agg_num_paths 3 -pc_hypre_boom eramg_max_levels 8

mpiexec -n 1 ./cuboid ../cuboid_settings.py 1 1 100 solver_scaling richardson pchypre -pc_type hypre -pc_hypre_type boomeramg -pc_hypre_boomeramg_coarsen_type Falgout -pc_hypre_boomeramg_agg_nl 0 -pc_hypre_boomeramg_strong_threshold 0.25 -pc_hypre_boomeramg_agg_num_paths 4
mpiexec -n 1 ./cuboid ../cuboid_settings.py 1 1 1000 solver_scaling richardson pchypre -pc_type hypre -pc_hypre_type boomeramg -pc_hypre_boomeramg_coarsen_type Falgout -pc_hypre_boomeramg_agg_nl 0 -pc_hypre_boomeramg_strong_threshold 0.25 -pc_hypre_boomeramg_agg_num_paths 4
mpiexec -n 1 ./cuboid ../cuboid_settings.py 1 1 10000 solver_scaling richardson pchypre -pc_type hypre -pc_hypre_type boomeramg -pc_hypre_boomeramg_coarsen_type Falgout -pc_hypre_boomeramg_agg_nl 0 -pc_hypre_boomeramg_strong_threshold 0.25 -pc_hypre_boomeramg_agg_num_paths 4
mpiexec -n 1 ./cuboid ../cuboid_settings.py 1 1 100000 solver_scaling richardson pchypre -pc_type hypre -pc_hypre_type boomeramg -pc_hypre_boomeramg_coarsen_type Falgout -pc_hypre_boomeramg_agg_nl 0 -pc_hypre_boomeramg_strong_threshold 0.25 -pc_hypre_boomeramg_agg_num_paths 4

mpiexec -n 1 ./cuboid ../cuboid_settings.py 1 1 100 solver_scaling cg gamg -pc_gamg_square_graph 4 -pc_gamg_threshold 0.06 -pc_mg_levels 8 
mpiexec -n 1 ./cuboid ../cuboid_settings.py 1 1 1000 solver_scaling cg gamg -pc_gamg_square_graph 4 -pc_gamg_threshold 0.06 -pc_mg_levels 8
mpiexec -n 1 ./cuboid ../cuboid_settings.py 1 1 10000 solver_scaling cg gamg -pc_gamg_square_graph 4 -pc_gamg_threshold 0.06 -pc_mg_levels 8
mpiexec -n 1 ./cuboid ../cuboid_settings.py 1 1 100000 solver_scaling cg gamg -pc_gamg_square_graph 4 -pc_gamg_threshold 0.06 -pc_mg_levels 8

mpiexec -n 1 ./cuboid ../cuboid_settings.py 1 1 100 solver_scaling preonly gamg -pc_gamg_square_graph 4 -pc_gamg_threshold 0.06 -pc_mg_levels 8
mpiexec -n 1 ./cuboid ../cuboid_settings.py 1 1 1000 solver_scaling preonly gamg -pc_gamg_square_graph 4 -pc_gamg_threshold 0.06 -pc_mg_levels 8
mpiexec -n 1 ./cuboid ../cuboid_settings.py 1 1 10000 solver_scaling preonly gamg -pc_gamg_square_graph 4 -pc_gamg_threshold 0.06 -pc_mg_levels 8
mpiexec -n 1 ./cuboid ../cuboid_settings.py 1 1 100000 solver_scaling preonly gamg -pc_gamg_square_graph 4 -pc_gamg_threshold 0.06 -pc_mg_levels 8
