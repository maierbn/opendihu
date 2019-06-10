#!/bin/bash

#../run_tests.py  # this generates the following calls
#<n_processes_per_fiber> <n_fibers> <n_nodes_per_fiber> <scenario_name> <solver_type> <preconditioner_type>

#cg-hypre settings
#coarsentype
mpiexec -n 1 ./cuboid ../cuboid_settings.py 1 1 10000 cg_boomeramg_Falgout cg pchypre -pc_type hypre -pc_hypre_type boomeramg -pc_hypre_boomeramg_coarsen_type Falgout
mpiexec -n 1 ./cuboid ../cuboid_settings.py 1 1 10000 cg_boomeramg_HMIS cg pchypre -pc_type hypre -pc_hypre_type boomeramg -pc_hypre_boomeramg_coarsen_type HMIS
mpiexec -n 1 ./cuboid ../cuboid_settings.py 1 1 10000 cg_boomeramg_PMIS cg pchypre -pc_type hypre -pc_hypre_type boomeramg -pc_hypre_boomeramg_coarsen_type PMIS
#cycle_type
mpiexec -n 1 ./cuboid ../cuboid_settings.py 1 1 10000 cg_boomeramg_V cg pchypre -pc_type hypre -pc_hypre_type boomeramg -pc_hypre_boomeramg_cycle_type V
mpiexec -n 1 ./cuboid ../cuboid_settings.py 1 1 10000 cg_boomeramg_W cg pchypre -pc_type hypre -pc_hypre_type boomeramg -pc_hypre_boomeramg_cycle_type W
#agg coarsening
mpiexec -n 1 ./cuboid ../cuboid_settings.py 1 1 10000 cg_boomeramg_agg_path1 cg pchypre -pc_type hypre -pc_hypre_type boomeramg -pc_hypre_boomeramg_agg_num_paths 1
mpiexec -n 1 ./cuboid ../cuboid_settings.py 1 1 10000 cg_boomeramg_agg_path2 cg pchypre -pc_type hypre -pc_hypre_type boomeramg -pc_hypre_boomeramg_agg_num_paths 2
mpiexec -n 1 ./cuboid ../cuboid_settings.py 1 1 10000 cg_boomeramg_agg_path3 cg pchypre -pc_type hypre -pc_hypre_type boomeramg -pc_hypre_boomeramg_agg_num_paths 3
mpiexec -n 1 ./cuboid ../cuboid_settings.py 1 1 10000 cg_boomeramg_agg_path4 cg pchypre -pc_type hypre -pc_hypre_type boomeramg -pc_hypre_boomeramg_agg_num_paths 4
mpiexec -n 1 ./cuboid ../cuboid_settings.py 1 1 10000 cg_boomeramg_agg_path5 cg pchypre -pc_type hypre -pc_hypre_type boomeramg -pc_hypre_boomeramg_agg_num_paths 5

#boomeramg settings
#coarsentype
mpiexec -n 1 ./cuboid ../cuboid_settings.py 1 1 10000 boomeramg_classical richardson pchypre -pc_type hypre -pc_hypre_type boomeramg -pc_hypre_boomeramg_interp_type classical
mpiexec -n 1 ./cuboid ../cuboid_settings.py 1 1 10000 boomeramg_standard richardson pchypre -pc_type hypre -pc_hypre_type boomeramg -pc_hypre_boomeramg_interp_type standard
mpiexec -n 1 ./cuboid ../cuboid_settings.py 1 1 10000 boomeramg_exti richardson pchypre -pc_type hypre -pc_hypre_type boomeramg -pc_hypre_boomeramg_interp_type ext+i
#cycle_type
mpiexec -n 1 ./cuboid ../cuboid_settings.py 1 1 10000 boomeramg_V richardson pchypre -pc_type hypre -pc_hypre_type boomeramg -pc_hypre_boomeramg_cycle_type V
mpiexec -n 1 ./cuboid ../cuboid_settings.py 1 1 10000 boomeramg_W richardson pchypre -pc_type hypre -pc_hypre_type boomeramg -pc_hypre_boomeramg_cycle_type W
#agg coarsening
mpiexec -n 1 ./cuboid ../cuboid_settings.py 1 1 10000 boomeramg_agg1 richardson pchypre -pc_type hypre -pc_hypre_type boomeramg -pc_hypre_boomeramg_agg_num_paths 1
mpiexec -n 1 ./cuboid ../cuboid_settings.py 1 1 10000 boomeramg_agg2 richardson pchypre -pc_type hypre -pc_hypre_type boomeramg -pc_hypre_boomeramg_agg_num_paths 2
mpiexec -n 1 ./cuboid ../cuboid_settings.py 1 1 10000 boomeramg_agg3 richardson pchypre -pc_type hypre -pc_hypre_type boomeramg -pc_hypre_boomeramg_agg_num_paths 3
mpiexec -n 1 ./cuboid ../cuboid_settings.py 1 1 10000 boomeramg_agg4 richardson pchypre -pc_type hypre -pc_hypre_type boomeramg -pc_hypre_boomeramg_agg_num_paths 4
mpiexec -n 1 ./cuboid ../cuboid_settings.py 1 1 10000 boomeramg_agg5 richardson pchypre -pc_type hypre -pc_hypre_type boomeramg -pc_hypre_boomeramg_agg_num_paths 5

mpiexec -n 1 ./cuboid ../cuboid_settings.py 1 1 10000 boomeramg_2 richardson pchypre -pc_type hypre -pc_hypre_type boomeramg -pc_hypre_boomeramg_max_levels 2
mpiexec -n 1 ./cuboid ../cuboid_settings.py 1 1 10000 boomeramg_4 richardson pchypre -pc_type hypre -pc_hypre_type boomeramg -pc_hypre_boomeramg_max_levels 4
mpiexec -n 1 ./cuboid ../cuboid_settings.py 1 1 10000 boomeramg_6 richardson pchypre -pc_type hypre -pc_hypre_type boomeramg -pc_hypre_boomeramg_max_levels 6
mpiexec -n 1 ./cuboid ../cuboid_settings.py 1 1 10000 boomeramg_8 richardson pchypre -pc_type hypre -pc_hypre_type boomeramg -pc_hypre_boom eramg_max_levels 8
mpiexec -n 1 ./cuboid ../cuboid_settings.py 1 1 10000 boomeramg_10 richardson pchypre -pc_type hypre -pc_hypre_type boomeramg -pc_hypre_boomeramg_max_levels 10

mpiexec -n 1 ./cuboid ../cuboid_settings.py 1 1 10000 cg_boomeramg_2 cg pchypre -pc_type hypre -pc_hypre_type boomeramg -pc_hypre_boomeramg_max_levels 2
mpiexec -n 1 ./cuboid ../cuboid_settings.py 1 1 10000 cg_boomeramg_4 cg pchypre -pc_type hypre -pc_hypre_type boomeramg -pc_hypre_boomeramg_max_levels 4
mpiexec -n 1 ./cuboid ../cuboid_settings.py 1 1 10000 cg_boomeramg_6 cg pchypre -pc_type hypre -pc_hypre_type boomeramg -pc_hypre_boomeramg_max_levels 6
mpiexec -n 1 ./cuboid ../cuboid_settings.py 1 1 10000 cg_boomeramg_8 cg pchypre -pc_type hypre -pc_hypre_type boomeramg -pc_hypre_boomeramg_max_levels 8
mpiexec -n 1 ./cuboid ../cuboid_settings.py 1 1 10000 cg_boomeramg_10 cg pchypre -pc_type hypre -pc_hypre_type boomeramg -pc_hypre_boomeramg_max_levels 10