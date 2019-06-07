#!/bin/bash

#../run_tests.py  # this generates the following calls
#<n_processes_per_fiber> <n_fibers> <n_nodes_per_fiber> <scenario_name> <solver_type> <preconditioner_type>

#cg-hypre settings
#threshold
mpiexec -n 1 ./cuboid ../cuboid_settings.py 1 1 10000 cg_boomeramg_th0 cg pchypre -pc_type hypre -pc_hypre_type boomeramg -pc_hypre_boomeramg_strong_threshold 0.0
mpiexec -n 1 ./cuboid ../cuboid_settings.py 1 1 10000 cg_boomeramg_th1 cg pchypre -pc_type hypre -pc_hypre_type boomeramg -pc_hypre_boomeramg_strong_threshold 0.1
mpiexec -n 1 ./cuboid ../cuboid_settings.py 1 1 10000 cg_boomeramg_th2 cg pchypre -pc_type hypre -pc_hypre_type boomeramg -pc_hypre_boomeramg_strong_threshold 0.2
mpiexec -n 1 ./cuboid ../cuboid_settings.py 1 1 10000 cg_boomeramg_th3 cg pchypre -pc_type hypre -pc_hypre_type boomeramg -pc_hypre_boomeramg_strong_threshold 0.3
mpiexec -n 1 ./cuboid ../cuboid_settings.py 1 1 10000 cg_boomeramg_th4 cg pchypre -pc_type hypre -pc_hypre_type boomeramg -pc_hypre_boomeramg_strong_threshold 0.4
mpiexec -n 1 ./cuboid ../cuboid_settings.py 1 1 10000 cg_boomeramg_th5 cg pchypre -pc_type hypre -pc_hypre_type boomeramg -pc_hypre_boomeramg_strong_threshold 0.5
mpiexec -n 1 ./cuboid ../cuboid_settings.py 1 1 10000 cg_boomeramg_th6 cg pchypre -pc_type hypre -pc_hypre_type boomeramg -pc_hypre_boomeramg_strong_threshold 0.6
mpiexec -n 1 ./cuboid ../cuboid_settings.py 1 1 10000 cg_boomeramg_th7 cg pchypre -pc_type hypre -pc_hypre_type boomeramg -pc_hypre_boomeramg_strong_threshold 0.7
mpiexec -n 1 ./cuboid ../cuboid_settings.py 1 1 10000 cg_boomeramg_th8 cg pchypre -pc_type hypre -pc_hypre_type boomeramg -pc_hypre_boomeramg_strong_threshold 0.8
#coarsentype
mpiexec -n 1 ./cuboid ../cuboid_settings.py 1 1 10000 cg_boomeramg_Falgout cg pchypre -pc_type hypre -pc_hypre_type boomeramg -pc_hypre_boomeramg_coarsen_type Falgout
mpiexec -n 1 ./cuboid ../cuboid_settings.py 1 1 10000 cg_boomeramg_HMIS cg pchypre -pc_type hypre -pc_hypre_type boomeramg -pc_hypre_boomeramg_coarsen_type HMIS
mpiexec -n 1 ./cuboid ../cuboid_settings.py 1 1 10000 cg_boomeramg_PMIS cg pchypre -pc_type hypre -pc_hypre_type boomeramg -pc_hypre_boomeramg_coarsen_type PMIS
#cycle_type
mpiexec -n 1 ./cuboid ../cuboid_settings.py 1 1 10000 cg_boomeramg_V cg pchypre -pc_type hypre -pc_hypre_type boomeramg -pc_hypre_boomeramg_cycle_type V
mpiexec -n 1 ./cuboid ../cuboid_settings.py 1 1 10000 cg_boomeramg_W cg pchypre -pc_type hypre -pc_hypre_type boomeramg -pc_hypre_boomeramg_cycle_type W
#agg coarsening
mpiexec -n 1 ./cuboid ../cuboid_settings.py 1 1 10000 cg_boomeramg_agg0 cg pchypre -pc_type hypre -pc_hypre_type boomeramg -pc_hypre_boomeramg_agg_nl 0
mpiexec -n 1 ./cuboid ../cuboid_settings.py 1 1 10000 cg_boomeramg_agg1 cg pchypre -pc_type hypre -pc_hypre_type boomeramg -pc_hypre_boomeramg_agg_nl 1
mpiexec -n 1 ./cuboid ../cuboid_settings.py 1 1 10000 cg_boomeramg_agg2 cg pchypre -pc_type hypre -pc_hypre_type boomeramg -pc_hypre_boomeramg_agg_nl 2
mpiexec -n 1 ./cuboid ../cuboid_settings.py 1 1 10000 cg_boomeramg_agg3 cg pchypre -pc_type hypre -pc_hypre_type boomeramg -pc_hypre_boomeramg_agg_nl 3
mpiexec -n 1 ./cuboid ../cuboid_settings.py 1 1 10000 cg_boomeramg_agg4 cg pchypre -pc_type hypre -pc_hypre_type boomeramg -pc_hypre_boomeramg_agg_nl 4

#boomeramg settings
#threshold
mpiexec -n 1 ./cuboid ../cuboid_settings.py 1 1 10000 boomeramg_th0 richardson pchypre -pc_type hypre -pc_hypre_type boomeramg -pc_hypre_boomeramg_strong_threshold 0.0
mpiexec -n 1 ./cuboid ../cuboid_settings.py 1 1 10000 boomeramg_th1 richardson pchypre -pc_type hypre -pc_hypre_type boomeramg -pc_hypre_boomeramg_strong_threshold 0.1
mpiexec -n 1 ./cuboid ../cuboid_settings.py 1 1 10000 boomeramg_th2 richardson pchypre -pc_type hypre -pc_hypre_type boomeramg -pc_hypre_boomeramg_strong_threshold 0.2
mpiexec -n 1 ./cuboid ../cuboid_settings.py 1 1 10000 boomeramg_th3 richardson pchypre -pc_type hypre -pc_hypre_type boomeramg -pc_hypre_boomeramg_strong_threshold 0.3
mpiexec -n 1 ./cuboid ../cuboid_settings.py 1 1 10000 boomeramg_th4 richardson pchypre -pc_type hypre -pc_hypre_type boomeramg -pc_hypre_boomeramg_strong_threshold 0.4
mpiexec -n 1 ./cuboid ../cuboid_settings.py 1 1 10000 boomeramg_th5 richardson pchypre -pc_type hypre -pc_hypre_type boomeramg -pc_hypre_boomeramg_strong_threshold 0.5
mpiexec -n 1 ./cuboid ../cuboid_settings.py 1 1 10000 boomeramg_th6 richardson pchypre -pc_type hypre -pc_hypre_type boomeramg -pc_hypre_boomeramg_strong_threshold 0.6
mpiexec -n 1 ./cuboid ../cuboid_settings.py 1 1 10000 boomeramg_th7 richardson pchypre -pc_type hypre -pc_hypre_type boomeramg -pc_hypre_boomeramg_strong_threshold 0.7
mpiexec -n 1 ./cuboid ../cuboid_settings.py 1 1 10000 boomeramg_th8 richardson pchypre -pc_type hypre -pc_hypre_type boomeramg -pc_hypre_boomeramg_strong_threshold 0.8
#coarsentype
mpiexec -n 1 ./cuboid ../cuboid_settings.py 1 1 10000 boomeramg_Falgout richardson pchypre -pc_type hypre -pc_hypre_type boomeramg -pc_hypre_boomeramg_coarsen_type Falgout
mpiexec -n 1 ./cuboid ../cuboid_settings.py 1 1 10000 boomeramg_HMIS richardson pchypre -pc_type hypre -pc_hypre_type boomeramg -pc_hypre_boomeramg_coarsen_type HMIS
mpiexec -n 1 ./cuboid ../cuboid_settings.py 1 1 10000 boomeramg_PMIS richardson pchypre -pc_type hypre -pc_hypre_type boomeramg -pc_hypre_boomeramg_coarsen_type PMIS
#cycle_type
mpiexec -n 1 ./cuboid ../cuboid_settings.py 1 1 10000 boomeramg_V richardson pchypre -pc_type hypre -pc_hypre_type boomeramg -pc_hypre_boomeramg_cycle_type V
mpiexec -n 1 ./cuboid ../cuboid_settings.py 1 1 10000 boomeramg_W richardson pchypre -pc_type hypre -pc_hypre_type boomeramg -pc_hypre_boomeramg_cycle_type W
#agg coarsening
mpiexec -n 1 ./cuboid ../cuboid_settings.py 1 1 10000 boomeramg_agg0 richardson pchypre -pc_type hypre -pc_hypre_type boomeramg -pc_hypre_boomeramg_agg_nl 0
mpiexec -n 1 ./cuboid ../cuboid_settings.py 1 1 10000 boomeramg_agg1 richardson pchypre -pc_type hypre -pc_hypre_type boomeramg -pc_hypre_boomeramg_agg_nl 1
mpiexec -n 1 ./cuboid ../cuboid_settings.py 1 1 10000 boomeramg_agg2 richardson pchypre -pc_type hypre -pc_hypre_type boomeramg -pc_hypre_boomeramg_agg_nl 2
mpiexec -n 1 ./cuboid ../cuboid_settings.py 1 1 10000 boomeramg_agg3 richardson pchypre -pc_type hypre -pc_hypre_type boomeramg -pc_hypre_boomeramg_agg_nl 3
mpiexec -n 1 ./cuboid ../cuboid_settings.py 1 1 10000 boomeramg_agg4 richardson pchypre -pc_type hypre -pc_hypre_type boomeramg -pc_hypre_boomeramg_agg_nl 4
#cg-gamg settingsFileName
#threshold
mpiexec -n 1 ./cuboid ../cuboid_settings.py 1 1 10000 cg_gamg_th cg gamg -pc_gamg_threshold 0.0
mpiexec -n 1 ./cuboid ../cuboid_settings.py 1 1 10000 cg_gamg_th cg gamg -pc_gamg_threshold 0.01
mpiexec -n 1 ./cuboid ../cuboid_settings.py 1 1 10000 cg_gamg_th cg gamg -pc_gamg_threshold 0.02
mpiexec -n 1 ./cuboid ../cuboid_settings.py 1 1 10000 cg_gamg_th cg gamg -pc_gamg_threshold 0.03
mpiexec -n 1 ./cuboid ../cuboid_settings.py 1 1 10000 cg_gamg_th cg gamg -pc_gamg_threshold 0.04
mpiexec -n 1 ./cuboid ../cuboid_settings.py 1 1 10000 cg_gamg_th cg gamg -pc_gamg_threshold 0.05
mpiexec -n 1 ./cuboid ../cuboid_settings.py 1 1 10000 cg_gamg_th cg gamg -pc_gamg_threshold 0.06
mpiexec -n 1 ./cuboid ../cuboid_settings.py 1 1 10000 cg_gamg_th cg gamg -pc_gamg_threshold 0.07
mpiexec -n 1 ./cuboid ../cuboid_settings.py 1 1 10000 cg_gamg_th cg gamg -pc_gamg_threshold 0.08
#squaregraph
mpiexec -n 1 ./cuboid ../cuboid_settings.py 1 1 10000 cg_gamg_sq0 cg gamg -pc_gamg_square_graph 0
mpiexec -n 1 ./cuboid ../cuboid_settings.py 1 1 10000 cg_gamg_sq1 cg gamg -pc_gamg_square_graph 1
mpiexec -n 1 ./cuboid ../cuboid_settings.py 1 1 10000 cg_gamg_sq2 cg gamg -pc_gamg_square_graph 2
mpiexec -n 1 ./cuboid ../cuboid_settings.py 1 1 10000 cg_gamg_sq3 cg gamg -pc_gamg_square_graph 3
mpiexec -n 1 ./cuboid ../cuboid_settings.py 1 1 10000 cg_gamg_sq4 cg gamg -pc_gamg_square_graph 4
#Cycle type
mpiexec -n 1 ./cuboid ../cuboid_settings.py 1 1 10000 cg_gamg_V cg gamg -pc_mg_cycles v
mpiexec -n 1 ./cuboid ../cuboid_settings.py 1 1 10000 cg_gamg_W cg gamg -pc_mg_cycles w


#gamg settingsFileName
#threshold
mpiexec -n 1 ./cuboid ../cuboid_settings.py 1 1 10000 gamg_th0 preonly gamg -pc_gamg_threshold 0.0
mpiexec -n 1 ./cuboid ../cuboid_settings.py 1 1 10000 gamg_th1 preonly gamg -pc_gamg_threshold 0.01
mpiexec -n 1 ./cuboid ../cuboid_settings.py 1 1 10000 gamg_th2 preonly gamg -pc_gamg_threshold 0.02
mpiexec -n 1 ./cuboid ../cuboid_settings.py 1 1 10000 gamg_th3 preonly gamg -pc_gamg_threshold 0.03
mpiexec -n 1 ./cuboid ../cuboid_settings.py 1 1 10000 gamg_th4 preonly gamg -pc_gamg_threshold 0.04
mpiexec -n 1 ./cuboid ../cuboid_settings.py 1 1 10000 gamg_th5 preonly gamg -pc_gamg_threshold 0.05
mpiexec -n 1 ./cuboid ../cuboid_settings.py 1 1 10000 gamg_th6 preonly gamg -pc_gamg_threshold 0.06
mpiexec -n 1 ./cuboid ../cuboid_settings.py 1 1 10000 gamg_th7 preonly gamg -pc_gamg_threshold 0.07
mpiexec -n 1 ./cuboid ../cuboid_settings.py 1 1 10000 gamg_th8 preonly gamg -pc_gamg_threshold 0.08

#squaregraph
mpiexec -n 1 ./cuboid ../cuboid_settings.py 1 1 10000 gamg_sq0 preonly gamg -pc_gamg_square_graph 0
mpiexec -n 1 ./cuboid ../cuboid_settings.py 1 1 10000 gamg_sq1 preonly gamg -pc_gamg_square_graph 1
mpiexec -n 1 ./cuboid ../cuboid_settings.py 1 1 10000 gamg_sq2 preonly gamg -pc_gamg_square_graph 2
mpiexec -n 1 ./cuboid ../cuboid_settings.py 1 1 10000 gamg_sq3 preonly gamg -pc_gamg_square_graph 3
mpiexec -n 1 ./cuboid ../cuboid_settings.py 1 1 10000 gamg_sq4 preonly gamg -pc_gamg_square_graph 4
#cycle type
mpiexec -n 1 ./cuboid ../cuboid_settings.py 1 1 10000 gamg_V cg preonly -pc_mg_cycles v
mpiexec -n 1 ./cuboid ../cuboid_settings.py 1 1 10000 gamg_W cg preonly -pc_mg_cycles w