#!/bin/bash

#../run_tests.py  # this generates the following calls
#<n_processes_per_fiber> <n_fibers> <n_nodes_per_fiber> <scenario_name> <solver_type> <preconditioner_type>

#threshold
mpiexec -n 1 ./cuboid ../cuboid_settings.py 1 1 10000 cg_gamg_th5 cg gamg V 25 -pc_gamg_threshold 0.05
mpiexec -n 1 ./cuboid ../cuboid_settings.py 1 1 10000 cg_gamg_th52 cg gamg V 25 -pc_gamg_threshold 0.052
mpiexec -n 1 ./cuboid ../cuboid_settings.py 1 1 10000 cg_gamg_th54 cg gamg V 25 -pc_gamg_threshold 0.054
mpiexec -n 1 ./cuboid ../cuboid_settings.py 1 1 10000 cg_gamg_th56 cg gamg V 25 -pc_gamg_threshold 0.056
mpiexec -n 1 ./cuboid ../cuboid_settings.py 1 1 10000 cg_gamg_th58 cg gamg V 25 -pc_gamg_threshold 0.058
mpiexec -n 1 ./cuboid ../cuboid_settings.py 1 1 10000 cg_gamg_th6 cg gamg V 25 -pc_gamg_threshold 0.06
mpiexec -n 1 ./cuboid ../cuboid_settings.py 1 1 10000 cg_gamg_th62 cg gamg V 25 -pc_gamg_threshold 0.062
mpiexec -n 1 ./cuboid ../cuboid_settings.py 1 1 10000 cg_gamg_th64 cg gamg V 25 -pc_gamg_threshold 0.064
#Cycle type

mpiexec -n 1 ./cuboid ../cuboid_settings.py 1 1 10000 cg_gamg_V2 cg gamg V 2
mpiexec -n 1 ./cuboid ../cuboid_settings.py 1 1 10000 cg_gamg_V4 cg gamg V 4
mpiexec -n 1 ./cuboid ../cuboid_settings.py 1 1 10000 cg_gamg_V6 cg gamg V 6
mpiexec -n 1 ./cuboid ../cuboid_settings.py 1 1 10000 cg_gamg_V8 cg gamg V 8
mpiexec -n 1 ./cuboid ../cuboid_settings.py 1 1 10000 cg_gamg_V10 cg gamg V 10

mpiexec -n 1 ./cuboid ../cuboid_settings.py 1 1 10000 cg_gamg_W2 cg gamg W 2
mpiexec -n 1 ./cuboid ../cuboid_settings.py 1 1 10000 cg_gamg_W4 cg gamg W 4
mpiexec -n 1 ./cuboid ../cuboid_settings.py 1 1 10000 cg_gamg_W6 cg gamg W 6
mpiexec -n 1 ./cuboid ../cuboid_settings.py 1 1 10000 cg_gamg_W8 cg gamg W 8
mpiexec -n 1 ./cuboid ../cuboid_settings.py 1 1 10000 cg_gamg_W10 cg gamg W 10

#gamg settingsFileName
#threshold
mpiexec -n 1 ./cuboid ../cuboid_settings.py 1 1 10000 gamg_th5 preonly gamg V 25 -pc_gamg_threshold 0.05
mpiexec -n 1 ./cuboid ../cuboid_settings.py 1 1 10000 gamg_th52 preonly gamg V 25 -pc_gamg_threshold 0.052
mpiexec -n 1 ./cuboid ../cuboid_settings.py 1 1 10000 gamg_th54 preonly gamg V 25 -pc_gamg_threshold 0.054
mpiexec -n 1 ./cuboid ../cuboid_settings.py 1 1 10000 gamg_th56 preonly gamg V 25 -pc_gamg_threshold 0.056
mpiexec -n 1 ./cuboid ../cuboid_settings.py 1 1 10000 gamg_th58 preonly gamg V 25 -pc_gamg_threshold 0.058
mpiexec -n 1 ./cuboid ../cuboid_settings.py 1 1 10000 gamg_th6 preonly gamg V 25 -pc_gamg_threshold 0.06
mpiexec -n 1 ./cuboid ../cuboid_settings.py 1 1 10000 gamg_th62 preonly gamg V 25 -pc_gamg_threshold 0.062
mpiexec -n 1 ./cuboid ../cuboid_settings.py 1 1 10000 gamg_th64 preonly gamg V 25 -pc_gamg_threshold 0.064

#cycle type

mpiexec -n 1 ./cuboid ../cuboid_settings.py 1 1 10000 gamg_V2 richardson preonly V 2
mpiexec -n 1 ./cuboid ../cuboid_settings.py 1 1 10000 gamg_V4 richardson preonly V 4
mpiexec -n 1 ./cuboid ../cuboid_settings.py 1 1 10000 gamg_V6 richardson preonly V 6
mpiexec -n 1 ./cuboid ../cuboid_settings.py 1 1 10000 gamg_V8 richardson preonly V 8
mpiexec -n 1 ./cuboid ../cuboid_settings.py 1 1 10000 gamg_V10 richardson preonly V 10

mpiexec -n 1 ./cuboid ../cuboid_settings.py 1 1 10000 gamg_W2 richardson preonly W 2
mpiexec -n 1 ./cuboid ../cuboid_settings.py 1 1 10000 gamg_W4 richardson preonly W 4
mpiexec -n 1 ./cuboid ../cuboid_settings.py 1 1 10000 gamg_W6 richardson preonly W 6
mpiexec -n 1 ./cuboid ../cuboid_settings.py 1 1 10000 gamg_W8 richardson preonly W 8
mpiexec -n 1 ./cuboid ../cuboid_settings.py 1 1 10000 gamg_W10 richardson preonly W 10