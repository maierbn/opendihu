#!/bin/bash

#../run_tests.py  # this generates the following calls
#<n_processes_per_fiber> <n_fibers> <n_nodes_per_fiber> <scenario_name> <solver_type> <preconditioner_type>

#cg-hypre settings
#threshold
mpiexec -n 1 ./cuboid ../cuboid_settings.py 1 1 10000 cg_boomeramg cg pchypre -pc_type hypre -pc_hypre_type boomeramg

#coarsentype

#cycle_type

#

#boomeramg settings
#threshold
mpiexec -n 1 ./cuboid ../cuboid_settings.py 1 1 10000 boomeramg richardson pchypre -pc_type hypre -pc_hypre_type boomeramg

#coarsentype

#cycle_type

#

#cg-gamg settingsFileName
#threshold
mpiexec -n 1 ./cuboid ../cuboid_settings.py 1 1 10000 cg_gamg cg gamg

#squaregraph

#cycle_type

#

#gamg settingsFileName
#threshold
mpiexec -n 1 ./cuboid ../cuboid_settings.py 1 1 10000 boomeramg preonly gamg

#squaregraph

#cycle_type

#