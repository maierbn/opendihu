#!/bin/bash

#../run_tests.py  # this generates the following calls

# hardcoded calls for hazel hen
# arguments: <n_processes_per_fiber> <n_fibers> <n_nodes_per_fiber> <scenario_name> <solver_type>

while true; do

# 1 node
if [ $nnodes = "1" ]; then
aprun -N 1 -n 1 ./cuboid ../cuboid_settings.py 1 1 100 solver_scaling lu      | tee -a out_solver_scaling_${nnodes}.txt 2>&1 # 1 nodes
aprun -N 1 -n 1 ./cuboid ../cuboid_settings.py 1 1 100 solver_scaling gmres   | tee -a out_solver_scaling_${nnodes}.txt 2>&1 # 1 nodes
aprun -N 1 -n 1 ./cuboid ../cuboid_settings.py 1 1 100 solver_scaling cg      | tee -a out_solver_scaling_${nnodes}.txt 2>&1 # 1 nodes

aprun -N 4 -n 4 ./cuboid ../cuboid_settings.py 4 1 400 solver_scaling lu      | tee -a out_solver_scaling_${nnodes}.txt 2>&1 # 1 nodes
aprun -N 4 -n 4 ./cuboid ../cuboid_settings.py 4 1 400 solver_scaling gmres   | tee -a out_solver_scaling_${nnodes}.txt 2>&1 # 1 nodes
aprun -N 4 -n 4 ./cuboid ../cuboid_settings.py 4 1 400 solver_scaling cg      | tee -a out_solver_scaling_${nnodes}.txt 2>&1 # 1 nodes

aprun -N 10 -n 10 ./cuboid ../cuboid_settings.py 10 1 1000 solver_scaling lu      | tee -a out_solver_scaling_${nnodes}.txt 2>&1 # 1 nodes
aprun -N 10 -n 10 ./cuboid ../cuboid_settings.py 10 1 1000 solver_scaling gmres   | tee -a out_solver_scaling_${nnodes}.txt 2>&1 # 1 nodes
aprun -N 10 -n 10 ./cuboid ../cuboid_settings.py 10 1 1000 solver_scaling cg      | tee -a out_solver_scaling_${nnodes}.txt 2>&1 # 1 nodes
fi

# 2 nodes
if [ $nnodes = "2" ]; then
aprun -N 24 -n 32 ./cuboid ../cuboid_settings.py 32 1 3200 solver_scaling lu      | tee -a out_solver_scaling_${nnodes}.txt 2>&1 # 2 nodes
aprun -N 24 -n 32 ./cuboid ../cuboid_settings.py 32 1 3200 solver_scaling gmres   | tee -a out_solver_scaling_${nnodes}.txt 2>&1 # 2 nodes
aprun -N 24 -n 32 ./cuboid ../cuboid_settings.py 32 1 3200 solver_scaling cg      | tee -a out_solver_scaling_${nnodes}.txt 2>&1 # 2 nodes
fi

# 14 nodes
if [ $nnodes = "14" ]; then
aprun -N 24 -n 100 ./cuboid ../cuboid_settings.py 100 1 10000 solver_scaling lu      | tee -a out_solver_scaling_${nnodes}.txt 2>&1 # 5 nodes
aprun -N 24 -n 100 ./cuboid ../cuboid_settings.py 100 1 10000 solver_scaling gmres   | tee -a out_solver_scaling_${nnodes}.txt 2>&1 # 5 nodes
aprun -N 24 -n 100 ./cuboid ../cuboid_settings.py 100 1 10000 solver_scaling cg      | tee -a out_solver_scaling_${nnodes}.txt 2>&1 # 5 nodes

aprun -N 24 -n 317 ./cuboid ../cuboid_settings.py 317 1 31700 solver_scaling lu      | tee -a out_solver_scaling_${nnodes}.txt 2>&1 # 14 nodes
aprun -N 24 -n 317 ./cuboid ../cuboid_settings.py 317 1 31700 solver_scaling gmres   | tee -a out_solver_scaling_${nnodes}.txt 2>&1 # 14 nodes
aprun -N 24 -n 317 ./cuboid ../cuboid_settings.py 317 1 31700 solver_scaling cg      | tee -a out_solver_scaling_${nnodes}.txt 2>&1 # 14 nodes
fi

# 42 nodes
if [ $nnodes = "42" ]; then
aprun -N 24 -n 1000 ./cuboid ../cuboid_settings.py 1000 1 100000 solver_scaling lu      | tee -a out_solver_scaling_${nnodes}.txt 2>&1 # 42 nodes
aprun -N 24 -n 1000 ./cuboid ../cuboid_settings.py 1000 1 100000 solver_scaling gmres   | tee -a out_solver_scaling_${nnodes}.txt 2>&1 # 42 nodes
aprun -N 24 -n 1000 ./cuboid ../cuboid_settings.py 1000 1 100000 solver_scaling cg      | tee -a out_solver_scaling_${nnodes}.txt 2>&1 # 42 nodes
fi

# 132 nodes
if [ $nnodes = "132" ]; then
aprun -N 24 -n 3163 ./cuboid ../cuboid_settings.py 3163 1 316300 solver_scaling lu      | tee -a out_solver_scaling_${nnodes}.txt 2>&1 # 132 nodes
aprun -N 24 -n 3163 ./cuboid ../cuboid_settings.py 3163 1 316300 solver_scaling gmres   | tee -a out_solver_scaling_${nnodes}.txt 2>&1 # 132 nodes
aprun -N 24 -n 3163 ./cuboid ../cuboid_settings.py 3163 1 316300 solver_scaling cg      | tee -a out_solver_scaling_${nnodes}.txt 2>&1 # 132 nodes
fi

# 417 nodes
if [ $nnodes = "417" ]; then
aprun -N 24 -n 10000 ./cuboid ../cuboid_settings.py 10000 1 1000000 solver_scaling lu      | tee -a out_solver_scaling_${nnodes}.txt 2>&1 # 417 nodes
aprun -N 24 -n 10000 ./cuboid ../cuboid_settings.py 10000 1 1000000 solver_scaling gmres   | tee -a out_solver_scaling_${nnodes}.txt 2>&1 # 417 nodes
aprun -N 24 -n 10000 ./cuboid ../cuboid_settings.py 10000 1 1000000 solver_scaling cg      | tee -a out_solver_scaling_${nnodes}.txt 2>&1 # 417 nodes
fi

done

cd $workdir
