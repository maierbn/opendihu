# MG weak scaling
# args: <n_elements> <solver> <preconditioner> <name> <petsc options>


#CG-Boomeramg settigns testing

# threshold
mpiexec -n 1 ./diffusion_1d_mg settings_scaling.py 10000 cg pchypre cg_boomearmg_th0 -pc_type hypre -pc_hypre_type boomeramg -pc_hypre_boomeramg_strong_threshold 0.
mpiexec -n 1 ./diffusion_1d_mg settings_scaling.py 10000 cg pchypre cg_boomearmg_th1 -pc_type hypre -pc_hypre_type boomeramg -pc_hypre_boomeramg_strong_threshold 0.1
mpiexec -n 1 ./diffusion_1d_mg settings_scaling.py 10000 cg pchypre cg_boomearmg_th2 -pc_type hypre -pc_hypre_type boomeramg  -pc_hypre_boomeramg_strong_threshold 0.2
mpiexec -n 1 ./diffusion_1d_mg settings_scaling.py 10000 cg pchypre cg_boomearmg_th3 -pc_type hypre -pc_hypre_type boomeramg  -pc_hypre_boomeramg_strong_threshold 0.3
mpiexec -n 1 ./diffusion_1d_mg settings_scaling.py 10000 cg pchypre cg_boomearmg_th4 -pc_type hypre -pc_hypre_type boomeramg  -pc_hypre_boomeramg_strong_threshold 0.4
mpiexec -n 1 ./diffusion_1d_mg settings_scaling.py 10000 cg pchypre cg_boomearmg_th5 -pc_type hypre -pc_hypre_type boomeramg  -pc_hypre_boomeramg_strong_threshold 0.5
mpiexec -n 1 ./diffusion_1d_mg settings_scaling.py 10000 cg pchypre cg_boomearmg_th6 -pc_type hypre -pc_hypre_type boomeramg  -pc_hypre_boomeramg_strong_threshold 0.6
mpiexec -n 1 ./diffusion_1d_mg settings_scaling.py 10000 cg pchypre cg_boomearmg_th7 -pc_type hypre -pc_hypre_type boomeramg  -pc_hypre_boomeramg_strong_threshold 0.7
mpiexec -n 1 ./diffusion_1d_mg settings_scaling.py 10000 cg pchypre cg_boomearmg_th8 -pc_type hypre -pc_hypre_type boomeramg  -pc_hypre_boomeramg_strong_threshold 0.8



#coarsentype
mpiexec -n 1 ./diffusion_1d_mg settings_scaling.py 10000 cg pchypre cg_boomearmg_PMIS -pc_type hypre -pc_hypre_type boomeramg -pc_hypre_boomeramg_coarsen_type PMIS
mpiexec -n 1 ./diffusion_1d_mg settings_scaling.py 10000 cg pchypre cg_boomearmg_HMIS -pc_type hypre -pc_hypre_type boomeramg -pc_hypre_boomeramg_coarsen_type HMIS
mpiexec -n 1 ./diffusion_1d_mg settings_scaling.py 10000 cg pchypre cg_boomearmg_Falgout -pc_type hypre -pc_hypre_type boomeramg -pc_hypre_boomeramg_coarsen_type Falgout

#cycle_type
mpiexec -n 1 ./diffusion_1d_mg settings_scaling.py 10000 cg pchypre cg_boomeramg_V -pc_type hypre -pc_hypre_type boomeramg  -pc_hypre_boomeramg_cycle_type V
mpiexec -n 1 ./diffusion_1d_mg settings_scaling.py 10000 cg pchypre cg_boomeramg
_W -pc_type hypre -pc_hypre_type boomeramg -pc_hypre_boomeramg_cycle_type W

#aggressive Coarsening
mpiexec -n 1 ./diffusion_1d_mg settings_scaling.py 10000 cg pchypre cg_boomeramg_agglvl0  -pc_type hypre -pc_hypre_type boomeramg -pc_hypre_boomeramg_agg_nl 0
mpiexec -n 1 ./diffusion_1d_mg settings_scaling.py 10000 cg pchypre cg_boomeramg_agglvl1  -pc_type hypre -pc_hypre_type boomeramg -pc_hypre_boomeramg_agg_nl 1
mpiexec -n 1 ./diffusion_1d_mg settings_scaling.py 10000 cg pchypre cg_boomeramg_agglvl2  -pc_type hypre -pc_hypre_type boomeramg -pc_hypre_boomeramg_agg_nl 2
mpiexec -n 1 ./diffusion_1d_mg settings_scaling.py 10000 cg pchypre cg_boomeramg_agglvl3  -pc_type hypre -pc_hypre_type boomeramg -pc_hypre_boomeramg_agg_nl 3
mpiexec -n 1 ./diffusion_1d_mg settings_scaling.py 10000 cg pchypre cg_boomearmg_agglvl4  -pc_type hypre -pc_hypre_type boomeramg -pc_hypre_boomeramg_agg_nl 4

#Boomeramg settings testing

#threshold
mpiexec -n 1 ./diffusion_1d_mg settings_scaling.py 10000 richardson pchypre boomeramg_th0 -pc_type hypre -pc_hypre_type boomeramg  -pc_hypre_boomeramg_strong_threshold 0.
mpiexec -n 1 ./diffusion_1d_mg settings_scaling.py 10000 richardson pchypre boomeramg_th1 -pc_type hypre -pc_hypre_type boomeramg  -pc_hypre_boomeramg_strong_threshold 0.1
mpiexec -n 1 ./diffusion_1d_mg settings_scaling.py 10000 richardson pchypre boomeramg_th2 -pc_type hypre -pc_hypre_type boomeramg  -pc_hypre_boomeramg_strong_threshold 0.2
mpiexec -n 1 ./diffusion_1d_mg settings_scaling.py 10000 richardson pchypre boomeramg_th3 -pc_type hypre -pc_hypre_type boomeramg  -pc_hypre_boomeramg_strong_threshold 0.3
mpiexec -n 1 ./diffusion_1d_mg settings_scaling.py 10000 richardson pchypre boomeramg_th4 -pc_type hypre -pc_hypre_type boomeramg  -pc_hypre_boomeramg_strong_threshold 0.4
mpiexec -n 1 ./diffusion_1d_mg settings_scaling.py 10000 richardson pchypre boomeramg_th5 -pc_type hypre -pc_hypre_type boomeramg  -pc_hypre_boomeramg_strong_threshold 0.5
mpiexec -n 1 ./diffusion_1d_mg settings_scaling.py 10000 richardson pchypre boomeramg_th6 -pc_type hypre -pc_hypre_type boomeramg  -pc_hypre_boomeramg_strong_threshold 0.6
mpiexec -n 1 ./diffusion_1d_mg settings_scaling.py 10000 richardson pchypre boomeramg_th7 -pc_type hypre -pc_hypre_type boomeramg  -pc_hypre_boomeramg_strong_threshold 0.7
mpiexec -n 1 ./diffusion_1d_mg settings_scaling.py 10000 richardson pchypre boomeramg_th8 -pc_type hypre -pc_hypre_type boomeramg  -pc_hypre_boomeramg_strong_threshold 0.8

#coarsentype
mpiexec -n 1 ./diffusion_1d_mg settings_scaling.py 10000 richardson pchypre boomeramg_PMIS -pc_type hypre -pc_hypre_type boomeramg -pc_hypre_boomeramg_coarsen_type PMIS
mpiexec -n 1 ./diffusion_1d_mg settings_scaling.py 10000 richardson pchypre boomeramg_HMIS -pc_type hypre -pc_hypre_type boomeramg -pc_hypre_boomeramg_coarsen_type HMIS
mpiexec -n 1 ./diffusion_1d_mg settings_scaling.py 10000 richardson pchypre boomeramg_Falgout -pc_type hypre -pc_hypre_type boomeramg -pc_hypre_boomeramg_coarsen_type Falgout

#cycle_type
mpiexec -n 1 ./diffusion_1d_mg settings_scaling.py 10000 richardson pchypre boomeramg_V -pc_type hypre -pc_hypre_type boomeramg  -pc_hypre_boomeramg_cycle_type V
mpiexec -n 1 ./diffusion_1d_mg settings_scaling.py 10000 richardson pchypre boomeramg_W -pc_type hypre -pc_hypre_type boomeramg -pc_hypre_boomeramg_cycle_type W

#aggressive coarsening levels
mpiexec -n 1 ./diffusion_1d_mg settings_scaling.py 10000 richardson pchypre boomeramg_aggLvl0 -pc_type hypre -pc_hypre_type boomeramg -pc_hypre_boomeramg_agg_nl 0
mpiexec -n 1 ./diffusion_1d_mg settings_scaling.py 10000 richardson pchypre boomeramg_aggLvl1 -pc_type hypre -pc_hypre_type boomeramg -pc_hypre_boomeramg_agg_nl 1
mpiexec -n 1 ./diffusion_1d_mg settings_scaling.py 10000 richardson pchypre boomeramg_aggLvl2 -pc_type hypre -pc_hypre_type boomeramg -pc_hypre_boomeramg_agg_nl 2
mpiexec -n 1 ./diffusion_1d_mg settings_scaling.py 10000 richardson pchypre boomeramg_aggLvl3 -pc_type hypre -pc_hypre_type boomeramg -pc_hypre_boomeramg_agg_nl 3
mpiexec -n 1 ./diffusion_1d_mg settings_scaling.py 10000 richardson pchypre boomeramg_aggLvl4 -pc_type hypre -pc_hypre_type boomeramg -pc_hypre_boomeramg_agg_nl 4
#CG-gamg settings testing

#threshold
#mpiexec -n 1 ./diffusion_1d_mg settings_scaling.py 10000 cg gamg cg_gamg

#squaregraph

#symmetry

#cycle_type

#



#gamg settings testing

#threshold
#mpiexec -n 1 ./diffusion_1d_mg settings_scaling.py 10000 preonly gamg gamg

#squaregraph

#symmetry

#cycle_type

#

