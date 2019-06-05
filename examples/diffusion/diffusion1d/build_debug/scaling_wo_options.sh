# MG strong scaling
# args: <n_elements> <solver> <preconditioner> <name> <petsc options>

mpiexec -n 1 ./diffusion_1d_mg settings_scaling.py 10000 cg pchypre cg_boomeramg -ksp_type cg -pc_type hypre -pc_hypre_type boomeramg
mpiexec -n 2 ./diffusion_1d_mg settings_scaling.py 10000 cg pchypre cg_boomeramg -ksp_type cg -pc_type hypre -pc_hypre_type boomeramg
mpiexec -n 4 ./diffusion_1d_mg settings_scaling.py 10000 cg pchypre cg_boomeramg -ksp_type cg -pc_type hypre -pc_hypre_type boomeramg
mpiexec -n 8 ./diffusion_1d_mg settings_scaling.py 10000 cg pchypre cg_boomeramg -ksp_type cg -pc_type hypre -pc_hypre_type boomeramg
mpiexec -n 18 ./diffusion_1d_mg settings_scaling.py 10000 cg pchypre cg_boomeramg -ksp_type cg -pc_type hypre -pc_hypre_type boomeramg
mpiexec -n 24 ./diffusion_1d_mg settings_scaling.py 10000 cg pchypre cg_boomeramg -ksp_type cg -pc_type hypre -pc_hypre_type boomeramg
mpiexec -n 30 ./diffusion_1d_mg settings_scaling.py 10000 cg pchypre cg_boomeramg -ksp_type cg -pc_type hypre -pc_hypre_type boomeramg
mpiexec -n 36 ./diffusion_1d_mg settings_scaling.py 10000 cg pchypre cg_boomeramg -ksp_type cg -pc_type hypre -pc_hypre_type boomeramg
mpiexec -n 42 ./diffusion_1d_mg settings_scaling.py 10000 cg pchypre cg_boomeramg -ksp_type cg -pc_type hypre -pc_hypre_type boomeramg
mpiexec -n 48 ./diffusion_1d_mg settings_scaling.py 10000 cg pchypre cg_boomeramg -ksp_type cg -pc_type hypre -pc_hypre_type boomeramg
mpiexec -n 54 ./diffusion_1d_mg settings_scaling.py 10000 cg pchypre cg_boomeramg -ksp_type cg -pc_type hypre -pc_hypre_type boomeramg
mpiexec -n 60 ./diffusion_1d_mg settings_scaling.py 10000 cg pchypre cg_boomeramg -ksp_type cg -pc_type hypre -pc_hypre_type boomeramg
mpiexec -n 66 ./diffusion_1d_mg settings_scaling.py 10000 cg pchypre cg_boomeramg -ksp_type cg -pc_type hypre -pc_hypre_type boomeramg
mpiexec -n 72 ./diffusion_1d_mg settings_scaling.py 10000 cg pchypre cg_boomeramg -ksp_type cg -pc_type hypre -pc_hypre_type boomeramg

###########################################################

mpiexec -n 1 ./diffusion_1d_mg settings_scaling.py 10000 richardson pchypre boomeramg -ksp_type richardson -pc_type hypre -pc_hypre_type boomeramg
mpiexec -n 2 ./diffusion_1d_mg settings_scaling.py 10000 richardson pchypre boomeramg -ksp_type richardson -pc_type hypre -pc_hypre_type boomeramg
mpiexec -n 4 ./diffusion_1d_mg settings_scaling.py 10000 richardson pchypre boomeramg -ksp_type richardson -pc_type hypre -pc_hypre_type boomeramg
mpiexec -n 8 ./diffusion_1d_mg settings_scaling.py 10000 richardson pchypre boomeramg -ksp_type richardson -pc_type hypre -pc_hypre_type boomeramg
mpiexec -n 18 ./diffusion_1d_mg settings_scaling.py 10000 richardson pchypre boomeramg -ksp_type richardson -pc_type hypre -pc_hypre_type boomeramg
mpiexec -n 24 ./diffusion_1d_mg settings_scaling.py 10000 richardson pchypre boomeramg -ksp_type richardson -pc_type hypre -pc_hypre_type boomeramg
mpiexec -n 30 ./diffusion_1d_mg settings_scaling.py 10000 richardson pchypre boomeramg -ksp_type richardson -pc_type hypre -pc_hypre_type boomeramg
mpiexec -n 36 ./diffusion_1d_mg settings_scaling.py 10000 richardson pchypre boomeramg -ksp_type richardson -pc_type hypre -pc_hypre_type boomeramg
mpiexec -n 42 ./diffusion_1d_mg settings_scaling.py 10000 richardson pchypre boomeramg -ksp_type richardson -pc_type hypre -pc_hypre_type boomeramg
mpiexec -n 48 ./diffusion_1d_mg settings_scaling.py 10000 richardson pchypre boomeramg -ksp_type richardson -pc_type hypre -pc_hypre_type boomeramg
mpiexec -n 54 ./diffusion_1d_mg settings_scaling.py 10000 richardson pchypre boomeramg -ksp_type richardson -pc_type hypre -pc_hypre_type boomeramg
mpiexec -n 60 ./diffusion_1d_mg settings_scaling.py 10000 richardson pchypre boomeramg -ksp_type richardson -pc_type hypre -pc_hypre_type boomeramg
mpiexec -n 66 ./diffusion_1d_mg settings_scaling.py 10000 richardson pchypre boomeramg -ksp_type richardson -pc_type hypre -pc_hypre_type boomeramg
mpiexec -n 72 ./diffusion_1d_mg settings_scaling.py 10000 richardson pchypre boomeramg -ksp_type richardson -pc_type hypre -pc_hypre_type boomeramg

############################################################

mpiexec -n 1 ./diffusion_1d_mg settings_scaling.py 10000 cg gamg cg_gamg -ksp_type cg -pc_type gamg
mpiexec -n 2 ./diffusion_1d_mg settings_scaling.py 10000 cg gamg cg_gamg -ksp_type cg -pc_type gamg
mpiexec -n 4 ./diffusion_1d_mg settings_scaling.py 10000 cg gamg cg_gamg -ksp_type cg -pc_type gamg
mpiexec -n 8 ./diffusion_1d_mg settings_scaling.py 10000 cg gamg cg_gamg -ksp_type cg -pc_type gamg
mpiexec -n 18 ./diffusion_1d_mg settings_scaling.py 10000 cg gamg cg_gamg -ksp_type cg -pc_type gamg
mpiexec -n 24 ./diffusion_1d_mg settings_scaling.py 10000 cg gamg cg_gamg -ksp_type cg -pc_type gamg
mpiexec -n 30 ./diffusion_1d_mg settings_scaling.py 10000 cg gamg cg_gamg -ksp_type cg -pc_type gamg
mpiexec -n 36 ./diffusion_1d_mg settings_scaling.py 10000 cg gamg cg_gamg -ksp_type cg -pc_type gamg
mpiexec -n 42 ./diffusion_1d_mg settings_scaling.py 10000 cg gamg cg_gamg -ksp_type cg -pc_type gamg
mpiexec -n 48 ./diffusion_1d_mg settings_scaling.py 10000 cg gamg cg_gamg -ksp_type cg -pc_type gamg
mpiexec -n 54 ./diffusion_1d_mg settings_scaling.py 10000 cg gamg cg_gamg -ksp_type cg -pc_type gamg
mpiexec -n 60 ./diffusion_1d_mg settings_scaling.py 10000 cg gamg cg_gamg -ksp_type cg -pc_type gamg
mpiexec -n 66 ./diffusion_1d_mg settings_scaling.py 10000 cg gamg cg_gamg -ksp_type cg -pc_type gamg
mpiexec -n 72 ./diffusion_1d_mg settings_scaling.py 10000 cg gamg cg_gamg -ksp_type cg -pc_type gamg

#############################################################

mpiexec -n 1 ./diffusion_1d_mg settings_scaling.py 10000 preonly gamg gamg -ksp_type preonly -pc_type gamg
mpiexec -n 2 ./diffusion_1d_mg settings_scaling.py 10000 preonly gamg gamg -ksp_type preonly -pc_type gamg
mpiexec -n 4 ./diffusion_1d_mg settings_scaling.py 10000 preonly gamg gamg -ksp_type preonly -pc_type gamg
mpiexec -n 8 ./diffusion_1d_mg settings_scaling.py 10000 preonly gamg gamg -ksp_type preonly -pc_type gamg
mpiexec -n 18 ./diffusion_1d_mg settings_scaling.py 10000 preonly gamg gamg -ksp_type preonly -pc_type gamg
mpiexec -n 24 ./diffusion_1d_mg settings_scaling.py 10000 preonly gamg gamg -ksp_type preonly -pc_type gamg
mpiexec -n 30 ./diffusion_1d_mg settings_scaling.py 10000 preonly gamg gamg -ksp_type preonly -pc_type gamg
mpiexec -n 36 ./diffusion_1d_mg settings_scaling.py 10000 preonly gamg gamg -ksp_type preonly -pc_type gamg
mpiexec -n 42 ./diffusion_1d_mg settings_scaling.py 10000 preonly gamg gamg -ksp_type preonly -pc_type gamg
mpiexec -n 48 ./diffusion_1d_mg settings_scaling.py 10000 preonly gamg gamg -ksp_type preonly -pc_type gamg
mpiexec -n 54 ./diffusion_1d_mg settings_scaling.py 10000 preonly gamg gamg -ksp_type preonly -pc_type gamg
mpiexec -n 60 ./diffusion_1d_mg settings_scaling.py 10000 preonly gamg gamg -ksp_type preonly -pc_type gamg
mpiexec -n 66 ./diffusion_1d_mg settings_scaling.py 10000 preonly gamg gamg -ksp_type preonly -pc_type gamg
mpiexec -n 72 ./diffusion_1d_mg settings_scaling.py 10000 preonly gamg gamg -ksp_type preonly -pc_type gamg

#############################################################

mpiexec -n 1 ./diffusion_1d_mg settings_scaling.py 10000 gmres none gmres -ksp_type gmres -pc_type none
mpiexec -n 2 ./diffusion_1d_mg settings_scaling.py 10000 gmres none gmres -ksp_type gmres -pc_type none
mpiexec -n 4 ./diffusion_1d_mg settings_scaling.py 10000 gmres none gmres -ksp_type gmres -pc_type none
mpiexec -n 8 ./diffusion_1d_mg settings_scaling.py 10000 gmres none gmres -ksp_type gmres -pc_type none
mpiexec -n 18 ./diffusion_1d_mg settings_scaling.py 10000 gmres none gmres -ksp_type gmres -pc_type none
mpiexec -n 24 ./diffusion_1d_mg settings_scaling.py 10000 gmres none gmres -ksp_type gmres -pc_type none
mpiexec -n 30 ./diffusion_1d_mg settings_scaling.py 10000 gmres none gmres -ksp_type gmres -pc_type none
mpiexec -n 36 ./diffusion_1d_mg settings_scaling.py 10000 gmres none gmres -ksp_type gmres -pc_type none
mpiexec -n 42 ./diffusion_1d_mg settings_scaling.py 10000 gmres none gmres -ksp_type gmres -pc_type none
mpiexec -n 48 ./diffusion_1d_mg settings_scaling.py 10000 gmres none gmres -ksp_type gmres -pc_type none
mpiexec -n 54 ./diffusion_1d_mg settings_scaling.py 10000 gmres none gmres -ksp_type gmres -pc_type none
mpiexec -n 60 ./diffusion_1d_mg settings_scaling.py 10000 gmres none gmres -ksp_type gmres -pc_type none
mpiexec -n 66 ./diffusion_1d_mg settings_scaling.py 10000 gmres none gmres -ksp_type gmres -pc_type none
mpiexec -n 72 ./diffusion_1d_mg settings_scaling.py 10000 gmres none gmres -ksp_type gmres -pc_type none
