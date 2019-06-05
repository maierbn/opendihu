# MG strong scaling
# args: <n_elements> <solver> <preconditioner> <name> <petsc options>

mpiexec -n 1 ./diffusion_1d_mg settings_scaling.py 100 cg pchypre cg_boomeramg -ksp_type cg -pc_type hypre -pc_hypre_type boomeramg
mpiexec -n 2 ./diffusion_1d_mg settings_scaling.py 200 cg pchypre cg_boomeramg -ksp_type cg -pc_type hypre -pc_hypre_type boomeramg
mpiexec -n 4 ./diffusion_1d_mg settings_scaling.py 400 cg pchypre cg_boomeramg -ksp_type cg -pc_type hypre -pc_hypre_type boomeramg
mpiexec -n 8 ./diffusion_1d_mg settings_scaling.py 800 cg pchypre cg_boomeramg -ksp_type cg -pc_type hypre -pc_hypre_type boomeramg
mpiexec -n 18 ./diffusion_1d_mg settings_scaling.py 1800 cg pchypre cg_boomeramg -ksp_type cg -pc_type hypre -pc_hypre_type boomeramg
mpiexec -n 24 ./diffusion_1d_mg settings_scaling.py 2400 cg pchypre cg_boomeramg -ksp_type cg -pc_type hypre -pc_hypre_type boomeramg
mpiexec -n 30 ./diffusion_1d_mg settings_scaling.py 3000 cg pchypre cg_boomeramg -ksp_type cg -pc_type hypre -pc_hypre_type boomeramg
mpiexec -n 36 ./diffusion_1d_mg settings_scaling.py 3600 cg pchypre cg_boomeramg -ksp_type cg -pc_type hypre -pc_hypre_type boomeramg
mpiexec -n 42 ./diffusion_1d_mg settings_scaling.py 4200 cg pchypre cg_boomeramg -ksp_type cg -pc_type hypre -pc_hypre_type boomeramg
mpiexec -n 48 ./diffusion_1d_mg settings_scaling.py 4800 cg pchypre cg_boomeramg -ksp_type cg -pc_type hypre -pc_hypre_type boomeramg
mpiexec -n 54 ./diffusion_1d_mg settings_scaling.py 5400 cg pchypre cg_boomeramg -ksp_type cg -pc_type hypre -pc_hypre_type boomeramg
mpiexec -n 60 ./diffusion_1d_mg settings_scaling.py 6000 cg pchypre cg_boomeramg -ksp_type cg -pc_type hypre -pc_hypre_type boomeramg
mpiexec -n 66 ./diffusion_1d_mg settings_scaling.py 6600 cg pchypre cg_boomeramg -ksp_type cg -pc_type hypre -pc_hypre_type boomeramg
mpiexec -n 72 ./diffusion_1d_mg settings_scaling.py 7200 cg pchypre cg_boomeramg -ksp_type cg -pc_type hypre -pc_hypre_type boomeramg

###########################################################

mpiexec -n 1 ./diffusion_1d_mg settings_scaling.py 100 richardson pchypre boomeramg -ksp_type richardson -pc_type hypre -pc_hypre_type boomeramg
mpiexec -n 2 ./diffusion_1d_mg settings_scaling.py 200 richardson pchypre boomeramg -ksp_type richardson -pc_type hypre -pc_hypre_type boomeramg
mpiexec -n 4 ./diffusion_1d_mg settings_scaling.py 400 richardson pchypre boomeramg -ksp_type richardson -pc_type hypre -pc_hypre_type boomeramg
mpiexec -n 8 ./diffusion_1d_mg settings_scaling.py 800 richardson pchypre boomeramg -ksp_type richardson -pc_type hypre -pc_hypre_type boomeramg
mpiexec -n 18 ./diffusion_1d_mg settings_scaling.py 1800 richardson pchypre boomeramg -ksp_type richardson -pc_type hypre -pc_hypre_type boomeramg
mpiexec -n 24 ./diffusion_1d_mg settings_scaling.py 2400 richardson pchypre boomeramg -ksp_type richardson -pc_type hypre -pc_hypre_type boomeramg
mpiexec -n 30 ./diffusion_1d_mg settings_scaling.py 3000 richardson pchypre boomeramg -ksp_type richardson -pc_type hypre -pc_hypre_type boomeramg
mpiexec -n 36 ./diffusion_1d_mg settings_scaling.py 3600 richardson pchypre boomeramg -ksp_type richardson -pc_type hypre -pc_hypre_type boomeramg
mpiexec -n 42 ./diffusion_1d_mg settings_scaling.py 4200 richardson pchypre boomeramg -ksp_type richardson -pc_type hypre -pc_hypre_type boomeramg
mpiexec -n 48 ./diffusion_1d_mg settings_scaling.py 4800 richardson pchypre boomeramg -ksp_type richardson -pc_type hypre -pc_hypre_type boomeramg
mpiexec -n 54 ./diffusion_1d_mg settings_scaling.py 5400 richardson pchypre boomeramg -ksp_type richardson -pc_type hypre -pc_hypre_type boomeramg
mpiexec -n 60 ./diffusion_1d_mg settings_scaling.py 6000 richardson pchypre boomeramg -ksp_type richardson -pc_type hypre -pc_hypre_type boomeramg
mpiexec -n 66 ./diffusion_1d_mg settings_scaling.py 6600 richardson pchypre boomeramg -ksp_type richardson -pc_type hypre -pc_hypre_type boomeramg
mpiexec -n 72 ./diffusion_1d_mg settings_scaling.py 7200 richardson pchypre boomeramg -ksp_type richardson -pc_type hypre -pc_hypre_type boomeramg

############################################################

mpiexec -n 1 ./diffusion_1d_mg settings_scaling.py 100 cg gamg cg_gamg -ksp_type cg -pc_type gamg
mpiexec -n 2 ./diffusion_1d_mg settings_scaling.py 200 cg gamg cg_gamg -ksp_type cg -pc_type gamg
mpiexec -n 4 ./diffusion_1d_mg settings_scaling.py 400 cg gamg cg_gamg -ksp_type cg -pc_type gamg
mpiexec -n 8 ./diffusion_1d_mg settings_scaling.py 800 cg gamg cg_gamg -ksp_type cg -pc_type gamg
mpiexec -n 18 ./diffusion_1d_mg settings_scaling.py 1800 cg gamg cg_gamg -ksp_type cg -pc_type gamg
mpiexec -n 24 ./diffusion_1d_mg settings_scaling.py 2400 cg gamg cg_gamg -ksp_type cg -pc_type gamg
mpiexec -n 30 ./diffusion_1d_mg settings_scaling.py 3000 cg gamg cg_gamg -ksp_type cg -pc_type gamg
mpiexec -n 36 ./diffusion_1d_mg settings_scaling.py 3600 cg gamg cg_gamg -ksp_type cg -pc_type gamg
mpiexec -n 42 ./diffusion_1d_mg settings_scaling.py 4200 cg gamg cg_gamg -ksp_type cg -pc_type gamg
mpiexec -n 48 ./diffusion_1d_mg settings_scaling.py 4800 cg gamg cg_gamg -ksp_type cg -pc_type gamg
mpiexec -n 54 ./diffusion_1d_mg settings_scaling.py 5400 cg gamg cg_gamg -ksp_type cg -pc_type gamg
mpiexec -n 60 ./diffusion_1d_mg settings_scaling.py 6000 cg gamg cg_gamg -ksp_type cg -pc_type gamg
mpiexec -n 66 ./diffusion_1d_mg settings_scaling.py 6600 cg gamg cg_gamg -ksp_type cg -pc_type gamg
mpiexec -n 72 ./diffusion_1d_mg settings_scaling.py 7200 cg gamg cg_gamg -ksp_type cg -pc_type gamg

#############################################################

mpiexec -n 1 ./diffusion_1d_mg settings_scaling.py 100 preonly gamg gamg -ksp_type preonly -pc_type gamg
mpiexec -n 2 ./diffusion_1d_mg settings_scaling.py 200 preonly gamg gamg -ksp_type preonly -pc_type gamg
mpiexec -n 4 ./diffusion_1d_mg settings_scaling.py 400 preonly gamg gamg -ksp_type preonly -pc_type gamg
mpiexec -n 8 ./diffusion_1d_mg settings_scaling.py 800 preonly gamg gamg -ksp_type preonly -pc_type gamg
mpiexec -n 18 ./diffusion_1d_mg settings_scaling.py 1800 preonly gamg gamg -ksp_type preonly -pc_type gamg
mpiexec -n 24 ./diffusion_1d_mg settings_scaling.py 2400 preonly gamg gamg -ksp_type preonly -pc_type gamg
mpiexec -n 30 ./diffusion_1d_mg settings_scaling.py 3000 preonly gamg gamg -ksp_type preonly -pc_type gamg
mpiexec -n 36 ./diffusion_1d_mg settings_scaling.py 3600 preonly gamg gamg -ksp_type preonly -pc_type gamg
mpiexec -n 42 ./diffusion_1d_mg settings_scaling.py 4200 preonly gamg gamg -ksp_type preonly -pc_type gamg
mpiexec -n 48 ./diffusion_1d_mg settings_scaling.py 4800 preonly gamg gamg -ksp_type preonly -pc_type gamg
mpiexec -n 54 ./diffusion_1d_mg settings_scaling.py 5400 preonly gamg gamg -ksp_type preonly -pc_type gamg
mpiexec -n 60 ./diffusion_1d_mg settings_scaling.py 6000 preonly gamg gamg -ksp_type preonly -pc_type gamg
mpiexec -n 66 ./diffusion_1d_mg settings_scaling.py 6600 preonly gamg gamg -ksp_type preonly -pc_type gamg
mpiexec -n 72 ./diffusion_1d_mg settings_scaling.py 7200 preonly gamg gamg -ksp_type preonly -pc_type gamg

#############################################################

mpiexec -n 1 ./diffusion_1d_mg settings_scaling.py 100 gmres none gmres -ksp_type gmres -pc_type none
mpiexec -n 2 ./diffusion_1d_mg settings_scaling.py 200 gmres none gmres -ksp_type gmres -pc_type none
mpiexec -n 4 ./diffusion_1d_mg settings_scaling.py 400 gmres none gmres -ksp_type gmres -pc_type none
mpiexec -n 8 ./diffusion_1d_mg settings_scaling.py 800 gmres none gmres -ksp_type gmres -pc_type none
mpiexec -n 18 ./diffusion_1d_mg settings_scaling.py 1800 gmres none gmres -ksp_type gmres -pc_type none
mpiexec -n 24 ./diffusion_1d_mg settings_scaling.py 2400 gmres none gmres -ksp_type gmres -pc_type none
mpiexec -n 30 ./diffusion_1d_mg settings_scaling.py 3000 gmres none gmres -ksp_type gmres -pc_type none
mpiexec -n 36 ./diffusion_1d_mg settings_scaling.py 3600 gmres none gmres -ksp_type gmres -pc_type none
mpiexec -n 42 ./diffusion_1d_mg settings_scaling.py 4200 gmres none gmres -ksp_type gmres -pc_type none
mpiexec -n 48 ./diffusion_1d_mg settings_scaling.py 4800 gmres none gmres -ksp_type gmres -pc_type none
mpiexec -n 54 ./diffusion_1d_mg settings_scaling.py 5400 gmres none gmres -ksp_type gmres -pc_type none
mpiexec -n 60 ./diffusion_1d_mg settings_scaling.py 6000 gmres none gmres -ksp_type gmres -pc_type none
mpiexec -n 66 ./diffusion_1d_mg settings_scaling.py 6600 gmres none gmres -ksp_type gmres -pc_type none
mpiexec -n 72 ./diffusion_1d_mg settings_scaling.py 7200 gmres none gmres -ksp_type gmres -pc_type none
