# MG strong scaling
# args: <n_elements> <solver> <preconditioner> <name> <petsc options>

mpiexec -n 1 ./diffusion_1d_mg settings_scaling.py 10000 gmres pchypre gmres_boomearmg -ksp_type gmres -pc_type hypre -pc_hypre_type boomeramg
mpiexec -n 2 ./diffusion_1d_mg settings_scaling.py 10000 gmres pchypre gmres_boomearmg -ksp_type gmres -pc_type hypre -pc_hypre_type boomeramg
mpiexec -n 4 ./diffusion_1d_mg settings_scaling.py 10000 gmres pchypre gmres_boomearmg -ksp_type gmres -pc_type hypre -pc_hypre_type boomeramg
mpiexec -n 8 ./diffusion_1d_mg settings_scaling.py 10000 gmres pchypre gmres_boomearmg -ksp_type gmres -pc_type hypre -pc_hypre_type boomeramg
mpiexec -n 16 ./diffusion_1d_mg settings_scaling.py 10000 gmres pchypre gmres_boomearmg -ksp_type gmres -pc_type hypre -pc_hypre_type boomeramg
mpiexec -n 32 ./diffusion_1d_mg settings_scaling.py 10000 gmres pchypre gmres_boomearmg -ksp_type gmres -pc_type hypre -pc_hypre_type boomeramg
mpiexec -n 64 ./diffusion_1d_mg settings_scaling.py 10000 gmres pchypre gmres_boomearmg -ksp_type gmres -pc_type hypre -pc_hypre_type boomeramg

###########################################################

mpiexec -n 1 ./diffusion_1d_mg settings_scaling.py 10000 richardson pchypre richardson_boomeramg -ksp_type richardson -pc_type hypre -pc_hypre_type boomeramg
mpiexec -n 2 ./diffusion_1d_mg settings_scaling.py 10000 richardson pchypre richardson_boomeramg -ksp_type richardson -pc_type hypre -pc_hypre_type boomeramg
mpiexec -n 4 ./diffusion_1d_mg settings_scaling.py 10000 richardson pchypre richardson_boomeramg -ksp_type richardson -pc_type hypre -pc_hypre_type boomeramg
mpiexec -n 8 ./diffusion_1d_mg settings_scaling.py 10000 richardson pchypre richardson_boomeramg -ksp_type richardson -pc_type hypre -pc_hypre_type boomeramg
mpiexec -n 16 ./diffusion_1d_mg settings_scaling.py 10000 richardson pchypre richardson_boomeramg -ksp_type richardson -pc_type hypre -pc_hypre_type boomeramg
mpiexec -n 32 ./diffusion_1d_mg settings_scaling.py 10000 richardson pchypre richardson_boomeramg -ksp_type richardson -pc_type hypre -pc_hypre_type boomeramg
mpiexec -n 64 ./diffusion_1d_mg settings_scaling.py 10000 richardson pchypre richardson_boomeramg -ksp_type richardson -pc_type hypre -pc_hypre_type boomeramg

############################################################

mpiexec -n 1 ./diffusion_1d_mg settings_scaling.py 10000 fgmres gamg fgmres_gamg -ksp_type fgmres -pc_type gamg
mpiexec -n 2 ./diffusion_1d_mg settings_scaling.py 10000 fgmres gamg fgmres_gamg -ksp_type fgmres -pc_type gamg
mpiexec -n 4 ./diffusion_1d_mg settings_scaling.py 10000 fgmres gamg fgmres_gamg -ksp_type fgmres -pc_type gamg
mpiexec -n 8 ./diffusion_1d_mg settings_scaling.py 10000 fgmres gamg fgmres_gamg -ksp_type fgmres -pc_type gamg
mpiexec -n 16 ./diffusion_1d_mg settings_scaling.py 10000 fgmres gamg fgmres_gamg -ksp_type fgmres -pc_type gamg
mpiexec -n 32 ./diffusion_1d_mg settings_scaling.py 10000 fgmres gamg fgmres_gamg -ksp_type fgmres -pc_type gamg
mpiexec -n 64 ./diffusion_1d_mg settings_scaling.py 10000 fgmres gamg fgmres_gamg -ksp_type fgmres -pc_type gamg

#############################################################

mpiexec -n 1 ./diffusion_1d_mg settings_scaling.py 10000 preonly gamg preonly_gamg -ksp_type preonly -pc_type gamg
mpiexec -n 2 ./diffusion_1d_mg settings_scaling.py 10000 preonly gamg preonly_gamg -ksp_type preonly -pc_type gamg
mpiexec -n 4 ./diffusion_1d_mg settings_scaling.py 10000 preonly gamg preonly_gamg -ksp_type preonly -pc_type gamg
mpiexec -n 8 ./diffusion_1d_mg settings_scaling.py 10000 preonly gamg preonly_gamg -ksp_type preonly -pc_type gamg
mpiexec -n 16 ./diffusion_1d_mg settings_scaling.py 10000 preonly gamg preonly_gamg -ksp_type preonly -pc_type gamg
mpiexec -n 32 ./diffusion_1d_mg settings_scaling.py 10000 preonly gamg preonly_gamg -ksp_type preonly -pc_type gamg
mpiexec -n 64 ./diffusion_1d_mg settings_scaling.py 10000 preonly gamg preonly_gamg -ksp_type preonly -pc_type gamg

mpiexec -n 1 ./diffusion_1d_mg settings_scaling.py 10000 gmres none gmres -ksp_type gmres -pc_type none
mpiexec -n 2 ./diffusion_1d_mg settings_scaling.py 10000 gmres none gmres -ksp_type gmres -pc_type none
mpiexec -n 4 ./diffusion_1d_mg settings_scaling.py 10000 gmres none gmres -ksp_type gmres -pc_type none
mpiexec -n 8 ./diffusion_1d_mg settings_scaling.py 10000 gmres none gmres -ksp_type gmres -pc_type none
mpiexec -n 16 ./diffusion_1d_mg settings_scaling.py 10000 gmres none gmres -ksp_type gmres -pc_type none
mpiexec -n 32 ./diffusion_1d_mg settings_scaling.py 10000 gmres none gmres -ksp_type gmres -pc_type none
mpiexec -n 64 ./diffusion_1d_mg settings_scaling.py 10000 gmres none gmres -ksp_type gmres -pc_type none
