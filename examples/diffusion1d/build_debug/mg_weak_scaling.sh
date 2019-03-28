# MG weak scaling
# args: <n_elements> <solver> <preconditioner> <name> <petsc options>

mpiexec -n 1 ./diffusion_1d_mg settings_scaling.py 500 gmres pchypre gmres_boomearmg -ksp_type gmres -pc_type hypre -pc_hypre_type boomeramg
mpiexec -n 2 ./diffusion_1d_mg settings_scaling.py 1000 gmres pchypre gmres_boomearmg -ksp_type gmres -pc_type hypre -pc_hypre_type boomeramg
mpiexec -n 4 ./diffusion_1d_mg settings_scaling.py 2000 gmres pchypre gmres_boomearmg -ksp_type gmres -pc_type hypre -pc_hypre_type boomeramg
mpiexec -n 8 ./diffusion_1d_mg settings_scaling.py 4000 gmres pchypre gmres_boomearmg -ksp_type gmres -pc_type hypre -pc_hypre_type boomeramg
mpiexec -n 16 ./diffusion_1d_mg settings_scaling.py 8000 gmres pchypre gmres_boomearmg -ksp_type gmres -pc_type hypre -pc_hypre_type boomeramg
mpiexec -n 32 ./diffusion_1d_mg settings_scaling.py 16000 gmres pchypre gmres_boomearmg -ksp_type gmres -pc_type hypre -pc_hypre_type boomeramg
mpiexec -n 64 ./diffusion_1d_mg settings_scaling.py 32000 gmres pchypre gmres_boomearmg -ksp_type gmres -pc_type hypre -pc_hypre_type boomeramg

###########################################################

mpiexec -n 1 ./diffusion_1d_mg settings_scaling.py 500 richardson pchypre richardson_boomeramg -ksp_type richardson -pc_type hypre -pc_hypre_type boomeramg
mpiexec -n 2 ./diffusion_1d_mg settings_scaling.py 1000 richardson pchypre richardson_boomeramg -ksp_type richardson -pc_type hypre -pc_hypre_type boomeramg
mpiexec -n 4 ./diffusion_1d_mg settings_scaling.py 2000 richardson pchypre richardson_boomeramg -ksp_type richardson -pc_type hypre -pc_hypre_type boomeramg
mpiexec -n 8 ./diffusion_1d_mg settings_scaling.py 4000 richardson pchypre richardson_boomeramg -ksp_type richardson -pc_type hypre -pc_hypre_type boomeramg
mpiexec -n 16 ./diffusion_1d_mg settings_scaling.py 8000 richardson pchypre richardson_boomeramg -ksp_type richardson -pc_type hypre -pc_hypre_type boomeramg
mpiexec -n 32 ./diffusion_1d_mg settings_scaling.py 16000 richardson pchypre richardson_boomeramg -ksp_type richardson -pc_type hypre -pc_hypre_type boomeramg
mpiexec -n 64 ./diffusion_1d_mg settings_scaling.py 32000 richardson pchypre richardson_boomeramg -ksp_type richardson -pc_type hypre -pc_hypre_type boomeramg

############################################################

mpiexec -n 1 ./diffusion_1d_mg settings_scaling.py 500 fgmres gamg fgmres_gamg -ksp_type fgmres -pc_type gamg
mpiexec -n 2 ./diffusion_1d_mg settings_scaling.py 1000 fgmres gamg fgmres_gamg -ksp_type fgmres -pc_type gamg
mpiexec -n 4 ./diffusion_1d_mg settings_scaling.py 2000 fgmres gamg fgmres_gamg -ksp_type fgmres -pc_type gamg
mpiexec -n 8 ./diffusion_1d_mg settings_scaling.py 4000 fgmres gamg fgmres_gamg -ksp_type fgmres -pc_type gamg
mpiexec -n 16 ./diffusion_1d_mg settings_scaling.py 8000 fgmres gamg fgmres_gamg -ksp_type fgmres -pc_type gamg
mpiexec -n 32 ./diffusion_1d_mg settings_scaling.py 16000 fgmres gamg fgmres_gamg -ksp_type fgmres -pc_type gamg
mpiexec -n 64 ./diffusion_1d_mg settings_scaling.py 32000 fgmres gamg fgmres_gamg -ksp_type fgmres -pc_type gamg

#############################################################

mpiexec -n 1 ./diffusion_1d_mg settings_scaling.py 500 preonly gamg preonly_gamg -ksp_type preonly -pc_type gamg
mpiexec -n 2 ./diffusion_1d_mg settings_scaling.py 1000 preonly gamg preonly_gamg -ksp_type preonly -pc_type gamg
mpiexec -n 4 ./diffusion_1d_mg settings_scaling.py 2000 preonly gamg preonly_gamg -ksp_type preonly -pc_type gamg
mpiexec -n 8 ./diffusion_1d_mg settings_scaling.py 4000 preonly gamg preonly_gamg -ksp_type preonly -pc_type gamg
mpiexec -n 16 ./diffusion_1d_mg settings_scaling.py 8000 preonly gamg preonly_gamg -ksp_type preonly -pc_type gamg
mpiexec -n 32 ./diffusion_1d_mg settings_scaling.py 16000 preonly gamg preonly_gamg -ksp_type preonly -pc_type gamg
mpiexec -n 64 ./diffusion_1d_mg settings_scaling.py 32000 preonly gamg preonly_gamg -ksp_type preonly -pc_type gamg

#############################################################

mpiexec -n 1 ./diffusion_1d_mg settings_scaling.py 500 gmres none gmres -ksp_type gmres -pc_type none
mpiexec -n 2 ./diffusion_1d_mg settings_scaling.py 1000 gmres none gmres -ksp_type gmres -pc_type none
mpiexec -n 4 ./diffusion_1d_mg settings_scaling.py 2000 gmres none gmres -ksp_type gmres -pc_type none
mpiexec -n 8 ./diffusion_1d_mg settings_scaling.py 4000 gmres none gmres -ksp_type gmres -pc_type none
mpiexec -n 16 ./diffusion_1d_mg settings_scaling.py 8000 gmres none gmres -ksp_type gmres -pc_type none
mpiexec -n 32 ./diffusion_1d_mg settings_scaling.py 16000 gmres none gmres -ksp_type gmres -pc_type none
mpiexec -n 64 ./diffusion_1d_mg settings_scaling.py 32000 gmres none gmres -ksp_type gmres -pc_type none
