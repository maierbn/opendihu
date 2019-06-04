# MG weak scaling
# args: <n_elements> <solver> <preconditioner> <name> <petsc options>


#CG-Boomeramg settigns testing

# threshold
mpiexec -n 1 ./diffusion_1d_mg settings_scaling.py 10000 cg pchypre cg_boomearmg -pc_type hypre -pc_hypre_type boomeramg

#coarsentype

#cycle_type

#


#Boomeramg settings testing

#threshold
mpiexec -n 1 ./diffusion_1d_mg settings_scaling.py 10000 richardson pchypre boomeramg -pc_type hypre -pc_hypre_type boomeramg

#coarsentype

#cycle_type

#


#CG-gamg settings testing

#threshold
mpiexec -n 1 ./diffusion_1d_mg settings_scaling.py 10000 cg gamg cg_gamg

#squaregraph

#symmetry

#cycle_type

#



#gamg settings testing

#threshold
mpiexec -n 1 ./diffusion_1d_mg settings_scaling.py 10000 preonly gamg gamg

#squaregraph

#symmetry

#cycle_type

#


