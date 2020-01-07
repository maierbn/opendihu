./fibers_linear_elasticity ../settings_fibers_emg.py \
--n_subdomains 1 1 1 \
--use_elasticity \
--fiber_file="../../input/2x2fibers.bin" \
--paraview_output \
--dt_0D=2e-4 \
--dt_1D=2e-3 \
--dt_splitting=2e-3 \
--dt_3D=1e-2 \
--output_timestep=1e-1 
#-vmodule=heun.tpp=1 
#-vmodule=09_function_space_structured_check_nei*=1 

