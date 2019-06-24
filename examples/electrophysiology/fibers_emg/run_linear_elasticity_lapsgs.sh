./fibers_linear_elasticity ../settings_fibers_emg.py \
--n_subdomains 1 1 1 \
--linear_elasticity \
--fiber_file="../../input/2x2fibers.bin" \
--paraview_output \
--dt_0D=2e-3 \
--dt_1D=4e-3 \
--dt_splitting=4e-3 \
--dt_3D=1e-2 \
--output_timestep=1e-1 

