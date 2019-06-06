mpirun -n 2 \
./fibers_linear_elasticity ../settings_fibers_emg.py \
--n_subdomains 2 1 1 \
--linear_elasticity \
--fiber_file="../../input/7x7fibers.bin" \
--paraview_output \
--dt_0D=2e-3 \
--dt_1D=4e-3 \
--dt_3D=4e-3 \
--dt_bidomain=1e-2 \
--output_timestep=1e-1 

