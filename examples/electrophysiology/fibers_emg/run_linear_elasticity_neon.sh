
nice -n 10 mpirun -n 144 \
./fibers_linear_elasticity ../settings_fibers_emg.py \
--n_subdomains 6 6 4 \
--use_elasticity \
--fiber_file="../../input/13x13fibers.bin" \
--paraview_output \
--dt_0D=2e-4 \
--dt_1D=2e-3 \
--dt_splitting=2e-3 \
--dt_3D=1e-2 \
--output_timestep=2e0 \
--end_time=10000 
#-vmodule=mapping_between_meshes.tpp=1,*_n_com*=1 
#dt_0D=1e-4
#dt_1D=2e-3
#dt_splitting=2e-3
