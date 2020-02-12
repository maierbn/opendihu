
nice -n 10 mpirun -n 64 \
./fast_fibers_emg ../settings_fibers_emg.py \
--scenario_name fast_fibers \
--n_subdomains 4 4 4 \
--fiber_file="../../input/13x13fibers.bin" \
--paraview_output \
--emg_solver_type lu \
--dt_0D=2e-3 \
--dt_1D=4e-3 \
--dt_splitting=4e-3 \
--dt_3D=1e-3 \
--output_timestep=2e0 \
--end_time=100 

#--emg_initial_guess_nonzero \


#-vmodule=mapping_between_meshes.tpp=1,*_n_com*=1 
#dt_0D=1e-4
#dt_1D=2e-3
#dt_splitting=2e-3
