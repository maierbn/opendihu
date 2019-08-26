mpirun -n 40 \
  ./fibers_emg_2d_output ../settings_fibers_emg.py \
  --fiber_file ../../input/31x31fibers.bin \
  --n_subdomains 4 2 5 \
  --emg_initial_guess_nonzero \
  --end_time=1000.0 \
  --scenario_name=2d \
  --paraview_output 

