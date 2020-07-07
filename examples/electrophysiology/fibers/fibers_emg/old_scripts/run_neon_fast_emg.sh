nice -n 10 mpirun -n 64 ./fast_fibers_emg ../settings_fibers_emg.py ramp_emg.py --n_subdomains 4 4 4 --emg_solver_type lu  --end_time 10000
