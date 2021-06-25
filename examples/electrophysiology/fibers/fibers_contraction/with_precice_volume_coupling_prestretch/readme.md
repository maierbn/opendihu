mpirun -n 64 ./fibers_with_3d_hh ../settings_muscle_contraction.py ramp.py
mpirun -n 2 ./muscle_contraction ../settings_muscle_contraction.py ramp.py
# e.g., on neon
