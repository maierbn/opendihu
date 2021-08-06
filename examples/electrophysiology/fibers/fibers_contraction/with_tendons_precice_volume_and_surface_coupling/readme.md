mpirun -n 64 ./fibers_with_3d_hh ../settings_fibers_with_3d.py ramp.py
mpirun -n 2 ./muscle ../settings_muscle_contraction.py ramp.py
./tendon ../settings_tendon_bottom.py 
./tendon ../settings_tendon_top_a.py 
./tendon ../settings_tendon_top_b.py 
