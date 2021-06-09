#./tendon_linear_precice_dynamic settings_tendon_top_a.py  | tee out_tendon_top_a.txt  &
#./tendon_linear_precice_dynamic settings_tendon_top_b.py  | tee out_tendon_top_b.txt  &
#./tendon_linear_precice_dynamic settings_tendon_bottom.py | tee out_tendon_bottom.txt &
mpirun -n 60 ./muscle_electrophysiology_precice settings_muscle.py hawk_medium.py --n_subdomains 2 2 15 | tee out_muscle.txt
