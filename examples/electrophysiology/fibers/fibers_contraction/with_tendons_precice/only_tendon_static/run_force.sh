# run tensile tests with different force
for force in `seq 0 200 2000`; do

  ./tendon_precice_quasistatic settings_only_tendon.py --fiber_file=../meshes/tendon_box.bin --force=$force

done
