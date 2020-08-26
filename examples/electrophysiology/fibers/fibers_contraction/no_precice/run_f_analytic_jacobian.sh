mpirun -n 24 ./biceps_contraction ../settings_biceps_contraction.py 15mus.py \
  --use_vc \
  --fast_monodomain_solver_optimizations \
  --use_analytic_jacobian \
  --end_time 100.0 \
  --scenario_name f



