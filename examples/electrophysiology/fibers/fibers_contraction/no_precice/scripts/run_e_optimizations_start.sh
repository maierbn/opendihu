mpirun -n 12 ./biceps_contraction ../settings_biceps_contraction.py 15mus.py \
  --use_vc \
  --fast_monodomain_solver_optimizations \
  --end_time 1.0 \
  --scenario_name e


#  --use_analytic_jacobian \

