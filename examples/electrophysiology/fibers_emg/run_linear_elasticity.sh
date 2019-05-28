
/data/scratch/maierbn/opendihu/dependencies/mpi/install/bin/mpirun -n 4 \
./fibers_linear_elasticity ../settings_fibers_emg.py \
--n_subdomains 2 2 1 \
--linear_elasticity \
--fiber_file="../../input/7x7fibers.bin" \
--paraview_output \
--dt_0D=5e-4 \
--dt_1D=1e-3 \
--dt_3D=1e-3 \
--dt_bidomain=1e-2

