#perf stat -e r5301c7,r5304c7,r5310c7,r5340c7 -o perf.txt --delay 3000
 ./fast_fibers_emg ../settings_fibers_emg.py review_bw_stiftung.py \
--n_subdomains 1 1 1 \
--fiber_file="../../input/2x2fibers.bin" \
--firing_times_file="../../input/MU_firing_times_immediately.txt" \
--paraview_output \
--dt_0D=2e-3 \
--dt_1D=4e-3 \
--dt_splitting=4e-3 \
--dt_3D=1e-1 \
--output_timestep=1e-1 \
--end_time=10.0 
#-vmodule=vector_operators*=1

# double: 0x5310c7
# simd:  r5301c7

# perf stat -e r5301c7,r5304c7,r5310c7,r5340c7  -I 1000 --per-core



#cellml_file = "../../input/hodgkin_huxley_1952.c"
