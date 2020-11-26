
# scenario name for log file
scenario_name = "biceps"

# timing parameters
# -----------------
end_time = 100.0                    # [ms] end time of the simulation
stimulation_frequency = 100*1e-3    # [ms^-1] sampling frequency of stimuli in firing_times_file, in stimulations per ms, number before 1e-3 factor is in Hertz.
dt_3D = 1e-1                           # [ms] time step width of coupling, when 3D should be performed, also sampling time of monopolar EMG
output_timestep = 1                 # [ms] timestep for output surface EMG, 0.5
output_timestep_fibers = 1          # [ms] timestep for fiber output, 0.5
output_timestep_big = 1             # [ms] timestep for output files of 3D intramuscular EMG data

# 3D mesh
# stride for sampling the 3D mesh from the fiber_file data, higher values mean a coarser EMG mesh
sampling_stride_x = 1
sampling_stride_y = 1
sampling_stride_z = 50

# other options
paraview_output = True
adios_output = False
exfile_output = False
python_output = False
enable_surface_emg = True
disable_firing_output = False

#fiber_file = "../../../input/left_biceps_brachii_13x13fibers.bin"
fiber_file = "../../../input/left_biceps_brachii_7x7fibers.bin"
firing_times_file = "../../../input/MU_firing_times_real.txt"
fiber_distribution_file = "../../../input/MU_fibre_distribution_10MUs.txt"
