
# scenario name for log file
scenario_name = "measure_flops"

# multiple fibers, new_slow_TK_2014_12_08.c this is the same as cuboid example with ModelType 0 (3a, "MultiPhysStrain", old tomo mechanics)

Conductivity = 3.828                # sigma, conductivity [mS/cm]
Am = 500.0                          # surface area to volume ratio [cm^-1]
Cm = 0.58                           # [uF/cm^2] membrane capacitance, (1 = fast twitch, 0.58 = slow twitch)
# diffusion prefactor = Conductivity/(Am*Cm)

# timing parameters
# -----------------
end_time = 10.0                     # [ms] end time of the simulation
stimulation_frequency = 100*1e-3    # [ms^-1] sampling frequency of stimuli in firing_times_file, in stimulations per ms, number before 1e-3 factor is in Hertz.
dt_0D = 1e-3                        # [ms] timestep width of ODEs
dt_1D = 1e-3                        # [ms] timestep width of diffusion
dt_splitting = dt_1D                 # [ms] overall timestep width of strang splitting
output_timestep = 0.1               # [ms] timestep for output files
activation_start_time = 1           # [ms] time when to start checking for stimulation

# other options
paraview_output = True

# input files
fiber_file              = "../../../input/13x13fibers.bin"
firing_times_file       = "../../../input/MU_firing_times_always.txt"
fiber_distribution_file = "../../../input/MU_fibre_distribution_3780.txt"
cellml_file             = "../../../input/new_slow_TK_2014_12_08.c"
