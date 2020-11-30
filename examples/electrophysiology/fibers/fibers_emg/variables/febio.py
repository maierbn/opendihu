
# scenario name for log file
scenario_name = "febio"

# material parameters
# --------------------

# Fixed units in cellMl models:
# These define the unit system.
# 1 cm = 1e-2 m
# 1 ms = 1e-3 s
# 1 uA = 1e-6 A
# 1 uF = 1e-6 F
# 
# derived units:
#   (F=s^4*A^2*m^-2*kg^-1) => 1 ms^4*uA^2*cm^-2*x*kg^-1 = (1e-3)^4 s^4 * (1e-6)^2 A^2 * (1e-2)^-2 m^-2 * (x)^-1 kg^-1 = 1e-12 * 1e-12 * 1e4 F = 1e-20 * x^-1 F := 1e-6 F => x = 1e-14
# 1e-14 kg = 10e-15 kg = 10e-12 g = 10 pg

# (N=kg*m*s^-2) => 1 10pg*cm*ms^2 = 1e-14 kg * 1e-2 m * (1e-3)^-2 s^-2 = 1e-14 * 1e-2 * 1e6 N = 1e-10 N = 10 nN
# (S=kg^-1*m^-2*s^3*A^2, Siemens not Sievert!) => (1e-14*kg)^-1*cm^-2*ms^3*uA^2 = (1e-14)^-1 kg^-1 * (1e-2)^-2 m^-2 * (1e-3)^3 s^3 * (1e-6)^2 A^2 = 1e14 * 1e4 * 1e-9 * 1e-12 S = 1e-3 S = 1 mS
# (V=kg*m^2*s^-3*A^-1) => 1 10pg*cm^2*ms^-3*uA^-1 = (1e-14) kg * (1e-2)^2 m^2 * (1e-3)^-3 s^-3 * (1e-6)^-1 A^-1 = 1e-14 * 1e-4 * 1e6 * 1e6 V = 1e-6 V = 1mV
# Hz=s^-1 => 1 ms^-1 = (1e-3)^-1 s^-1 = 1e3 Hz

# Hodgkin-Huxley
# t: ms
# STATES[0], Vm: mV
# CONSTANTS[1], Cm: uF*cm^-2
# CONSTANTS[2], I_Stim: uA*cm^-2
# -> all units are consistent

# Shorten
# t: ms
# CONSTANTS[0], Cm: uF*cm^-2
# STATES[0], Vm: mV
# ALGEBRAIC[32], I_Stim: uA*cm^-2
# -> all units are consistent

Conductivity = 3.828                # sigma, conductivity [mS/cm]
Am = 500.0                          # surface area to volume ratio [cm^-1]
Cm = 0.58                           # [uF/cm^2] membrane capacitance, (1 = fast twitch, 0.58 = slow twitch)
# diffusion prefactor = Conductivity/(Am*Cm)

# timing parameters
# -----------------
end_time = 10000.0                     # [ms] end time of the simulation
stimulation_frequency = 100*1e-3    # [ms^-1] sampling frequency of stimuli in firing_times_file, in stimulations per ms, number before 1e-3 factor is in Hertz.
dt_0D = 1e-3                        # [ms] timestep width of ODEs
dt_1D = 1e-3                        # [ms] timestep width of diffusion
dt_splitting = dt_1D                # [ms] overall timestep width of strang splitting
dt_3D = 1e-1                        # [ms] timestep width of elasticity problem
output_timestep = 1e0              # [ms] timestep for output files
activation_start_time = 0           # [ms] time when to start checking for stimulation

# other options
paraview_output = True
use_elasticity = True
generate_linear_3d_mesh = True
generate_quadratic_3d_mesh = False
states_output = False                # output states of the 0D problem in a separate file

# how to sample the 3D mesh for contraction from the fiber meshes, higher values mean coarser meshes
sampling_stride_x = 10
sampling_stride_y = 10
sampling_stride_z = 100

#n_fibers_x = 1
#n_fibers_y = 1
n_points_whole_fiber = 1480
    
# input files
fiber_file              = "../../../input/13x13fibers.bin"
firing_times_file       = "../../../input/MU_firing_times_always.txt"
fiber_distribution_file = "../../../input/MU_fibre_distribution_3780.txt"
cellml_file             = "../../../input/new_slow_TK_2014_12_08.c"

