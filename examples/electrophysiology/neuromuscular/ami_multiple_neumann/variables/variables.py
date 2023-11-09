# scenario name for log file
case_name = "default"
scenario_name = ""
precice_file = "../precice_config.xml"

# timing parameters
# -----------------
dt_elasticity = 0.01
dt_0D = 0.5e-3                        # [ms] timestep width of ODEs
dt_1D = 1e-3                      # [ms] timestep width of diffusion
dt_splitting_0D1D = dt_1D           # [ms] overall timestep width of strang splitting
end_time = 20.0

output_interval_fibers = dt_elasticity/dt_1D
output_interval_0D = dt_elasticity/dt_0D
output_timestep_elasticity = 10      # [ms] timestep for elasticity output files



stimulation_frequency = 100*1e-3    # [ms^-1] sampling frequency of stimuli in firing_times_file, in stimulations per ms, number before 1e-3 factor is in Hertz. This is not used here.
activation_start_time = 0.0           # [ms] time when to start checking for stimulation

# partitioning
# ------------
# this has to match the total number of processes
n_subdomains_x = 1
n_subdomains_y = 1
n_subdomains_z = 1

tendon_material = "nonLinear"
#tendon_material = "SaintVenantKirchoff"         

tendon_extent = [3.0, 3.0, 5.0] # [cm, cm, cm] 2.96
n_elements_tendon = [8,8,8] 

muscle_extent = [3.0, 3.0, 15.0] # [cm, cm, cm]
n_elements_muscle = [8, 8, 24] # linear elements. each qudaratic element uses the combined nodes of 8 linear elements
n_points_whole_fiber = 60
n_fibers_x = 8
n_fibers_y = 8

muscle_left_offset = [0.0, 0.0, 0.0]
tendon_offset = [0.0, 0.0, muscle_extent[2]]
muscle_right_offset = [0.0, 0.0, tendon_extent[2]+muscle_extent[2]]

rho = 10   ## [1e-4 kg/cm^3] density of the water

force = 1e5

elasticity_dirichlet_bc = {}
elasticity_neumann_bc = []
meshes = {}
# material parameters
# --------------------
Pmax = 7.3                          # maximum stress [N/cm^2]
Conductivity = 3.828                # sigma, conductivity [mS/cm]
Am = 500.0                          # surface area to volume ratio [cm^-1]
Cm = 0.58                           # membrane capacitance [uF/cm^2]
damping_factor = 0                  # velocity dependent damping factor

rho = 10                    # [1e-4 kg/cm^3] density of the muscle (density of water)
# Mooney-Rivlin parameters [c1,c2,b,d] of c1*(Ibar1 - 3) + c2*(Ibar2 - 3) + b/d (λ - 1) - b*ln(λ)
# Heidlauf13: [6.352e-10 kPa, 3.627 kPa, 2.756e-5 kPa, 43.373] = [6.352e-11 N/cm^2, 3.627e-1 N/cm^2, 2.756e-6 N/cm^2, 43.373], pmax = 73 kPa = 7.3 N/cm^2
# Heidlauf16: [3.176e-10 N/cm^2, 1.813 N/cm^2, 1.075e-2 N/cm^2, 9.1733], pmax = 7.3 N/cm^2
c1 = 3.176e-10              # [N/cm^2]
c2 = 1.813                  # [N/cm^2]
b  = 1.075e-2               # [N/cm^2] anisotropy parameter
d  = 9.1733                 # [-] anisotropy parameter
# for debugging, b = 0 leads to normal Mooney-Rivlin

muscle_material_parameters = [c1, c2, b, d]   # material parameters

import os
input_directory = os.path.join(os.environ.get('OPENDIHU_HOME', '../../../../../'), "examples/electrophysiology/input")

cellml_file = input_directory+"/hodgkin_huxley-razumova.cellml"
fiber_distribution_file = input_directory+"/MU_fibre_distribution_multidomain_67x67_100.txt"
firing_times_file = input_directory + "/MU_firing_times_real.txt"
no_firing_times_file = input_directory + "/MU_firing_times_never.txt" # no firing

maximum_number_of_threads = 1
use_aovs_memory_layout = True
fast_monodomain_solver_optimizations = True

own_subdomain_coordinate_x = 0 # TODO fix this for parallelization
own_subdomain_coordinate_y = 0 # TODO fix this for parallelization
own_subdomain_coordinate_z = 0 # TODO fix this for parallelization

states_output = False

# functions, here, Am, Cm and Conductivity are constant for all fibers and MU's
# These functions can be redefined differently in a custom variables script
def get_am(fiber_no, mu_no):
  return Am

def get_cm(fiber_no, mu_no):
  return Cm
  
def get_conductivity(fiber_no, mu_no):
  return Conductivity

def get_specific_states_call_frequency(fiber_no, mu_no):
  return stimulation_frequency

def get_specific_states_frequency_jitter(fiber_no, mu_no):
  return [0]

def get_specific_states_call_enable_begin(fiber_no, mu_no):
  return activation_start_time

# solvers
# -------
diffusion_solver_type = "cg"        # solver and preconditioner for the diffusion part of the Monodomain equation
diffusion_preconditioner_type = "none"      # preconditioner
diffusion_solver_maxit = 1e4
diffusion_solver_reltol = 1e-10

elasticity_solver_type = "preonly"
elasticity_preconditioner_type = "lu"
snes_max_iterations = 10                  # maximum number of iterations in the nonlinear solver
snes_rebuild_jacobian_frequency = 2       # how often the jacobian should be recomputed, -1 indicates NEVER rebuild, 1 means rebuild every time the Jacobian is computed within a single nonlinear solve, 2 means every second time the Jacobian is built etc. -2 means rebuild at next chance but then never again 
snes_relative_tolerance = 1e-5      # relative tolerance of the nonlinear solver
snes_absolute_tolerance = 1e-5      # absolute tolerance of the nonlinear solver
linear_relative_tolerance = 1e-5           # relative tolerance of the residual of the linear solver
linear_absolute_tolerance = 1e-10          # absolute tolerance of the residual of the linear solver

optimization_type = "vc"            # the optimization_type used in the cellml adapter, "vc" uses explicit vectorization
approximate_exponential_function = False   # if the exponential function should be approximated by a Taylor series with only 11 FLOPS
dynamic = True                      # if the dynamic hyperelasticity solver should be used

def get_from_obj(data, path):
    for elem in path:
        if type(elem) == str:
            data = data[elem]
        elif type(elem) == int:
            data = data[elem]
        elif type(elem) == tuple:
            # search for key == value with (key, value) = elem
            key, value = elem
            data = next(filter(lambda e: e[key] == value, data))
        else:
            raise KeyError(f"Unknown type of '{elem}': '{type(elem)}'. Path: '{'.'.join(path)}'")
    return data

# available slots for data fields

# geometry
# u
# v
# PK2-Stress (Voigt)
# active PK2-Stress (Voigt)
# fiberDirection
# t (current traction)
# T (material traction)
# F
# Fdot
# P (PK1 stress)
# σ (Cauchy stress)
# J

def muscle_left_write_to_file(data):
    t = get_from_obj(data, [0, 'currentTime'])
    z_data = get_from_obj(data, [0, 'data', ('name','geometry'), 'components', 2, 'values'])

    [mx, my, mz] = get_from_obj(data, [0, 'nElementsLocal'])
    nx = 2*mx + 1
    ny = 2*my + 1
    nz = 2*mz + 1
    # compute average z-value of end of muscle
    z_value_begin = 0.0
    z_value_end = 0.0

    for j in range(ny):
        for i in range(nx):
            z_value_begin += z_data[j*nx + i]
            z_value_end += z_data[(nz-1)*nx*ny + j*nx + i]


    z_value_begin /= ny*nx
    z_value_end /= ny*nx


    f = open(case_name + "muscle_left.txt", "a")
    f.write("{:6.2f} {:+2.8f} {:+2.8f}\n".format(t,z_value_begin, z_value_end))
    f.close()

def muscle_right_write_to_file(data):
    t = get_from_obj(data, [0, 'currentTime'])
    z_data = get_from_obj(data, [0, 'data', ('name','geometry'), 'components', 2, 'values'])

    [mx, my, mz] = get_from_obj(data, [0, 'nElementsLocal'])
    nx = 2*mx + 1
    ny = 2*my + 1
    nz = 2*mz + 1
    # compute average z-value of end of muscle
    z_value_begin = 0.0
    z_value_end = 0.0

    for j in range(ny):
        for i in range(nx):
            z_value_begin += z_data[j*nx + i]
            z_value_end += z_data[(nz-1)*nx*ny + j*nx + i]

    z_value_begin /= ny*nx
    z_value_end /= ny*nx

    f = open(case_name + "muscle_right.txt", "a")
    f.write("{:6.2f} {:+2.8f} {:+2.8f}\n".format(t,z_value_begin, z_value_end))
    f.close()

def tendon_write_to_file(data):
    t = get_from_obj(data, [0, 'currentTime'])
    z_data = get_from_obj(data, [0, 'data', ('name','geometry'), 'components', 2, 'values'])

    [mx, my, mz] = get_from_obj(data, [0, 'nElementsLocal'])
    nx = 2*mx + 1
    ny = 2*my + 1
    nz = 2*mz + 1
    # compute average z-value of end of muscle
    z_value_begin = 0.0
    z_value_end = 0.0

    for j in range(ny):
        for i in range(nx):
            z_value_begin += z_data[j*nx + i]
            z_value_end += z_data[(nz-1)*nx*ny + j*nx + i]

    z_value_begin /= ny*nx
    z_value_end /= ny*nx


    f = open(case_name + "tendon.txt", "a")
    f.write("{:6.2f} {:+2.8f} {:+2.8f}\n".format(t,z_value_begin, z_value_end))
    f.close()

