# Multiple 1D fibers (monodomain) with 3D EMG (static bidomain), biceps geometry
# arguments: -help
#
# if fiber_file=cuboid.bin, it uses a small cuboid test example

end_time = 100.0

import numpy as np
import pickle
import sys
import struct
import argparse

# global parameters
PMax = 7.3              # maximum stress [N/cm^2]
Conductivity = 3.828    # sigma, conductivity [mS/cm]
Am = 500.0              # surface area to volume ratio [cm^-1]
Cm = 0.58               # membrane capacitance [uF/cm^2]
innervation_zone_width = 1.  # cm
innervation_zone_width = 0.  # cm
diffusion_solver_type = "cg"
diffusion_preconditioner_type = "none"
potential_flow_solver_type = "gmres"
potential_flow_preconditioner_type = "none"
emg_solver_type = "gmres"
emg_preconditioner_type = "none"
emg_initial_guess_nonzero = False

# timing parameters
stimulation_frequency = 10000*1e-3   # sampling frequency of stimuli in firing_times_file, in stimulations per ms, number before 1e-3 factor is in Hertz.
dt_1D = 1e-3                      # timestep width of diffusion
dt_0D = 1.5e-3                    # timestep width of ODEs
dt_splitting = 3e-3               # overall timestep width of strang splitting
dt_3D = 1e0                 # time step width of coupling, when 3D should be performed
output_timestep = 1e0             # timestep for output files

# input files
#cellml_file = "../../input/shorten_ocallaghan_davidson_soboleva_2007.c"
#cellml_file = "../../input/shorten.cpp"
cellml_file = "../../input/hodgkin_huxley_1952.c"

#fiber_file = "../../input/3000fibers.bin"
#fiber_file = "../../input/7x7fibers.bin"
fiber_file = "../../input/13x13fibers.bin"
#fiber_file = "../../input/49fibers.bin"
load_data_from_file = False
debug_output = False
disable_firing_output = False

fiber_distribution_file = "../../input/MU_fibre_distribution_3780.txt"
firing_times_file = "../../input/MU_firing_times_real.txt"
#firing_times_file = "../../input/MU_firing_times_immediately.txt"

#print("prefactor: ",Conductivity/(Am*Cm))
#print("numpy path: ",np.__path__)

# partitioning
# this has to match the total number of processes
n_subdomains_x = 1   # example values for 4 processes
n_subdomains_y = 1
n_subdomains_z = 1

# stride for sampling the 3D elements from the fiber data
# here any number is possible
sampling_stride_x = 2
sampling_stride_y = 2
sampling_stride_z = 50

# parse arguments
scenario_name = ""
rank_no = (int)(sys.argv[-2])
n_ranks = (int)(sys.argv[-1])

# define command line arguments
parser = argparse.ArgumentParser(description='fibers_emg')
parser.add_argument('--scenario_name',            help='The name to identify this run in the log.', default=scenario_name)
parser.add_argument('--n_subdomains', nargs=3,    type=int, help='Number of subdomains in x,y,z direction.')
parser.add_argument('--n_subdomains_x', '-x',     type=int, help='Number of subdomains in x direction.', default=n_subdomains_x)
parser.add_argument('--n_subdomains_y', '-y',     type=int, help='Number of subdomains in y direction.', default=n_subdomains_y)
parser.add_argument('--n_subdomains_z', '-z',     type=int, help='Number of subdomains in z direction.', default=n_subdomains_z)
parser.add_argument('--diffusion_solver_type',    help='The solver for the diffusion.', default=diffusion_solver_type, choices=["gmres","cg","lu","gamg","richardson","chebyshev","cholesky","jacobi","sor","preonly"])
parser.add_argument('--diffusion_preconditioner_type',    help='The preconditioner for the diffusion.', default=diffusion_preconditioner_type, choices=["jacobi","sor","lu","ilu","gamg","none"])
parser.add_argument('--potential_flow_solver_type', help='The solver for the potential flow (non-spd matrix).', default=potential_flow_solver_type, choices=["gmres","cg","lu","gamg","richardson","chebyshev","cholesky","jacobi","sor","preonly"])
parser.add_argument('--potential_flow_preconditioner_type',    help='The preconditioner for the potential flow.', default=potential_flow_preconditioner_type, choices=["jacobi","sor","lu","ilu","gamg","none"])
parser.add_argument('--emg_solver_type',          help='The solver for the static bidomain.', default=emg_solver_type)
#parser.add_argument('--emg_solver_type',          help='The solver for the static bidomain.', default=emg_solver_type, choices=["gmres","cg","lu","gamg","richardson","chebyshev","cholesky","jacobi","sor","preonly"])
parser.add_argument('--emg_preconditioner_type',  help='The preconditioner for the static bidomain.', default=emg_preconditioner_type, choices=["jacobi","sor","lu","ilu","gamg","none"])
parser.add_argument('--emg_initial_guess_nonzero',  help='It the initial guess for the emg linear system should be set to the previous solution.', default=False, action='store_true')
parser.add_argument('--paraview_output',          help='Enable the paraview output writer.', default=False, action='store_true')
parser.add_argument('--fiber_file',               help='The filename of the file that contains the fiber data.', default=fiber_file)
parser.add_argument('--cellml_file',              help='The filename of the file that contains the cellml model.', default=cellml_file)
parser.add_argument('--fiber_distribution_file',  help='The filename of the file that contains the MU firing times.', default=fiber_distribution_file)
parser.add_argument('--firing_times_file',        help='The filename of the file that contains the cellml model.', default=firing_times_file)
parser.add_argument('--end_time', '--tend', '-t', type=float, help='The end simulation time.', default=end_time)
parser.add_argument('--output_timestep',          type=float, help='The timestep for writing outputs.', default=output_timestep)
parser.add_argument('--dt_0D',                    type=float, help='The timestep for the 0D model.', default=dt_0D)
parser.add_argument('--dt_1D',                    type=float, help='The timestep for the 1D model.', default=dt_1D)
parser.add_argument('--dt_splitting',             type=float, help='The timestep for the splitting.', default=dt_splitting)
parser.add_argument('--dt_3D',                    type=float, help='The timestep for the 3D model, either bidomain or mechanics.', default=dt_3D)
parser.add_argument('--disable_firing_output',    help='Disables the initial list of fiber firings.', default=False, action='store_true')
parser.add_argument('--v',                        help='Enable full verbosity in c++ code')
parser.add_argument('-v',                         help='Enable verbosity level in c++ code', action="store_true")
parser.add_argument('-vmodule',                   help='Enable verbosity level for given file in c++ code')
parser.add_argument('--rank_reordering',          help='Enable rank reordering in the c++ code', action="store_true")
parser.add_argument('--linear_elasticity',        help='Enable linear elasticity', action="store_true")
 
# parse arguments and assign values to global variables
args = parser.parse_args(args=sys.argv[:-2])
globals().update(args.__dict__)
if n_subdomains is not None:
  n_subdomains_x = n_subdomains[0]
  n_subdomains_y = n_subdomains[1]
  n_subdomains_z = n_subdomains[2]
  
if linear_elasticity:
  cellml_file = "../../input/shorten.cpp"
  emg_solver_type = "cg"
  
if rank_no == 0:
  print("scenario_name: {},  n_subdomains: {} {} {},  n_ranks: {},  end_time: {}".format(scenario_name, n_subdomains_x, n_subdomains_y, n_subdomains_z, n_ranks, end_time))
  print("dt_0D:           {:0.0e}, diffusion_solver_type:      {}".format(dt_0D, diffusion_solver_type))
  print("dt_1D:           {:0.0e}, potential_flow_solver_type: {}".format(dt_1D, potential_flow_solver_type))
  print("dt_splitting:    {:0.0e}, emg_solver_type:            {}, emg_initial_guess_nonzero: {}".format(dt_splitting, emg_solver_type, emg_initial_guess_nonzero))
  print("dt_3D:           {:0.0e}, paraview_output: {}".format(dt_3D, paraview_output))
  print("output_timestep: {:0.0e}  stimulation_frequency: {} 1/ms = {} Hz".format(output_timestep, stimulation_frequency, stimulation_frequency*1e3))
  print("fiber_file:              {}".format(fiber_file))
  print("cellml_file:             {}".format(cellml_file))
  print("fiber_distribution_file: {}".format(fiber_distribution_file))
  print("firing_times_file:       {}".format(firing_times_file))
  print("********************************************************************************")

if rank_no == 0:
  import timeit
  t_start_script = timeit.default_timer()
  
n_subdomains_xy = n_subdomains_x * n_subdomains_y
own_subdomain_coordinate_x = rank_no % n_subdomains_x
own_subdomain_coordinate_y = (int)(rank_no / n_subdomains_x) % n_subdomains_y
own_subdomain_coordinate_z = (int)(rank_no / n_subdomains_xy)

#print("rank: {}/{}".format(rank_no,n_ranks))

# generate cuboid fiber file
if fiber_file == "cuboid.bin":
  
  size_x = 1e-1
  size_y = 1e-1
  size_z = 0.2
  
  n_fibers_x = 4
  n_fibers_y = n_fibers_x
  n_points_whole_fiber = 20
  
  if rank_no == 0:
    print("create cuboid.bin with size [{},{},{}], n points [{},{},{}]".format(size_x, size_y, size_z, n_fibers_x, n_fibers_y, n_points_whole_fiber))
    
    # write header
    with open(fiber_file, "wb") as outfile:
      
      # write header
      header_str = "opendihu self-generated cuboid  "
      outfile.write(struct.pack('32s',bytes(header_str, 'utf-8')))   # 32 bytes
      outfile.write(struct.pack('i', 40))  # header length
      outfile.write(struct.pack('i', n_fibers_x*n_fibers_x))   # n_fibers
      outfile.write(struct.pack('i', n_points_whole_fiber))   # n_points_whole_fiber
      outfile.write(struct.pack('i', 0))   # nBorderPointsXNew
      outfile.write(struct.pack('i', 0))   # nBorderPointsZNew
      outfile.write(struct.pack('i', 0))   # nFineGridFibers_
      outfile.write(struct.pack('i', 1))   # nRanks
      outfile.write(struct.pack('i', 1))   # nRanksZ
      outfile.write(struct.pack('i', 0))   # nFibersPerRank
      outfile.write(struct.pack('i', 0))   # date
    
      # loop over points
      for y in range(n_fibers_x):
        for x in range(n_fibers_x):
          for z in range(n_points_whole_fiber):
            point = [x*(float)(size_x)/(n_fibers_x-1), y*(float)(size_y)/(n_fibers_y-1), z*(float)(size_z)/(n_points_whole_fiber-1)]
            outfile.write(struct.pack('3d', point[0], point[1], point[2]))   # data point

# set output writer    
output_writer_fibers = []
output_writer_emg = []
if paraview_output:
  output_writer_emg.append({"format": "Paraview", "outputInterval": int(1./dt_3D*output_timestep), "filename": "out/" + scenario_name + "/emg", "binary": True, "fixedFormat": False, "combineFiles": True})
  output_writer_fibers.append({"format": "Paraview", "outputInterval": int(1./dt_splitting*output_timestep), "filename": "out/" + scenario_name + "/fibers", "binary": True, "fixedFormat": False, "combineFiles": True})

# set values for cellml model
if "shorten" in cellml_file:
  parameters_used_as_intermediate = [32]
  parameters_used_as_constant = [65]
  parameters_initial_values = [0.0, 1.0]
  nodal_stimulation_current = 1200.
  
elif "hodgkin_huxley" in cellml_file:
  parameters_used_as_intermediate = []
  parameters_used_as_constant = [2]
  parameters_initial_values = [0.0]
  nodal_stimulation_current = 40.

def get_motor_unit_no(fiber_no):
  return int(fiber_distribution[fiber_no % len(fiber_distribution)]-1)

def fiber_gets_stimulated(fiber_no, frequency, current_time):
  """
  determine if fiber fiber_no gets stimulated at simulation time current_time
  """

  # determine motor unit
  alpha = 1.0   # 0.8
  mu_no = (int)(get_motor_unit_no(fiber_no)*alpha)
  
  # determine if fiber fires now
  index = int(np.round(current_time * frequency))
  n_firing_times = np.size(firing_times,0)
  
  #if firing_times[index % n_firing_times, mu_no] == 1:
    #print("{}: fiber {} is mu {}, t = {}, row: {}, stimulated: {} {}".format(rank_no, fiber_no, mu_no, current_time, (index % n_firing_times), firing_times[index % n_firing_times, mu_no], "true" if firing_times[index % n_firing_times, mu_no] == 1 else "false"))
  
  return firing_times[index % n_firing_times, mu_no] == 1
  
def set_parameters(n_nodes_global, time_step_no, current_time, parameters, dof_nos_global, fiber_no):
  
  # determine if fiber gets stimulated at the current time
  is_fiber_gets_stimulated = fiber_gets_stimulated(fiber_no, stimulation_frequency, current_time)
  
  # determine nodes to stimulate (center node, left and right neighbour)
  innervation_zone_width_n_nodes = innervation_zone_width*100  # 100 nodes per cm
  innervation_node_global = int(n_nodes_global / 2)  # + np.random.randint(-innervation_zone_width_n_nodes/2,innervation_zone_width_n_nodes/2+1)
  nodes_to_stimulate_global = [innervation_node_global]
  if innervation_node_global > 0:
    nodes_to_stimulate_global.insert(0, innervation_node_global-1)
  if innervation_node_global < n_nodes_global-1:
    nodes_to_stimulate_global.append(innervation_node_global+1)
  
  # stimulation value
  if is_fiber_gets_stimulated:
    stimulation_current = nodal_stimulation_current
  else:
    stimulation_current = 0.
  
  first_dof_global = dof_nos_global[0]
  last_dof_global = dof_nos_global[-1]
    
  for node_no_global in nodes_to_stimulate_global:
    if first_dof_global <= node_no_global <= last_dof_global:
      # get local no for global no (1D)
      dof_no_local = node_no_global - first_dof_global
      parameters[dof_no_local] = stimulation_current
 
# callback function that can set parameters, i.e. stimulation current
def set_specific_parameters(n_nodes_global, time_step_no, current_time, parameters, fiber_no):
  
  # determine if fiber gets stimulated at the current time
  is_fiber_gets_stimulated = fiber_gets_stimulated(fiber_no, stimulation_frequency, current_time)
  
  # determine nodes to stimulate (center node, left and right neighbour)
  innervation_zone_width_n_nodes = innervation_zone_width*100  # 100 nodes per cm
  innervation_node_global = int(n_nodes_global / 2)  # + np.random.randint(-innervation_zone_width_n_nodes/2,innervation_zone_width_n_nodes/2+1)
  nodes_to_stimulate_global = [innervation_node_global]
  
  for k in range(10):
    if innervation_node_global-k >= 0:
      nodes_to_stimulate_global.insert(0, innervation_node_global-k)
    if innervation_node_global+k <= n_nodes_global-1:
      nodes_to_stimulate_global.append(innervation_node_global+k)
  
  # stimulation value
  if is_fiber_gets_stimulated:
    stimulation_current = 40.
  else:
    stimulation_current = 0.

  for node_no_global in nodes_to_stimulate_global:
    parameters[(node_no_global,0)] = stimulation_current   # key: ((x,y,z),nodal_dof_index)

# callback function that can set states, i.e. prescribed values for stimulation
def set_specific_states(n_nodes_global, time_step_no, current_time, states, fiber_no):

  #print("call set_specific_states at time {}".format(current_time)) 
  
  # determine if fiber gets stimulated at the current time
  is_fiber_gets_stimulated = fiber_gets_stimulated(fiber_no, stimulation_frequency, current_time)

  if is_fiber_gets_stimulated:  
    # determine nodes to stimulate (center node, left and right neighbour)
    innervation_zone_width_n_nodes = innervation_zone_width*100  # 100 nodes per cm
    innervation_node_global = int(n_nodes_global / 2)  # + np.random.randint(-innervation_zone_width_n_nodes/2,innervation_zone_width_n_nodes/2+1)
    nodes_to_stimulate_global = [innervation_node_global]
    if innervation_node_global > 0:
      nodes_to_stimulate_global.insert(0, innervation_node_global-1)
    if innervation_node_global < n_nodes_global-1:
      nodes_to_stimulate_global.append(innervation_node_global+1)
    if rank_no == 0:
      print("t: {}, stimulate fiber {} at nodes {}".format(current_time, fiber_no, nodes_to_stimulate_global))

    for node_no_global in nodes_to_stimulate_global:
      states[(node_no_global,0,0)] = 20.0   # key: ((x,y,z),nodal_dof_index,state_no)

try:
  fiber_file_handle = open(fiber_file, "rb")
except:
  print("Error: Could not open fiber file \"{}\"".format(fiber_file))
  quit()

# parse fibers from a binary fiber file that was created by parallel_fiber_estimation
# parse file header to extract number of fibers
bytes_raw = fiber_file_handle.read(32)
header_str = struct.unpack('32s', bytes_raw)[0]
header_length_raw = fiber_file_handle.read(4)
header_length = struct.unpack('i', header_length_raw)[0]

parameters = []
for i in range(int(header_length/4.) - 1):
  double_raw = fiber_file_handle.read(4)
  value = struct.unpack('i', double_raw)[0]
  parameters.append(value)
  
n_fibers_total = parameters[0]
n_fibers_x = (int)(np.round(np.sqrt(n_fibers_total)))
n_fibers_y = n_fibers_x
n_points_whole_fiber = parameters[1]

if rank_no == 0:
  print("n fibers:              {} ({} x {})".format(n_fibers_total, n_fibers_x, n_fibers_y))
  print("n points per fiber:    {}".format(n_points_whole_fiber))
  
# parse whole fiber file, only if enabled
if load_data_from_file:
  fibers = []
  for fiber_no in range(n_fibers_total):
    fiber = []
    for point_no in range(n_points_whole_fiber):
      point = []
      for i in range(3):
        double_raw = fiber_file_handle.read(8)
        value = struct.unpack('d', double_raw)[0]
        point.append(value)
      fiber.append(point)
    fibers.append(fiber)
          
# load MU distribution and firing times
fiber_distribution = np.genfromtxt(fiber_distribution_file, delimiter=" ")
firing_times = np.genfromtxt(firing_times_file)

# for debugging output show when the first 20 fibers will fire
if rank_no == 0 and not disable_firing_output:
  print("Debugging output about fiber firing: Taking input from file \"{}\"".format(firing_times_file))
  import timeit
  t_start = timeit.default_timer()
  
  first_stimulation_info = []
  
  n_firing_times = np.size(firing_times,0)
  for fiber_no_index in range(n_fibers_total):
    if fiber_no_index % 100 == 0:
      t_intermediate = timeit.default_timer()
      if t_intermediate - t_start > 100:
        print("Note: break after {}/{} fibers ({:.0f}%) because it already took {:.3f}s".format(fiber_no_index,n_fibers_total,100.0*fiber_no_index/(n_fibers_total-1.),t_intermediate - t_start))
        break
    
    first_stimulation = None
    for current_time in np.linspace(0,1./stimulation_frequency*n_firing_times,n_firing_times):
      if fiber_gets_stimulated(fiber_no_index, stimulation_frequency, current_time):
        first_stimulation = current_time
        break
    mu_no = get_motor_unit_no(fiber_no_index)
    first_stimulation_info.append([fiber_no_index,mu_no,first_stimulation])
  
  first_stimulation_info.sort(key=lambda x: 1e6+1e-6*x[1]+1e-12*x[0] if x[2] is None else x[2]+1e-6*x[1]+1e-12*x[0])
  
  print("First stimulation times")
  print("    Time  MU fibers")
  n_stimulated_mus = 0
  n_not_stimulated_mus = 0
  fibers = []
  last_time = 0
  last_mu_no = first_stimulation_info[0][1]
  for stimulation_info in first_stimulation_info:
    mu_no = stimulation_info[1]
    fiber_no = stimulation_info[0]
    if mu_no == last_mu_no:
      fibers.append(fiber_no)
    else:
      if last_time is not None:
        if len(fibers) > 10:
          print("{:8.2f} {:3} {} (only showing first 10, {} total)".format(last_time,last_mu_no,str(fibers[0:10]),len(fibers)))
        else:
          print("{:8.2f} {:3} {}".format(last_time,last_mu_no,str(fibers)))
        n_stimulated_mus += 1
      else:
        if len(fibers) > 10:
          print("  never stimulated: MU {:3}, fibers {} (only showing first 10, {} total)".format(last_mu_no,str(fibers[0:10]),len(fibers)))
        else:
          print("  never stimulated: MU {:3}, fibers {}".format(last_mu_no,str(fibers)))
        n_not_stimulated_mus += 1
      fibers = [fiber_no]

    last_time = stimulation_info[2]
    last_mu_no = mu_no
    
  print("stimulated MUs: {}, not stimulated MUs: {}".format(n_stimulated_mus,n_not_stimulated_mus))

  t_end = timeit.default_timer()
  print("duration of assembling this list: {:.3f} s\n".format(t_end-t_start))  
  
# compute partitioning
if rank_no == 0:
  if n_ranks != n_subdomains_x*n_subdomains_y*n_subdomains_z:
    print("\n\nError! Number of ranks {} does not match given partitioning {} x {} x {} = {}.\n\n".format(n_ranks, n_subdomains_x, n_subdomains_y, n_subdomains_z, n_subdomains_x*n_subdomains_y*n_subdomains_z))
    quit()
  
n_fibers_per_subdomain_x = (int)(n_fibers_x / n_subdomains_x)
n_fibers_per_subdomain_y = (int)(n_fibers_y / n_subdomains_y)
n_points_per_subdomain_z = (int)(n_points_whole_fiber / n_subdomains_z)

# define helper functions for fiber numbering

# number of fibers that are handled inside the subdomain x
def n_fibers_in_subdomain_x(subdomain_coordinate_x):
  a1 = n_fibers_x - n_subdomains_x*n_fibers_per_subdomain_x              # number of subdomains with high number of fibers
  a2 = n_subdomains_x - a1                                               # number of subdomains with low number of fibers
  if subdomain_coordinate_x < a1:
    return n_fibers_per_subdomain_x+1      # high number of fibersr
  else:
    return n_fibers_per_subdomain_x    # low number of fibers
  
# number of fibers that are handled inside the subdomain y
def n_fibers_in_subdomain_y(subdomain_coordinate_y):
  a1 = n_fibers_y - n_subdomains_y*n_fibers_per_subdomain_y              # number of subdomains with high number of fibers
  a2 = n_subdomains_y - a1                                               # number of subdomains with low number of fibers
  if subdomain_coordinate_y < a1:
    return n_fibers_per_subdomain_y+1     # high number of fibers
  else:
    return n_fibers_per_subdomain_y     # low number of fibers

def n_fibers_in_previous_subdomains_y(subdomain_coordinate_y):
  # number of fibers handled in previous subdomains in y direction
  a1 = n_fibers_y - n_subdomains_y*n_fibers_per_subdomain_y              # number of subdomains with high number of fibers
  a2 = n_subdomains_y - a1                                               # number of subdomains with low number of fibers
  
  if subdomain_coordinate_y < a1:
    return subdomain_coordinate_y * (n_fibers_per_subdomain_y+1)
  else:
    return a1 * (n_fibers_per_subdomain_y+1) + (subdomain_coordinate_y-a1) * n_fibers_per_subdomain_y
    
def n_fibers_in_previous_subdomains_x(subdomain_coordinate_x):
  # number of fibers handled in previous subdomains in x direction
  a1 = n_fibers_x - n_subdomains_x*n_fibers_per_subdomain_x              # number of subdomains with high number of fibers
  a2 = n_subdomains_x - a1                                               # number of subdomains with low number of fibers
  
  if subdomain_coordinate_x < a1:
    return subdomain_coordinate_x * (n_fibers_per_subdomain_x+1)
  else:
    return a1 * (n_fibers_per_subdomain_x+1) + (subdomain_coordinate_x-a1) * n_fibers_per_subdomain_x

# global fiber no, from subdomain coordinate and coordinate inside the subdomain
def fiber_no(subdomain_coordinate_x, subdomain_coordinate_y, fiber_in_subdomain_coordinate_x, fiber_in_subdomain_coordinate_y):
  return (n_fibers_in_previous_subdomains_y(subdomain_coordinate_y) + fiber_in_subdomain_coordinate_y)*n_fibers_x \
    + n_fibers_in_previous_subdomains_x(subdomain_coordinate_x) + fiber_in_subdomain_coordinate_x

# number of points that are handled inside the subdomain z (without ghost points)
def n_points_in_subdomain_z(subdomain_coordinate_z):
  a1 = n_points_whole_fiber - n_subdomains_z*n_points_per_subdomain_z              # number of subdomains with high number of fibers
  a2 = n_subdomains_z - a1                                               # number of subdomains with low number of fibers
  if subdomain_coordinate_z < a1:
    return n_points_per_subdomain_z+1     # high number of points
  else:
    return n_points_per_subdomain_z     # low number of points
  
def n_points_in_previous_subdomains_z(subdomain_coordinate_z):
  # number of points handled in previous subdomains in z direction
  a1 = n_points_whole_fiber - n_subdomains_z*n_points_per_subdomain_z              # number of subdomains with high number of points
  a2 = n_subdomains_z - a1                                               # number of subdomains with low number of points
  
  if subdomain_coordinate_z < a1:
    return subdomain_coordinate_z * (n_points_per_subdomain_z+1)
  else:
    return a1 * (n_points_per_subdomain_z+1) + (subdomain_coordinate_z-a1) * n_points_per_subdomain_z

# number of points in 3D mesh in the subdomain, in x direction
def n_sampled_points_in_subdomain_x(subdomain_coordinate_x):
  result = (int)(np.ceil(n_fibers_in_subdomain_x(subdomain_coordinate_x) / sampling_stride_x))
  if subdomain_coordinate_x == n_subdomains_x-1 and (n_fibers_in_subdomain_x(subdomain_coordinate_x)-1) % sampling_stride_x != 0:
    result += 1
  return result

# number of points in 3D mesh in the subdomain, in y direction
def n_sampled_points_in_subdomain_y(subdomain_coordinate_y):
  result = (int)(np.ceil(n_fibers_in_subdomain_y(subdomain_coordinate_y) / sampling_stride_y))
  if subdomain_coordinate_y == n_subdomains_y-1 and (n_fibers_in_subdomain_y(subdomain_coordinate_y)-1) % sampling_stride_y != 0:
    result += 1
  return result

# number of points in 3D mesh in the subdomain, in z direction
def n_sampled_points_in_subdomain_z(subdomain_coordinate_z):
  result = (int)(np.ceil(n_points_in_subdomain_z(subdomain_coordinate_z) / sampling_stride_z))
  if subdomain_coordinate_z == n_subdomains_z-1 and (n_points_in_subdomain_z(subdomain_coordinate_z)-1) % sampling_stride_z != 0:
    result += 1
  return result

#####################
# define 3D mesh
# determine node positions of the 3D mesh
node_positions_3d_mesh = []
if not load_data_from_file:
  node_positions_3d_mesh = [fiber_file, []]

# range of points in z direction
z_point_index_start = n_points_in_previous_subdomains_z(own_subdomain_coordinate_z)
z_point_index_end = z_point_index_start + n_points_in_subdomain_z(own_subdomain_coordinate_z)

# loop over z point indices
for k in range(n_sampled_points_in_subdomain_z(own_subdomain_coordinate_z)):
  z_point_index = z_point_index_start + k*sampling_stride_z
  
  if own_subdomain_coordinate_z == n_subdomains_z-1 and k == n_sampled_points_in_subdomain_z(own_subdomain_coordinate_z)-1:
    z_point_index = z_point_index_end-1
    
  #print("{}: sampling_stride_z: {}, k: {}, z: {}/{}".format(rank_no, sampling_stride_z, k, z_point_index, z_point_index_end))
  
  # loop over fibers for own rank
  # loop over fiber in y-direction
  for j in range(n_sampled_points_in_subdomain_y(own_subdomain_coordinate_y)):
    fiber_in_subdomain_coordinate_y = j*sampling_stride_y
    
    # on border rank set last node positions to be the border nodes (it could be that they are not yet the outermost nodes because of sampling_stride)
    if own_subdomain_coordinate_y == n_subdomains_y-1 and j == n_sampled_points_in_subdomain_y(own_subdomain_coordinate_y)-1:
      fiber_in_subdomain_coordinate_y = n_fibers_in_subdomain_y(own_subdomain_coordinate_y)-1
    
    #print("{}: sampling_stride_y: {}, j: {}, y: {}/{}".format(rank_no, sampling_stride_y, j, fiber_in_subdomain_coordinate_y, n_fibers_in_subdomain_y(own_subdomain_coordinate_y)))
      
    # loop over fiber in x-direction
    for i in range(n_sampled_points_in_subdomain_x(own_subdomain_coordinate_x)):
      fiber_in_subdomain_coordinate_x = i*sampling_stride_x
      
      # on border rank set last node positions to be the border nodes (it could be that they are not yet the outermost nodes because of sampling_stride)
      if own_subdomain_coordinate_x == n_subdomains_x-1 and i == n_sampled_points_in_subdomain_x(own_subdomain_coordinate_x)-1:
        fiber_in_subdomain_coordinate_x = n_fibers_in_subdomain_x(own_subdomain_coordinate_x)-1
      
      #print("{}: sampling_stride_x: {}, i: {}, x: {}/{}".format(rank_no, sampling_stride_x, i, fiber_in_own_subdomain_coordinate_x, n_fibers_in_subdomain_x(own_subdomain_coordinate_x)))
      
      # get fiber no
      fiber_index = fiber_no(own_subdomain_coordinate_x, own_subdomain_coordinate_y, fiber_in_subdomain_coordinate_x, fiber_in_subdomain_coordinate_y)
      
      # read point from fiber file
      memory_size_fiber = n_points_whole_fiber*3*8
      offset = 32 + header_length + fiber_index*memory_size_fiber + z_point_index*3*8
      
      if load_data_from_file:
        
        fiber_file_handle.seek(offset)
      
        point = []
        for j in range(3):
          double_raw = fiber_file_handle.read(8)
          value = struct.unpack('d', double_raw)[0]
          point.append(value)
      
        reference_point = fibers[fiber_index][z_point_index]
        
        difference = np.linalg.norm(np.array(reference_point) - np.array(point))
        if difference > 1e-3:
          print("Error, point does not match: reference_point: ", reference_point, ", point: ", point)
          quit()
        node_positions_3d_mesh.append(point)
        
      else:
        node_positions_3d_mesh[1].append((offset, 1))   # command to read 1 point from offset
      #print("{}: {}".format(rank_no, point))
     
# set local number of elements for the 3D mesh
n_elements_3D_mesh = [
    n_sampled_points_in_subdomain_x(own_subdomain_coordinate_x),
    n_sampled_points_in_subdomain_y(own_subdomain_coordinate_y), 
    n_sampled_points_in_subdomain_z(own_subdomain_coordinate_z)
  ]

# border subdomains have one element less than fibers
if own_subdomain_coordinate_x == n_subdomains_x-1:
  n_elements_3D_mesh[0] -= 1
if own_subdomain_coordinate_y == n_subdomains_y-1:
  n_elements_3D_mesh[1] -= 1
if own_subdomain_coordinate_z == n_subdomains_z-1:
  n_elements_3D_mesh[2] -= 1

# set the entry for the config
meshes = {}
meshes["3Dmesh"] = {
  "nElements": n_elements_3D_mesh,
  "nRanks": [n_subdomains_x, n_subdomains_y, n_subdomains_z],
  "nodePositions": node_positions_3d_mesh,
  "inputMeshIsGlobal": False,
  "setHermiteDerivatives": False,
  "logKey": "3Dmesh"
}

####################################
# set Dirichlet BC for the flow problem

n_points_3D_mesh_global_x = sum([n_sampled_points_in_subdomain_x(subdomain_coordinate_x) for subdomain_coordinate_x in range(n_subdomains_x)])
n_points_3D_mesh_global_y = sum([n_sampled_points_in_subdomain_y(subdomain_coordinate_y) for subdomain_coordinate_y in range(n_subdomains_y)])
n_points_3D_mesh_global_z = sum([n_sampled_points_in_subdomain_z(subdomain_coordinate_z) for subdomain_coordinate_z in range(n_subdomains_z)])
n_points_3D_mesh_global = n_points_3D_mesh_global_x*n_points_3D_mesh_global_y*n_points_3D_mesh_global_z
 
# set Dirichlet BC values for bottom nodes to 0 and for top nodes to 1
potential_flow_dirichlet_bc = {}
for i in range(n_points_3D_mesh_global_x*n_points_3D_mesh_global_y):
  potential_flow_dirichlet_bc[i] = 0.0
  potential_flow_dirichlet_bc[(n_points_3D_mesh_global_z-1)*n_points_3D_mesh_global_x*n_points_3D_mesh_global_y + i] = 1.0
    
# set Dirichlet BC at top nodes for linear elasticity problem, fix muscle at top
linear_elasticity_dirichlet_bc = {}
for i in range(n_points_3D_mesh_global_x*n_points_3D_mesh_global_y):
  linear_elasticity_dirichlet_bc[(n_points_3D_mesh_global_z-1)*n_points_3D_mesh_global_x*n_points_3D_mesh_global_y + i] = 0.0
    
# Neumann BC at bottom nodes, traction downwards
nx = n_points_3D_mesh_global_x-1
ny = n_points_3D_mesh_global_y-1
nz = n_points_3D_mesh_global_z-1
#linear_elasticity_neumann_bc = [{"element": 0*nx*ny + j*nx + i, "constantVector": [0.0,0.0,-1.0], "face": "2-"} for j in range(ny) for i in range(nx)]
linear_elasticity_neumann_bc = []
    
if debug_output:
  print("{}: point sampling for elements, unsampled points: {} x {} x {}, sampling stride: {}, {}, {}".format(rank_no, n_fibers_x, n_fibers_y, n_points_whole_fiber, sampling_stride_x, sampling_stride_y, sampling_stride_z))
  print("{}: sampled points, n_points_3D_mesh_global: {} x {} x {} = sum({}) x sum({}) x sum({})".format(rank_no, n_points_3D_mesh_global_x, n_points_3D_mesh_global_y, n_points_3D_mesh_global_z, \
    [n_sampled_points_in_subdomain_x(subdomain_coordinate_x) for subdomain_coordinate_x in range(n_subdomains_x)],\
    [n_sampled_points_in_subdomain_y(subdomain_coordinate_y) for subdomain_coordinate_y in range(n_subdomains_y)],\
    [n_sampled_points_in_subdomain_z(subdomain_coordinate_z) for subdomain_coordinate_z in range(n_subdomains_z)]))

if rank_no == 0:
  print("diffusion solver type: {}".format(diffusion_solver_type))
  print("{} ranks, partitioning: x{} x y{} x z{}".format(n_ranks, n_subdomains_x, n_subdomains_y, n_subdomains_z))
  print("{} x {} = {} fibers, per partition: {} x {} = {}".format(n_fibers_x, n_fibers_y, n_fibers_total, n_fibers_per_subdomain_x, n_fibers_per_subdomain_y, n_fibers_per_subdomain_x*n_fibers_per_subdomain_y))
  print("{} points per fiber, per partition: {}".format(n_points_whole_fiber, n_points_per_subdomain_z))
  print("number of degrees of freedom:")
  print("                    1D fiber: {:10d}  (per process: {})".format(n_points_whole_fiber, n_points_per_subdomain_z))
  print("            0D-1D monodomain: {:10d}  (per process: {})".format(n_points_whole_fiber*4, n_points_per_subdomain_z*4))
  print(" all fibers 0D-1D monodomain: {:10d}  (per process: {})".format(n_fibers_total*n_points_whole_fiber*4, n_fibers_per_subdomain_x*n_fibers_per_subdomain_y*n_points_per_subdomain_z*4))
  print("                 3D bidomain: {:10d}  (per process: {})".format(n_points_3D_mesh_global, n_sampled_points_in_subdomain_x(own_subdomain_coordinate_x)*n_sampled_points_in_subdomain_y(own_subdomain_coordinate_y)*n_sampled_points_in_subdomain_z(own_subdomain_coordinate_z)))
  print("                       total: {:10d}  (per process: {})".format(n_fibers_total*n_points_whole_fiber*4+n_points_3D_mesh_global, n_fibers_per_subdomain_x*n_fibers_per_subdomain_y*n_points_per_subdomain_z*4+n_sampled_points_in_subdomain_x(own_subdomain_coordinate_x)*n_sampled_points_in_subdomain_y(own_subdomain_coordinate_y)*n_sampled_points_in_subdomain_z(own_subdomain_coordinate_z)))

###############################
# determine 1D meshes of fibers

# fiber nos of the fibers that are handled on the own subdomain
fibers_on_own_rank = [fiber_no(own_subdomain_coordinate_x, own_subdomain_coordinate_y, fiber_in_subdomain_coordinate_x, fiber_in_subdomain_coordinate_y) \
  for fiber_in_subdomain_coordinate_y in range(n_fibers_in_subdomain_y(own_subdomain_coordinate_y)) \
  for fiber_in_subdomain_coordinate_x in range(n_fibers_in_subdomain_x(own_subdomain_coordinate_x))]
  
if debug_output:
  print("{}: rank {}, n_elements_3D_mesh: {}, subdomain coordinate ({},{},{})/({},{},{})".format(rank_no, rank_no, n_elements_3D_mesh, own_subdomain_coordinate_x, own_subdomain_coordinate_y, own_subdomain_coordinate_z, n_subdomains_x, n_subdomains_y, n_subdomains_z))
  print("{}:    fibers x: [{}, {}]".format(rank_no, 0, n_fibers_in_subdomain_x(own_subdomain_coordinate_x)))
  print("{}:    fibers y: [{}, {}]".format(rank_no, 0, n_fibers_in_subdomain_y(own_subdomain_coordinate_y)))
  print("{}:       ({})".format(rank_no, fibers_on_own_rank))
  print("{}:    points z: [{}, {}] ({})".format(rank_no, z_point_index_start, z_point_index_end, n_points_in_subdomain_z(own_subdomain_coordinate_z)))
    
# determine number of nodes and elements of the local part of a fiber  
n_fiber_nodes_on_subdomain = n_points_in_subdomain_z(own_subdomain_coordinate_z)   # number of nodes without ghosts

fiber_start_node_no = n_points_in_previous_subdomains_z(own_subdomain_coordinate_z)

# loop over all fibers
for i in range(n_fibers_total):

  # if fiber is computed on own rank
  if i in fibers_on_own_rank:
    
    if debug_output:
      print("{}: fiber {} is in fibers on own rank, {}".format(rank_no, i, str(fibers_on_own_rank)))
    
    n_fiber_elements_on_subdomain = n_fiber_nodes_on_subdomain

    # top subdomain has one element less than nodes
    if own_subdomain_coordinate_z == n_subdomains_z-1:
      n_fiber_elements_on_subdomain -= 1

    # address fiber data
    memory_size_fiber = n_points_whole_fiber * 3 * 8
    offset = 32 + header_length + i*memory_size_fiber + fiber_start_node_no*3*8
    
    # if data was loaded directly here in the python script, assign the corresponding node positions
    if load_data_from_file:
      fiber_file_handle.seek(offset)
      
      fiber_node_positions = []
      for point_no in range(n_fiber_nodes_on_subdomain):
        point = []
        for j in range(3):
          double_raw = fiber_file_handle.read(8)
          value = struct.unpack('d', double_raw)[0]
          point.append(value)
        
        fiber_node_positions.append(point)
        
        if debug_output:
          print("------")  
          print("i: ",i)
          print("fiber_node_positions: ",fiber_node_positions[0:10])
          print("fiber_start_node_no: ",fiber_start_node_no,", n_fiber_elements_on_subdomain: ",n_fiber_elements_on_subdomain)
          print("fibers[i]: ",fibers[i][fiber_start_node_no:fiber_start_node_no+10])
          
        if np.linalg.norm(np.array(fibers[i][fiber_start_node_no:fiber_start_node_no+n_fiber_elements_on_subdomain]) - np.array(fiber_node_positions[0:n_fiber_elements_on_subdomain])) > 1e-3:
          print("mismatch fiber node positions!")
          quit()
          
    else:   # add command at which position the node data in the binary file can be found, the core loads the data
      fiber_node_positions = [fiber_file, [(offset, n_fiber_nodes_on_subdomain)]]
    
    if debug_output:
      print("{}: define mesh \"{}\", n_fiber_elements_on_subdomain: {}, fiber_node_positions: {}".format(rank_no, "MeshFiber_{}".format(i), \
        str(n_fiber_elements_on_subdomain), str(fiber_node_positions)))
      
    # define mesh
    meshes["MeshFiber_{}".format(i)] = {
      "nElements": n_fiber_elements_on_subdomain,
      "nodePositions": fiber_node_positions,
      "inputMeshIsGlobal": False,
      "nRanks": [n_subdomains_z],
      "setHermiteDerivatives": False
    }
    
  else:
    # for fibers that are not computed on own rank, set empty lists for node positions and number of elements
    fiber_node_positions = []
    n_fiber_elements_on_subdomain = []

  # only add log key for fiber 0 to prevent too much data in the log files
  if i == 0 and "MeshFiber_{}".format(i) in meshes:
    meshes["MeshFiber_{}".format(i)]["logKey"] = "Fiber{}".format(i)
  
if rank_no == 0 and n_ranks < 10 and False:
  print("rank configuration: ")
  
  for subdomain_coordinate_y in range(n_subdomains_y):
    for subdomain_coordinate_x in range(n_subdomains_x):
      print("subdomain (x,y)=({},{})".format(subdomain_coordinate_x, subdomain_coordinate_y))
      print("n_subdomains_z: {}".format(n_subdomains_z))
      for rankNo in range(subdomain_coordinate_y*n_subdomains_x + subdomain_coordinate_x, n_ranks, n_subdomains_x*n_subdomains_y):
        print("  rank {}".format(rankNo))
      for fiber_in_subdomain_coordinate_y in range(n_fibers_in_subdomain_y(subdomain_coordinate_y)):
        for fiber_in_subdomain_coordinate_x in range(n_fibers_in_subdomain_x(subdomain_coordinate_x)):
          print("  fiber {} ({},{}) in subdomain uses ranks {}".format(\
            fiber_no(subdomain_coordinate_x, subdomain_coordinate_y, fiber_in_subdomain_coordinate_x, fiber_in_subdomain_coordinate_y), \
            fiber_in_subdomain_coordinate_x, fiber_in_subdomain_coordinate_y, \
            list(range(subdomain_coordinate_y*n_subdomains_x + subdomain_coordinate_x, n_ranks, n_subdomains_x*n_subdomains_y))))

config = {
  "scenarioName": scenario_name,
  "Meshes": meshes,
  "MappingsBetweenMeshes": {"MeshFiber_{}".format(i) : "3Dmesh" for i in range(n_fibers_total)},
  "Solvers": {
    "implicitSolver": {
      "maxIterations": 1e4,
      "relativeTolerance": 1e-10,
      "solverType": diffusion_solver_type,
      "preconditionerType": diffusion_preconditioner_type,
    },
    "potentialFlowSolver": {
      "relativeTolerance": 1e-10,
      "maxIterations": 1e4,
      "solverType": potential_flow_solver_type,
      "preconditionerType": potential_flow_preconditioner_type,
    },
    "activationSolver": {
      "relativeTolerance": 1e-100,
      "maxIterations": 1e4,
      "solverType": emg_solver_type,
      "preconditionerType": emg_preconditioner_type,
    },
    "linearElasticitySolver": {
      "relativeTolerance": 1e-5,
      "maxIterations": 1e4,
      "solverType": emg_solver_type,
      "preconditionerType": emg_preconditioner_type,
    }, 
  },
  "Coupling": {
    "timeStepWidth": dt_3D,  # 1e-1
    "logTimeStepWidthAsKey": "dt_3D",
    "durationLogKey": "duration_total",
    "timeStepOutputInterval" : 10,
    "endTime": end_time,
    "transferSlotName": "intermediates" if linear_elasticity else "states",
    "Term1": {      # monodomain fibers
      "MultipleInstances": {
        "logKey": "duration_subdomains_xy",
        "ranksAllComputedInstances": list(range(n_ranks)),
        "nInstances": n_subdomains_xy,
        "instances": 
        [{
          "ranks": list(range(subdomain_coordinate_y*n_subdomains_x + subdomain_coordinate_x, n_ranks, n_subdomains_x*n_subdomains_y)),
          "StrangSplitting": {
            #"numberTimeSteps": 1,
            "timeStepWidth": dt_splitting,  # 1e-1
            "logTimeStepWidthAsKey": "dt_splitting",
            "durationLogKey": "duration_monodomain",
            "timeStepOutputInterval" : 100,
            "endTime": dt_splitting,
            "transferSlotName": "states",   # which output slot of the cellml adapter ("states" or "intermediates") to use for transfer to diffusion, in this case we need "states", states[0] which is Vm

            "Term1": {      # CellML
              "MultipleInstances": {
                "logKey": "duration_subdomains_z",
                "nInstances": n_fibers_in_subdomain_x(subdomain_coordinate_x)*n_fibers_in_subdomain_y(subdomain_coordinate_y),
                "instances": 
                [{
                  "ranks": list(range(n_subdomains_z)),    # these rank nos are local nos to the outer instance of MultipleInstances, i.e. from 0 to number of ranks in z direction
                  "Heun" : {
                    "timeStepWidth": dt_0D,  # 5e-5
                    "logTimeStepWidthAsKey": "dt_0D",
                    "durationLogKey": "duration_0D",
                    "initialValues": [],
                    "timeStepOutputInterval": 1e4,
                    "inputMeshIsGlobal": True,
                    "dirichletBoundaryConditions": {},
                      
                    "CellML" : {
                      "sourceFilename": cellml_file,             # input C++ source file, can be either generated by OpenCMISS or OpenCOR from cellml model
                      "compilerFlags": "-fPIC -O3 -shared ",
                      #"simdSourceFilename" : "simdcode.cpp",     # transformed C++ source file that gets generated from sourceFilename and is ready for multiple instances
                      #"libraryFilename": "cellml_simd_lib.so",   # compiled library
                      "useGivenLibrary": False,
                      #"statesInitialValues": [],
                      #"setSpecificParametersFunction": set_specific_parameters,    # callback function that sets parameters like stimulation current
                      #"setSpecificParametersCallInterval": int(1./stimulation_frequency/dt_0D),     # set_parameters should be called every 0.1, 5e-5 * 1e3 = 5e-2 = 0.05
                      "setSpecificStatesFunction": set_specific_states,    # callback function that sets states like Vm, activation can be implemented by using this method and directly setting Vm values, or by using setParameters/setSpecificParameters
                      #"setSpecificStatesCallInterval": 2*int(1./stimulation_frequency/dt_0D),     # set_specific_states should be called stimulation_frequency times per ms, the factor 2 is needed because every Heun step includes two calls to rhs
                      "setSpecificStatesCallInterval": 0,                  # 0 means disabled
                      "setSpecificStatesCallFrequency": stimulation_frequency, # set_specific_states should be called stimulation_frequency times per ms, the factor 2 is needed because every Heun step includes two calls to rhs
                      "setSpecificStatesRepeatAfterFirstCall": 0.01,  # simulation time span for which the setSpecificStates callback will be called after a call was triggered
                      "additionalArgument": fiber_no(subdomain_coordinate_x, subdomain_coordinate_y, fiber_in_subdomain_coordinate_x, fiber_in_subdomain_coordinate_y),
                      
                      "outputIntermediateIndex": 0,  # which intermediate value to use in further computation
                      "outputStateIndex": 0,     # Shorten / Hodgkin Huxley: state 0 = Vm, Shorten: rate 28 = gamma, intermediate 0 = gamma
                      "parametersUsedAsIntermediate": parameters_used_as_intermediate,  #[32],       # list of intermediate value indices, that will be set by parameters. Explicitely defined parameters that will be copied to intermediates, this vector contains the indices of the algebraic array. This is ignored if the input is generated from OpenCMISS generated c code.
                      "parametersUsedAsConstant": parameters_used_as_constant,          #[65],           # list of constant value indices, that will be set by parameters. This is ignored if the input is generated from OpenCMISS generated c code.
                      "parametersInitialValues": parameters_initial_values,            #[0.0, 1.0],      # initial values for the parameters: I_Stim, l_hs
                      "meshName": "MeshFiber_{}".format(fiber_no(subdomain_coordinate_x, subdomain_coordinate_y, fiber_in_subdomain_coordinate_x, fiber_in_subdomain_coordinate_y)),
                      "prefactor": 1.0,
                    },
                  },
                } for fiber_in_subdomain_coordinate_y in range(n_fibers_in_subdomain_y(subdomain_coordinate_y)) \
                    for fiber_in_subdomain_coordinate_x in range(n_fibers_in_subdomain_x(subdomain_coordinate_x))],
              }
            },
            "Term2": {     # Diffusion
              "MultipleInstances": {
                "nInstances": n_fibers_in_subdomain_x(subdomain_coordinate_x)*n_fibers_in_subdomain_y(subdomain_coordinate_y),
                "instances": 
                [{
                  "ranks": list(range(n_subdomains_z)),    # these rank nos are local nos to the outer instance of MultipleInstances, i.e. from 0 to number of ranks in z direction
                  "ImplicitEuler" : {
                    "initialValues": [],
                    #"numberTimeSteps": 1,
                    "timeStepWidth": dt_1D,  # 1e-5
                    "logTimeStepWidthAsKey": "dt_1D",
                    "durationLogKey": "duration_1D",
                    "timeStepOutputInterval": 1e4,
                    "dirichletBoundaryConditions": {0: -75.0036, -1: -75.0036},
                    "inputMeshIsGlobal": True,
                    "solverName": "implicitSolver",
                    "FiniteElementMethod" : {
                      "maxIterations": 1e4,
                      "relativeTolerance": 1e-10,
                      "inputMeshIsGlobal": True,
                      "meshName": "MeshFiber_{}".format(fiber_no(subdomain_coordinate_x, subdomain_coordinate_y, fiber_in_subdomain_coordinate_x, fiber_in_subdomain_coordinate_y)),
                      "prefactor": Conductivity/(Am*Cm),
                      "solverName": "implicitSolver",
                    },
                    "OutputWriter" : [
                      #{"format": "Paraview", "outputInterval": int(1./dt_1D*output_timestep), "filename": "out/fiber_"+str(fiber_no(subdomain_coordinate_x, subdomain_coordinate_y, fiber_in_subdomain_coordinate_x, fiber_in_subdomain_coordinate_y)), "binary": True, "fixedFormat": False, "combineFiles": True},
                      #{"format": "Paraview", "outputInterval": 1./dt_1D*output_timestep, "filename": "out/fiber_"+str(i)+"_txt", "binary": False, "fixedFormat": False},
                      #{"format": "ExFile", "filename": "out/fiber_"+str(i), "outputInterval": 1./dt_1D*output_timestep, "sphereSize": "0.02*0.02*0.02"},
                      #{"format": "PythonFile", "filename": "out/fiber_"+str(i), "outputInterval": 1./dt_1D*output_timestep, "binary":True, "onlyNodalValues":True},
                    ]
                  },
                } for fiber_in_subdomain_coordinate_y in range(n_fibers_in_subdomain_y(subdomain_coordinate_y)) \
                    for fiber_in_subdomain_coordinate_x in range(n_fibers_in_subdomain_x(subdomain_coordinate_x))],
                "OutputWriter" : output_writer_fibers,
              },
            },
          }
        } if (subdomain_coordinate_x,subdomain_coordinate_y) == (own_subdomain_coordinate_x,own_subdomain_coordinate_y) else None
        for subdomain_coordinate_y in range(n_subdomains_y)
            for subdomain_coordinate_x in range(n_subdomains_x)]
      }
    },
    "Term2": {        # Bidomain
      "StaticBidomainSolver": {       # version for fibers_emg
        "timeStepWidth": dt_3D,
        "timeStepOutputInterval": 50,
        "durationLogKey": "duration_bidomain",
        "solverName": "activationSolver",
        "initialGuessNonzero": emg_initial_guess_nonzero,
        "PotentialFlow": {
          "FiniteElementMethod" : {  
            "meshName": "3Dmesh",
            "solverName": "potentialFlowSolver",
            "prefactor": 1.0,
            "dirichletBoundaryConditions": potential_flow_dirichlet_bc,
            "inputMeshIsGlobal": True,
          },
        },
        "Activation": {
          "FiniteElementMethod" : {  
            "meshName": "3Dmesh",
            "solverName": "activationSolver",
            "prefactor": 1.0,
            "inputMeshIsGlobal": True,
            "dirichletBoundaryConditions": {},
            "diffusionTensor": [      # sigma_i           # fiber direction is (1,0,0)
              8.93, 0, 0,
              0, 0.893, 0,
              0, 0, 0.893
            ], 
            "extracellularDiffusionTensor": [      # sigma_e
              6.7, 0, 0,
              0, 6.7, 0,
              0, 0, 6.7,
            ],
          },
        },
        "OutputWriter" : output_writer_emg,
      },
      "OutputSurface": {        # version for fibers_emg_2d_output
        "OutputWriter": [
          {"format": "Paraview", "outputInterval": int(1./dt_3D*output_timestep), "filename": "out/" + scenario_name + "/surface_emg", "binary": True, "fixedFormat": False, "combineFiles": True},
        ],
        "face": "1-",
        "StaticBidomainSolver": {
          "timeStepWidth": dt_3D,
          "timeStepOutputInterval": 50,
          "durationLogKey": "duration_bidomain",
          "solverName": "activationSolver",
          "initialGuessNonzero": emg_initial_guess_nonzero,
          "PotentialFlow": {
            "FiniteElementMethod" : {  
              "meshName": "3Dmesh",
              "solverName": "potentialFlowSolver",
              "prefactor": 1.0,
              "dirichletBoundaryConditions": potential_flow_dirichlet_bc,
              "inputMeshIsGlobal": True,
            },
          },
          "Activation": {
            "FiniteElementMethod" : {  
              "meshName": "3Dmesh",
              "solverName": "activationSolver",
              "prefactor": 1.0,
              "inputMeshIsGlobal": True,
              "dirichletBoundaryConditions": {},
              "diffusionTensor": [      # sigma_i           # fiber direction is (1,0,0)
                8.93, 0, 0,
                0, 0.893, 0,
                0, 0, 0.893
              ], 
              "extracellularDiffusionTensor": [      # sigma_e
                6.7, 0, 0,
                0, 6.7, 0,
                0, 0, 6.7,
              ],
            },
          },
          "OutputWriter" : output_writer_emg,
        }
      },
      "QuasiStaticLinearElasticitySolver": {
        "PotentialFlow": {        # potential flow for fiber directions in the 3D mesh
          "FiniteElementMethod" : {  
            "meshName": "3Dmesh",
            "solverName": "potentialFlowSolver",
            "prefactor": 1.0,
            "dirichletBoundaryConditions": potential_flow_dirichlet_bc,
            "neumannBoundaryConditions": [],
            "inputMeshIsGlobal": True,
          },
        },
        "FiniteElementMethod" : {   # linear elasticity finite element method
          "meshName": "3Dmesh",
          "solverName": "linearElasticitySolver",
          "prefactor": 1.0,
          "inputMeshIsGlobal": True,
          "dirichletBoundaryConditions": linear_elasticity_dirichlet_bc,
          "neumannBoundaryConditions": linear_elasticity_neumann_bc,
          "bulkModulus": 40e3,  # https://www.researchgate.net/publication/230248067_Bulk_Modulus
          "shearModulus": 39e3, # https://onlinelibrary.wiley.com/doi/full/10.1002/mus.24104
        },
        "maximumActiveStress": 1.0,
        "strainScalingCurveWidth": 1.0,
        "OutputWriter" : [
          {"format": "Paraview", "outputInterval": int(1./dt_3D*output_timestep), "filename": "out/deformation", "binary": True, "fixedFormat": False, "onlyNodalValues":True, "combineFiles":True},
          #{"format": "PythonFile", "filename": "out/deformation", "outputInterval": 1, "binary":False, "onlyNodalValues":True},
        ]
      }
    }
  }
}

# sanity checking
if False:
  # check coupling instances
  multiple_instances = config["Coupling"]["Term1"]["MultipleInstances"]
  n_instances = multiple_instances["nInstances"]
  instances_size = len(multiple_instances["instances"])
  print("n_instances: {}".format(n_instances))

  print("n subdomains: {} x {}".format(n_subdomains_x, n_subdomains_y))
  print("n_fibers_per_subdomain_x: {} {}".format(n_fibers_per_subdomain_x, n_fibers_per_subdomain_y))

  for subdomain_coordinate_y in range(n_subdomains_y):
    print("n_fibers_in_subdomain_y({}) = {}".format(subdomain_coordinate_y, n_fibers_in_subdomain_y(subdomain_coordinate_y)))

  print("--")
  for subdomain_coordinate_x in range(n_subdomains_x):
    print("n_fibers_in_subdomain_x({}) = {}".format(subdomain_coordinate_x, n_fibers_in_subdomain_x(subdomain_coordinate_x)))
  print("--")

  # check fiber no
  counter = 0
  for subdomain_coordinate_y in range(n_subdomains_y):
    for fiber_in_subdomain_coordinate_y in range(n_fibers_in_subdomain_y(subdomain_coordinate_y)):
      for subdomain_coordinate_x in range(n_subdomains_x):
        for fiber_in_subdomain_coordinate_x in range(n_fibers_in_subdomain_x(subdomain_coordinate_x)):
          no = fiber_no(subdomain_coordinate_x, subdomain_coordinate_y, fiber_in_subdomain_coordinate_x, fiber_in_subdomain_coordinate_y)
          if no != counter:
            print("error: fiber_no({},{},{},{}) = {}, counter = {}".format(subdomain_coordinate_x, subdomain_coordinate_y, fiber_in_subdomain_coordinate_x, fiber_in_subdomain_coordinate_y,no,counter))
          else:
            print("   ok: fiber_no({},{},{},{}) = {}, counter = {}".format(subdomain_coordinate_x, subdomain_coordinate_y, fiber_in_subdomain_coordinate_x, fiber_in_subdomain_coordinate_y,no,counter))
          counter += 1
          

  if n_instances != instances_size or instances_size == 0:
    print("Error with top-level multiple instances: nInstances: {}, size of instances: {}".format(n_instances, instances_size))

  # loop over inner instances
  for i in range(n_instances):
    multiple_instances0 = multiple_instances["instances"][i]["StrangSplitting"]["Term1"]["MultipleInstances"]
    n_instances0 = multiple_instances0["nInstances"]
    instances_size0 = len(multiple_instances0["instances"])
    
    if n_instances0 != instances_size0 or instances_size0 == 0:
      print("Error with Term1 {} multiple instances: nInstances: {}, size of instances: {}".format(i ,n_instances0, instances_size0))
    
    multiple_instances1 = multiple_instances["instances"][i]["StrangSplitting"]["Term2"]["MultipleInstances"]
    n_instances1 = multiple_instances1["nInstances"]
    instances_size1 = len(multiple_instances1["instances"])

    if n_instances1 != instances_size1 or instances_size1 == 0:
      print("Error with Term2 {} multiple instances: nInstances: {}, size of instances: {}".format(i, n_instances1, instances_size1))

if rank_no == 0:
  t_stop_script = timeit.default_timer()
  print("Python config parsed in {:.1f}s.".format(t_stop_script - t_start_script))
