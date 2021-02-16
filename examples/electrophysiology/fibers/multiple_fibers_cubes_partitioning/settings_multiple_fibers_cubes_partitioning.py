# Multiple 1D fibers (monodomain), biceps geometry
# This is similar to the fibers_emg example, but without EMG.
# To see all available arguments, execute: ./multiple_fibers settings_multiple_fibers_cubes_partitioning.py -help
#
# if fiber_file=cuboid.bin, it uses a small cuboid test example (Contrary to the "cuboid" example, this produces a real cuboid).
#
# You have to set n_subdomains such that it matches the number of processes, e.g. 2x2x1 = 4 processes.
# Decomposition is in x,y,z direction, the fibers are aligned with the z axis.
# E.g. --n_subdomains 2 2 1 which is 2x2x1 means no subdivision per fiber, 
# --n_subdomains 8 8 4 means every fiber will be subdivided to 4 processes and all fibers will be computed by 8x8 processes.
#
# Example with 4 processes and end time 5, and otherwise default parameters:
#   mpirun -n 4 ./multiple_fibers ../settings_multiple_fibers_cubes_partitioning.py --n_subdomains 2 2 1 --end_time=5.0
#
# Three files contribute to the settings:
# A lot of variables are set by the helper.py script, the variables and their values are defined in variables.py and this file
# creates the composite config that is needed by opendihu.
# You can provide parameter values in a custom_variables.py file in the variables subfolder. (Instead of custom_variables.py you can choose any filename.)
# This custom variables file should be the next argument on the command line after settings_fibers_emg.py, e.g.:
#
#  ./multiple_fibers ../settings_multiple_fibers_cubes_partitioning.py custom_variables.py --n_subdomains 1 1 1 --end_time=5.0
#
# E.g. try
# ./multiple_fibers_57_states ../settings_multiple_fibers_cubes_partitioning.py compare_to_opencmiss.py
# mpirun -n 2 ./multiple_fibers_57_states ../settings_multiple_fibers_cubes_partitioning.py --n_subdomains 2 1 1

import sys, os
import timeit
import argparse
import importlib
import copy

# parse rank arguments
rank_no = (int)(sys.argv[-2])
n_ranks = (int)(sys.argv[-1])

# add variables subfolder to python path where the variables script is located
script_path = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, script_path)
sys.path.insert(0, os.path.join(script_path,'variables'))

import variables              # file variables.py, defined default values for all parameters, you can set the parameters there  
from create_partitioned_meshes_for_settings import *   # file create_partitioned_meshes_for_settings with helper functions about own subdomain

# if first argument contains "*.py", it is a custom variable definition file, load these values
if ".py" in sys.argv[0]:
  variables_path_and_filename = sys.argv[0]
  variables_path,variables_filename = os.path.split(variables_path_and_filename)  # get path and filename 
  sys.path.insert(0, os.path.join(script_path,variables_path))                    # add the directory of the variables file to python path
  variables_module,_ = os.path.splitext(variables_filename)                       # remove the ".py" extension to get the name of the module
  
  if rank_no == 0:
    print("Loading variables from \"{}\".".format(variables_path_and_filename))
    
  custom_variables = importlib.import_module(variables_module, package=variables_filename)    # import variables module
  variables.__dict__.update(custom_variables.__dict__)
  sys.argv = sys.argv[1:]     # remove first argument, which now has already been parsed

# define command line arguments
parser = argparse.ArgumentParser(description='fibers_emg')
parser.add_argument('--scenario_name',                       help='The name to identify this run in the log.',   default=variables.scenario_name)
parser.add_argument('--n_subdomains', nargs=3,               help='Number of subdomains in x,y,z direction.',    type=int)
parser.add_argument('--n_subdomains_x', '-x',                help='Number of subdomains in x direction.',        type=int, default=variables.n_subdomains_x)
parser.add_argument('--n_subdomains_y', '-y',                help='Number of subdomains in y direction.',        type=int, default=variables.n_subdomains_y)
parser.add_argument('--n_subdomains_z', '-z',                help='Number of subdomains in z direction.',        type=int, default=variables.n_subdomains_z)
parser.add_argument('--diffusion_solver_type',               help='The solver for the diffusion.',               default=variables.diffusion_solver_type, choices=["gmres","cg","lu","gamg","richardson","chebyshev","cholesky","jacobi","sor","preonly"])
parser.add_argument('--diffusion_preconditioner_type',       help='The preconditioner for the diffusion.',       default=variables.diffusion_preconditioner_type, choices=["jacobi","sor","lu","ilu","gamg","none"])
parser.add_argument('--paraview_output',                     help='Enable the paraview output writer.',          default=variables.paraview_output, action='store_true')
parser.add_argument('--adios_output',                        help='Enable the MegaMol/ADIOS output writer.',          default=variables.adios_output, action='store_true')
parser.add_argument('--fiber_file',                          help='The filename of the file that contains the fiber data.', default=variables.fiber_file)
parser.add_argument('--cellml_file',                         help='The filename of the file that contains the cellml model.', default=variables.cellml_file)
parser.add_argument('--fiber_distribution_file',             help='The filename of the file that contains the MU firing times.', default=variables.fiber_distribution_file)
parser.add_argument('--firing_times_file',                   help='The filename of the file that contains the cellml model.', default=variables.firing_times_file)
parser.add_argument('--end_time', '--tend', '-t',            help='The end simulation time.',                    type=float, default=variables.end_time)
parser.add_argument('--output_timestep',                     help='The timestep for writing outputs.',           type=float, default=variables.output_timestep)
parser.add_argument('--dt_0D',                               help='The timestep for the 0D model.',              type=float, default=variables.dt_0D)
parser.add_argument('--dt_1D',                               help='The timestep for the 1D model.',              type=float, default=variables.dt_1D)
parser.add_argument('--dt_splitting',                        help='The timestep for the splitting.',             type=float, default=variables.dt_splitting)
parser.add_argument('--disable_firing_output',               help='Disables the initial list of fiber firings.', default=variables.disable_firing_output, action='store_true')
parser.add_argument('--v',                                   help='Enable full verbosity in c++ code')
parser.add_argument('-v',                                    help='Enable verbosity level in c++ code', action="store_true")
parser.add_argument('-vmodule',                              help='Enable verbosity level for given file in c++ code')
parser.add_argument('-pause',                                help='Stop at parallel debugging barrier', action="store_true")
parser.add_argument('--n_fibers_y',                          help='Number of fibers when simulating a cuboid example.',        type=int, default=variables.n_fibers_y)

# parse command line arguments and assign values to variables module
args, other_args = parser.parse_known_args(args=sys.argv[:-2], namespace=variables)
if len(other_args) != 0 and rank_no == 0:
    print("Warning: These arguments were not parsed by the settings python file\n  " + "\n  ".join(other_args), file=sys.stderr)

# initialize some dependend variables
if variables.n_subdomains is not None:
  variables.n_subdomains_x = variables.n_subdomains[0]
  variables.n_subdomains_y = variables.n_subdomains[1]
  variables.n_subdomains_z = variables.n_subdomains[2]
  
variables.n_subdomains = variables.n_subdomains_x*variables.n_subdomains_y*variables.n_subdomains_z

# automatically initialize partitioning if it has not been set
if n_ranks != variables.n_subdomains:
  
  # create all possible partitionings to the given number of ranks
  optimal_value = n_ranks**(1/3)
  possible_partitionings = []
  for i in range(1,n_ranks+1):
    for j in range(1,n_ranks+1):
      if i*j <= n_ranks and n_ranks % (i*j) == 0:
        k = (int)(n_ranks / (i*j))
        performance = (k-optimal_value)**2 + (j-optimal_value)**2 + 1.1*(i-optimal_value)**2
        possible_partitionings.append([i,j,k,performance])
        
  # if no possible partitioning was found
  if len(possible_partitionings) == 0:
    if rank_no == 0:
      print("\n\n\033[0;31mError! Number of ranks {} does not match given partitioning {} x {} x {} = {} and no automatic partitioning could be done.\n\n\033[0m".format(n_ranks, variables.n_subdomains_x, variables.n_subdomains_y, variables.n_subdomains_z, variables.n_subdomains_x*variables.n_subdomains_y*variables.n_subdomains_z))
    quit()
    
  # select the partitioning with the lowest value of performance which is the best
  lowest_performance = possible_partitionings[0][3]+1
  for i in range(len(possible_partitionings)):
    if possible_partitionings[i][3] < lowest_performance:
      lowest_performance = possible_partitionings[i][3]
      variables.n_subdomains_x = possible_partitionings[i][0]
      variables.n_subdomains_y = possible_partitionings[i][1]
      variables.n_subdomains_z = possible_partitionings[i][2]

# output information of run
if rank_no == 0:
  print("scenario_name: {},  n_subdomains: {} {} {},  n_ranks: {},  end_time: {}".format(variables.scenario_name, variables.n_subdomains_x, variables.n_subdomains_y, variables.n_subdomains_z, n_ranks, variables.end_time))
  print("dt_0D:           {:0.0e}, diffusion_solver_type:      {}".format(variables.dt_0D, variables.diffusion_solver_type))
  print("dt_1D:           {:0.0e},".format(variables.dt_1D))
  print("dt_splitting:    {:0.0e}, paraview_output: {}".format(variables.dt_splitting, variables.paraview_output))
  print("output_timestep: {:0.0e}  stimulation_frequency: {} 1/ms = {} Hz".format(variables.output_timestep, variables.stimulation_frequency, variables.stimulation_frequency*1e3))
  print("fiber_file:              {}".format(variables.fiber_file))
  print("cellml_file:             {}".format(variables.cellml_file))
  print("fiber_distribution_file: {}".format(variables.fiber_distribution_file))
  print("firing_times_file:       {}".format(variables.firing_times_file))
  print("********************************************************************************")
  
  print("prefactor: sigma_eff/(Am*Cm) = {} = {} / ({}*{})".format(variables.Conductivity/(variables.Am*variables.Cm), variables.Conductivity, variables.Am, variables.Cm))
  
  # start timer to measure duration of parsing of this script  
  t_start_script = timeit.default_timer()
    
# initialize all helper variables
from helper import *

variables.n_subdomains_xy = variables.n_subdomains_x * variables.n_subdomains_y
variables.n_fibers_total = variables.n_fibers_x * variables.n_fibers_y

# define the config dict
config = {
  "scenarioName":                   variables.scenario_name,    # scenario name which will appear in the log file
  "logFormat":                      "csv",                      # "csv" or "json", format of the lines in the log file, csv gives smaller files
  "solverStructureDiagramFile":     "solver_structure.txt",     # output file of a diagram that shows data connection between solvers
  "mappingsBetweenMeshesLogFile":   "out/mappings_between_meshes.txt",     # output file that contains a log about creation of mappings between meshes
  "meta": {                 # additional fields that will appear in the log
    "partitioning": [variables.n_subdomains_x, variables.n_subdomains_y, variables.n_subdomains_z]
  },
  "Meshes": variables.meshes,
  "Solvers": {
    "implicitSolver": {     # solver for the implicit timestepping scheme of the diffusion time step
      "maxIterations":      1e4,
      "relativeTolerance":  1e-10,
      "absoluteTolerance":  1e-10,         # 1e-10 absolute tolerance of the residual          
      "solverType":         variables.diffusion_solver_type,
      "preconditionerType": variables.diffusion_preconditioner_type,
      "dumpFilename":       "",   # "out/dump_"
      "dumpFormat":         "matlab",
    }, 
  },
  "MultipleInstances": {
    "logKey":                     "duration_subdomains_xy",
    "ranksAllComputedInstances":  list(range(n_ranks)),
    "nInstances":                 variables.n_subdomains_xy,
    "instances": 
    [{
      "ranks": list(range(subdomain_coordinate_y*variables.n_subdomains_x + subdomain_coordinate_x, n_ranks, variables.n_subdomains_x*variables.n_subdomains_y)),
      "StrangSplitting": {
        #"numberTimeSteps": 1,
        "timeStepWidth":          variables.dt_splitting,  # 1e-1
        "logTimeStepWidthAsKey":  "dt_splitting",
        "durationLogKey":         "duration_monodomain",
        "timeStepOutputInterval": 100,
        "endTime":                variables.end_time,
        "connectedSlotsTerm1To2": [0,1],   # transfer slot 0 = state Vm from Term1 (CellML) to Term2 (Diffusion)
        "connectedSlotsTerm2To1": [0,1],   # transfer the same back

        "Term1": {      # CellML, i.e. reaction term of Monodomain equation
          "MultipleInstances": {
            "logKey":             "duration_subdomains_z",
            "nInstances":         n_fibers_in_subdomain_x(subdomain_coordinate_x)*n_fibers_in_subdomain_y(subdomain_coordinate_y),
            "instances": 
            [{
              "ranks":                          list(range(variables.n_subdomains_z)),    # these rank nos are local nos to the outer instance of MultipleInstances, i.e. from 0 to number of ranks in z direction
              "Heun" : {
                "timeStepWidth":                variables.dt_0D,  # 5e-5
                "logTimeStepWidthAsKey":        "dt_0D",
                "durationLogKey":               "duration_0D",
                "initialValues":                [],
                "timeStepOutputInterval":       1e4,
                "inputMeshIsGlobal":            True,
                "checkForNanInf":               False,
                "dirichletBoundaryConditions":  {},
                "dirichletOutputFilename":      None,                                 # filename for a vtp file that contains the Dirichlet boundary condition nodes and their values, set to None to disable
                "nAdditionalFieldVariables":    0,
                "additionalSlotNames":          [],
                  
                "CellML" : {
                  "modelFilename":                          variables.cellml_file,                                    # input C++ source file or cellml XML file
                  "statesInitialValues":                   [],                                             # if given, the initial values for the the states of one instance
                  #"statesInitialValues": [-81.5938, -81.5407, 7.166, 150.916, 6.03857, 12.6618, 131.577, 132.928, 0.00751472, 0.996183, 0.0292367, 0.569413, 0.731601, 0.0075721, 0.996084, 0.0294341, 0.567101, 0.730931, 1.75811e-06, 5.75735e-06, 7.07019e-06, 3.85884e-06, 7.89791e-07, 0.879053, 0.115147, 0.00565615, 0.000123483, 1.01094e-06, -916.582, 0.0284792, 56.5564, 0.0284779, 1687.31, 2.98725, 615, 615, 811, 811, 1342.65, 17807.7, 0.107772, 0.10777, 7243.03, 7243.03, 756.867, 756.867, 956.975, 956.975, 0.0343398, 0.0102587, 0.0136058, 0.0314258, 0.0031226, 0.00249808, 0.223377, 0.264145, 1.74046e-06],
                  "initializeStatesToEquilibrium":          False,                                           # if the equilibrium values of the states should be computed before the simulation starts
                  "initializeStatesToEquilibriumTimestepWidth": 1e-4,                                       # if initializeStatesToEquilibrium is enable, the timestep width to use to solve the equilibrium equation
                  
                  # optimization parameters
                  "optimizationType":                       "vc",                                           # "vc", "simd", "openmp" type of generated optimizated source file
                  "approximateExponentialFunction":         True,                                           # if optimizationType is "vc", whether the exponential function exp(x) should be approximate by (1+x/n)^n with n=1024
                  "compilerFlags":                          "-fPIC -O3 -march=native -shared ",             # compiler flags used to compile the optimized model code
                  "maximumNumberOfThreads":                 0,                                              # if optimizationType is "openmp", the maximum number of threads to use. Default value 0 means no restriction.
                  
                  # stimulation callbacks
                  #"libraryFilename":                       "cellml_simd_lib.so",                           # compiled library
                  #"setParametersFunction":                 set_parameters,                                 # callback function that sets parameters like stimulation current
                  #"setParametersCallInterval":             int(1./stimulation_frequency/dt_0D),            # set_parameters should be called every 0.1, 5e-5 * 1e3 = 5e-2 = 0.05
                  #"setSpecificParametersFunction":         set_specific_parameters,                        # callback function that sets parameters like stimulation current
                  #"setSpecificParametersCallInterval":     int(1./stimulation_frequency/dt_0D),            # set_parameters should be called every 0.1, 5e-5 * 1e3 = 5e-2 = 0.05
                  "setSpecificStatesFunction":              set_specific_states,                            # callback function that sets states like Vm, activation can be implemented by using this method and directly setting Vm values, or by using setParameters/setSpecificParameters
                  "setSpecificStatesCallInterval":          0,                                              # int(1./stimulation_frequency/dt_0D),     # set_specific_states should be called every 0.1, 5e-5 * 1e3 = 5e-2 = 0.05
                  "setSpecificStatesCallFrequency":         variables.get_specific_states_call_frequency(fiber_no, motor_unit_no),   # set_specific_states should be called variables.stimulation_frequency times per ms
                  "setSpecificStatesFrequencyJitter":       variables.get_specific_states_frequency_jitter(fiber_no, motor_unit_no), # random value to add or substract to setSpecificStatesCallFrequency every stimulation, this is to add random jitter to the frequency
                  "setSpecificStatesRepeatAfterFirstCall":  0.01,                                           # simulation time span for which the setSpecificStates callback will be called after a call was triggered
                  "setSpecificStatesCallEnableBegin":       variables.get_specific_states_call_enable_begin(fiber_no, motor_unit_no),# [ms] first time when to call setSpecificStates
                  "additionalArgument":                     fiber_no,
                  #"handleResultFunction": handleResult,
                  #"handleResultCallInterval": 2e3,
                  
                  # parameters to the cellml model
                  "mappings":                               variables.mappings,
                  "parametersInitialValues":                variables.parameters_initial_values,
                  
                  "meshName":                               "MeshFiber_{}".format(fiber_no),
                  "stimulationLogFilename":                 "out/stimulation.log",
                },
              },
            } for fiber_in_subdomain_coordinate_y in range(n_fibers_in_subdomain_y(subdomain_coordinate_y)) \
                for fiber_in_subdomain_coordinate_x in range(n_fibers_in_subdomain_x(subdomain_coordinate_x)) \
                  for fiber_no in [get_fiber_no(subdomain_coordinate_x, subdomain_coordinate_y, fiber_in_subdomain_coordinate_x, fiber_in_subdomain_coordinate_y)] \
                    for motor_unit_no in [get_motor_unit_no(fiber_no)]],
          }
        },
        "Term2": {     # Diffusion
          "MultipleInstances": {
            "nInstances": n_fibers_in_subdomain_x(subdomain_coordinate_x)*n_fibers_in_subdomain_y(subdomain_coordinate_y),
            "instances": 
            [{
              "ranks":                         list(range(variables.n_subdomains_z)),    # these rank nos are local nos to the outer instance of MultipleInstances, i.e. from 0 to number of ranks in z direction
              "CrankNicolson" : {
                "initialValues":               [],
                #"numberTimeSteps":            1,
                "timeStepWidth":               variables.dt_1D,  # 1e-5
                "timeStepWidthRelativeTolerance": 1e-10,
                "logTimeStepWidthAsKey":       "dt_1D",
                "durationLogKey":              "duration_1D",
                "timeStepOutputInterval":      1e4,
                "timeStepWidthRelativeTolerance": 1e-10,
                "dirichletBoundaryConditions": {},                                       # old Dirichlet BC that are not used in FastMonodomainSolver: {0: -75.0036, -1: -75.0036},
                "dirichletOutputFilename":     None,                                 # filename for a vtp file that contains the Dirichlet boundary condition nodes and their values, set to None to disable
                "inputMeshIsGlobal":           True,
                "solverName":                  "implicitSolver",
                "nAdditionalFieldVariables":   1,     # for stress that will be transferred from CellML and then written with output writer
                "additionalSlotNames":         [],
                "checkForNanInf":              False,
                    
                "FiniteElementMethod" : {
                  "inputMeshIsGlobal":         True,
                  "meshName":                  "MeshFiber_{}".format(fiber_no),
                  "prefactor":                 get_diffusion_prefactor(fiber_no, motor_unit_no),  # resolves to Conductivity / (Am * Cm)
                  "solverName":                "implicitSolver",
                  "slotName":                  "",
                },
                "OutputWriter" : [
                  #{"format": "Paraview", "outputInterval": int(1./variables.dt_1D*variables.output_timestep), "filename": "out/fiber_"+str(fiber_no), "binary": True, "fixedFormat": False, "combineFiles": True},
                  #{"format": "Paraview", "outputInterval": 1./variables.dt_1D*variables.output_timestep, "filename": "out/fiber_"+str(i)+"_txt", "binary": False, "fixedFormat": False},
                  #{"format": "ExFile", "filename": "out/fiber_"+str(i), "outputInterval": 1./variables.dt_1D*variables.output_timestep, "sphereSize": "0.02*0.02*0.02"},
                  #{"format": "PythonFile", "filename": "out/fiber_"+str(i), "outputInterval": 1./variables.dt_1D*variables.output_timestep, "binary":True, "onlyNodalValues":True},
                ]
              },
            } for fiber_in_subdomain_coordinate_y in range(n_fibers_in_subdomain_y(subdomain_coordinate_y)) \
                for fiber_in_subdomain_coordinate_x in range(n_fibers_in_subdomain_x(subdomain_coordinate_x)) \
                  for fiber_no in [get_fiber_no(subdomain_coordinate_x, subdomain_coordinate_y, fiber_in_subdomain_coordinate_x, fiber_in_subdomain_coordinate_y)] \
                    for motor_unit_no in [get_motor_unit_no(fiber_no)]],
            "OutputWriter" : variables.output_writer_fibers,
          },
        },
      }
    } if (subdomain_coordinate_x,subdomain_coordinate_y) == (variables.own_subdomain_coordinate_x,variables.own_subdomain_coordinate_y) else None
    for subdomain_coordinate_y in range(variables.n_subdomains_y)
        for subdomain_coordinate_x in range(variables.n_subdomains_x)]
  },
}

# add entry for when fast_fibers with "RepeatedCall" as top solver is used
config["RepeatedCall"] = {
  "timeStepWidth":          variables.output_timestep,  # 1e-1
  "logTimeStepWidthAsKey":  "dt_output_timestep",
  "durationLogKey":         "duration_repeated_call",
  "timeStepOutputInterval": 1,
  "endTime":                variables.end_time,
  "MultipleInstances": copy.deepcopy(config["MultipleInstances"]),
  "fiberDistributionFile":    variables.fiber_distribution_file,   # for FastMonodomainSolver, e.g. MU_fibre_distribution_3780.txt
  "firingTimesFile":          variables.firing_times_file,         # for FastMonodomainSolver, e.g. MU_firing_times_real.txt
  "onlyComputeIfHasBeenStimulated": True,                          # only compute fibers after they have been stimulated for the first time
  "disableComputationWhenStatesAreCloseToEquilibrium": True,       # optimization where states that are close to their equilibrium will not be computed again
  "valueForStimulatedPoint":  variables.vm_value_stimulated,       # to which value of Vm the stimulated node should be set
}

# loop over instances (fibers)
for i in range(len(config["RepeatedCall"]["MultipleInstances"]["instances"])):
  #config["RepeatedCall"]["MultipleInstances"]["instances"][i]["StrangSplitting"]["endTime"] = variables.output_timestep
  
  if config["RepeatedCall"]["MultipleInstances"]["instances"][i] is None:
    continue
    
  # loop over output writers
  if config["RepeatedCall"]["MultipleInstances"]["instances"][i]:
    for j in range(len(config["RepeatedCall"]["MultipleInstances"]["instances"][i]["StrangSplitting"]["Term2"]["MultipleInstances"]["OutputWriter"])):
      
      # set outputInterval to 1
      config["RepeatedCall"]["MultipleInstances"]["instances"][i]["StrangSplitting"]["Term2"]["MultipleInstances"]["OutputWriter"][j]["outputInterval"] = 1

#print(config["RepeatedCall"]["MultipleInstances"]["instances"][0]["StrangSplitting"]["Term1"]["MultipleInstances"]["instances"][0]["Heun"]["CellML"]);

# stop timer and calculate how long parsing lasted
if rank_no == 0:
  t_stop_script = timeit.default_timer()
  print("Python config parsed in {:.1f}s.".format(t_stop_script - t_start_script))
