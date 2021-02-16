# Multiple 1D fibers (monodomain) with 3D intra-muscular EMG (static bidomain) and 3D fat layer (anisotropic diffusion), on biceps geometry

import sys, os
import timeit
import argparse
import importlib

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
else:
  print("Error, please specify variables file, e.g:\n ./static_biceps_emg ../settings_static_biceps_emg.py ramp.py")
  exit(0)

# define command line arguments
parser = argparse.ArgumentParser(description='static_biceps_emg')
parser.add_argument('--scenario_name',                       help='The name to identify this run in the log.',   default=variables.scenario_name)
parser.add_argument('--n_subdomains', nargs=3,               help='Number of subdomains in x,y,z direction.',    type=int)
parser.add_argument('--n_subdomains_x', '-x',                help='Number of subdomains in x direction.',        type=int, default=variables.n_subdomains_x)
parser.add_argument('--n_subdomains_y', '-y',                help='Number of subdomains in y direction.',        type=int, default=variables.n_subdomains_y)
parser.add_argument('--n_subdomains_z', '-z',                help='Number of subdomains in z direction.',        type=int, default=variables.n_subdomains_z)
parser.add_argument('--diffusion_solver_type',               help='The solver for the diffusion.',               default=variables.diffusion_solver_type, choices=["gmres","cg","lu","gamg","richardson","chebyshev","cholesky","jacobi","sor","preonly"])
parser.add_argument('--diffusion_preconditioner_type',       help='The preconditioner for the diffusion.',       default=variables.diffusion_preconditioner_type, choices=["jacobi","sor","lu","ilu","gamg","none"])
parser.add_argument('--potential_flow_solver_type',          help='The solver for the potential flow (non-spd matrix).', default=variables.potential_flow_solver_type, choices=["gmres","cg","lu","gamg","richardson","chebyshev","cholesky","jacobi","sor","preonly"])
parser.add_argument('--potential_flow_preconditioner_type',  help='The preconditioner for the potential flow.',  default=variables.potential_flow_preconditioner_type, choices=["jacobi","sor","lu","ilu","gamg","none"])
parser.add_argument('--emg_solver_type',                     help='The solver for the static bidomain.',         default=variables.emg_solver_type)
#parser.add_argument('--emg_solver_type',                    help='The solver for the static bidomain.',         default=variables.emg_solver_type, choices=["gmres","cg","lu","gamg","richardson","chebyshev","cholesky","jacobi","sor","preonly"])
parser.add_argument('--emg_preconditioner_type',             help='The preconditioner for the static bidomain.', default=variables.emg_preconditioner_type, choices=["jacobi","sor","lu","ilu","gamg","none"])
parser.add_argument('--emg_initial_guess_nonzero',           help='If the initial guess for the emg linear system should be set to the previous solution.', default=variables.emg_initial_guess_nonzero, action='store_true')
parser.add_argument('--paraview_output',                     help='Enable the paraview output writer.',          default=variables.paraview_output, action='store_true')
parser.add_argument('--adios_output',                        help='Enable the MegaMol/ADIOS output writer.',          default=variables.adios_output, action='store_true')
parser.add_argument('--fiber_file',                          help='The filename of the file that contains the fiber data.', default=variables.fiber_file)
parser.add_argument('--fiber_distribution_file',             help='The filename of the file that contains the MU firing times.', default=variables.fiber_distribution_file)
parser.add_argument('--firing_times_file',                   help='The filename of the file that contains the cellml model.', default=variables.firing_times_file)
parser.add_argument('--end_time', '--tend', '-t',            help='The end simulation time.',                    type=float, default=variables.end_time)
parser.add_argument('--output_timestep',                     help='The timestep for writing outputs.',           type=float, default=variables.output_timestep)
parser.add_argument('--dt_0D',                               help='The timestep for the 0D model.',              type=float, default=variables.dt_0D)
parser.add_argument('--dt_1D',                               help='The timestep for the 1D model.',              type=float, default=variables.dt_1D)
parser.add_argument('--dt_splitting',                        help='The timestep for the splitting.',             type=float, default=variables.dt_splitting)
parser.add_argument('--dt_3D',                               help='The timestep for the 3D model, either bidomain or mechanics.', type=float, default=variables.dt_3D)
parser.add_argument('--disable_firing_output',               help='Disables the initial list of fiber firings.', default=variables.disable_firing_output, action='store_true')
parser.add_argument('--v',                                   help='Enable full verbosity in c++ code')
parser.add_argument('-v',                                    help='Enable verbosity level in c++ code', action="store_true")
parser.add_argument('-vmodule',                              help='Enable verbosity level for given file in c++ code')
parser.add_argument('-pause',                                help='Stop at parallel debugging barrier', action="store_true")

# parse command line arguments and assign values to variables module
args, other_args = parser.parse_known_args(args=sys.argv[:-2], namespace=variables)
if len(other_args) != 0 and rank_no == 0:
    print("Warning: These arguments were not parsed by the settings python file\n  " + "\n  ".join(other_args), file=sys.stderr)

# initialize some dependend variables
if variables.n_subdomains is not None:
  variables.n_subdomains_x = variables.n_subdomains[0]
  variables.n_subdomains_y = variables.n_subdomains[1]
  variables.n_subdomains_z = variables.n_subdomains[2]
  
# output information of run
if rank_no == 0:
  print("scenario_name: {},  n_subdomains: {} {} {},  n_ranks: {},  end_time: {}".format(variables.scenario_name, variables.n_subdomains_x, variables.n_subdomains_y, variables.n_subdomains_z, n_ranks, variables.end_time))
  print("dt_0D:           {:0.0e}, diffusion_solver_type:      {}".format(variables.dt_0D, variables.diffusion_solver_type))
  print("dt_1D:           {:0.0e}, potential_flow_solver_type: {}".format(variables.dt_1D, variables.potential_flow_solver_type))
  print("dt_splitting:    {:0.0e}, emg_solver_type:            {}, emg_initial_guess_nonzero: {}".format(variables.dt_splitting, variables.emg_solver_type, variables.emg_initial_guess_nonzero))
  print("dt_3D:           {:0.0e}, paraview_output: {}".format(variables.dt_3D, variables.paraview_output))
  print("output_timestep: {:0.0e}  stimulation_frequency: {} 1/ms = {} Hz".format(variables.output_timestep, variables.stimulation_frequency, variables.stimulation_frequency*1e3))
  print("fiber_file:              {}".format(variables.fiber_file))
  print("fat_mesh_file:           {}".format(variables.fat_mesh_file))
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
  "scenarioName":          variables.scenario_name,
  "logFormat":             "csv",
  "Meshes":                variables.meshes,
  "MappingsBetweenMeshes": variables.mappings_between_meshes,
  "Solvers": {
    "diffusionTermSolver": {# solver for the implicit timestepping scheme of the diffusion time step
      "maxIterations":      1e4,
      "relativeTolerance":  1e-10,
      "absoluteTolerance":  1e-10,         # 1e-10 absolute tolerance of the residual    
      "solverType":         variables.diffusion_solver_type,
      "preconditionerType": variables.diffusion_preconditioner_type,
      "dumpFilename":       "",   # "out/dump_"
      "dumpFormat":         "matlab",
    },
    "potentialFlowSolver": {# solver for the initial potential flow, that is needed to estimate fiber directions for the bidomain equation
      "relativeTolerance":  1e-10,
      "absoluteTolerance":  1e-10,         # 1e-10 absolute tolerance of the residual    
      "maxIterations":      1e4,
      "solverType":         variables.potential_flow_solver_type,
      "preconditionerType": variables.potential_flow_preconditioner_type,
      "dumpFilename":       "",
      "dumpFormat":         "matlab",
    },
    "muscularEMGSolver": {   # solver for the static Bidomain equation and the EMG
      "relativeTolerance":  1e-5,
      "absoluteTolerance":  1e-10,         # 1e-10 absolute tolerance of the residual    
      "maxIterations":      1e4,
      "solverType":         variables.emg_solver_type,
      "preconditionerType": variables.emg_preconditioner_type,
      "dumpFilename":       "",
      "dumpFormat":         "matlab",
    },
    "fatEMGSolver": {       # solver for the Laplace equation, conduction through fat layer
      "relativeTolerance":  1e-5,
      "absoluteTolerance":  1e-10,         # 1e-10 absolute tolerance of the residual    
      "maxIterations":      1e4,
      "solverType":         "lu",
      "preconditionerType": "none",
      #"dumpFilename":       "out/ramp/fat_emg",
      "dumpFilename":       "",
      "dumpFormat":         "matlab",
    },
  },
  "Coupling": {
    "timeStepWidth":          variables.dt_3D,  # 1e-1
    "logTimeStepWidthAsKey":  "dt_3D",
    "durationLogKey":         "duration_total",
    "timeStepOutputInterval": 1,
    "endTime":                variables.end_time,
    "connectedSlotsTerm1To2": [0],          # transfer Vm to StaticBidomainSolver
    "connectedSlotsTerm2To1": [None],       # transfer nothing back
    "Term1": {        # monodomain, fibers
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
            "endTime":                variables.dt_splitting,
            "connectedSlotsTerm1To2": [0],   # transfer slot 0 = state Vm from Term1 (CellML) to Term2 (Diffusion)
            "connectedSlotsTerm2To1": [0],   # transfer the same back, this avoids data copy

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
                    "dirichletBoundaryConditions":  {},
                    "dirichletOutputFilename":      None,             # filename for a vtp file that contains the Dirichlet boundary condition nodes and their values, set to None to disable
                    "nAdditionalFieldVariables":    0,
                      
                    "CellML" : {
                      "modelFilename":                         variables.cellml_file,                          # input C++ source file, can be either generated by OpenCMISS or OpenCOR from cellml model
                      "compilerFlags":                          "-fPIC -O3 -shared ",
                      #"simdSourceFilename" :                   "simdcode.cpp",                                 # transformed C++ source file that gets generated from modelFilename and is ready for multiple instances
                      #"libraryFilename":                       "cellml_simd_lib.so",                           # compiled library
                      "useGivenLibrary":                        False,
                      #"statesInitialValues":                   [],
                      #"setSpecificParametersFunction":         set_specific_parameters,                        # callback function that sets parameters like stimulation current
                      #"setSpecificParametersCallInterval":     int(1./variables.stimulation_frequency/variables.dt_0D),         # set_specific_parameters should be called every 0.1, 5e-5 * 1e3 = 5e-2 = 0.05
                      "setSpecificStatesFunction":              set_specific_states,                                             # callback function that sets states like Vm, activation can be implemented by using this method and directly setting Vm values, or by using setParameters/setSpecificParameters
                      #"setSpecificStatesCallInterval":         2*int(1./variables.stimulation_frequency/variables.dt_0D),       # set_specific_states should be called variables.stimulation_frequency times per ms, the factor 2 is needed because every Heun step includes two calls to rhs
                      "setSpecificStatesCallInterval":          0,                                                               # 0 means disabled
                      "setSpecificStatesCallFrequency":         variables.get_specific_states_call_frequency(fiber_no, motor_unit_no),   # set_specific_states should be called variables.stimulation_frequency times per ms
                      "setSpecificStatesFrequencyJitter":       variables.get_specific_states_frequency_jitter(fiber_no, motor_unit_no), # random value to add or substract to setSpecificStatesCallFrequency every stimulation, this is to add random jitter to the frequency
                      "setSpecificStatesRepeatAfterFirstCall":  0.01,                                                            # [ms] simulation time span for which the setSpecificStates callback will be called after a call was triggered
                      "setSpecificStatesCallEnableBegin":       variables.get_specific_states_call_enable_begin(fiber_no, motor_unit_no),# [ms] first time when to call setSpecificStates
                      "additionalArgument":                     fiber_no,
                      "algebraicsForTransfer":               variables.output_algebraic_index,            # which algebraic values to use in further computation
                      "statesForTransfer":                      variables.output_state_index,                   # which state values to use in further computation, Shorten / Hodgkin Huxley: state 0 = Vm
                      
                      "parametersUsedAsAlgebraic":           variables.parameters_used_as_algebraic,      #[32],       # list of algebraic value indices, that will be set by parameters. Explicitely defined parameters that will be copied to algebraics, this vector contains the indices of the algebraic array. This is ignored if the input is generated from OpenCMISS generated c code.
                      "parametersUsedAsConstant":               variables.parameters_used_as_constant,          #[65],           # list of constant value indices, that will be set by parameters. This is ignored if the input is generated from OpenCMISS generated c code.
                      "parametersInitialValues":                variables.parameters_initial_values,            #[0.0, 1.0],      # initial values for the parameters: I_Stim, l_hs
                      "meshName":                               "MeshFiber_{}".format(fiber_no),
                      "stimulationLogFilename":                 "out/stimulation.log",
                    },      
                    "OutputWriter" : [
                      {"format": "Paraview", "outputInterval": 1, "filename": "out/" + variables.scenario_name + "/0D_states({},{})".format(fiber_in_subdomain_coordinate_x,fiber_in_subdomain_coordinate_y), "binary": True, "fixedFormat": False, "combineFiles": True}
                    ] if variables.states_output else []
                    
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
                  "ImplicitEuler" : {
                    "initialValues":               [],
                    #"numberTimeSteps":            1,
                    "timeStepWidth":               variables.dt_1D,  # 1e-5
                    "timeStepWidthRelativeTolerance": 1e-10,
                    "logTimeStepWidthAsKey":       "dt_1D",
                    "durationLogKey":              "duration_1D",
                    "timeStepOutputInterval":      1e4,
                    "dirichletBoundaryConditions": {},            # old Dirichlet BC that are not used in FastMonodomainSolver: {0: -75.0036, -1: -75.0036},
                    "dirichletOutputFilename":     None,          # filename for a vtp file that contains the Dirichlet boundary condition nodes and their values, set to None to disable
                    "inputMeshIsGlobal":           True,
                    "solverName":                  "diffusionTermSolver",
                    "nAdditionalFieldVariables":   0,
                    "FiniteElementMethod" : {
                      "inputMeshIsGlobal":         True,
                      "meshName":                  "MeshFiber_{}".format(fiber_no),
                      "prefactor":                 get_diffusion_prefactor(fiber_no, motor_unit_no),  # resolves to Conductivity / (Am * Cm)
                      "solverName":                "diffusionTermSolver",
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
      "fiberDistributionFile":    variables.fiber_distribution_file,   # for FastMonodomainSolver, e.g. MU_fibre_distribution_3780.txt
      "firingTimesFile":          variables.firing_times_file,         # for FastMonodomainSolver, e.g. MU_firing_times_real.txt
    },
    "Term2": {        # Bidomain, EMG
      "Coupling": {
        "numberTimeSteps":        1,
        "logTimeStepWidthAsKey":  "dt_3D",
        "durationLogKey":         "duration_emg",
        "timeStepOutputInterval": 1,
        "endTime":                1,                   # not relevant here, will be set by the outer coupling
        "connectedSlotsTerm1To2": [0],                 # transfer Vm from StaticBidomainSolver to fat diffusion
        "connectedSlotsTerm2To1": [None],              # transfer nothing back
        "Term1": {        # bidomain
          "OutputSurface": {        # version for fibers_emg_2d_output
            "OutputWriter": [
              {"format": "Paraview", "outputInterval": int(1./variables.dt_3D*variables.output_timestep_smaller_files), "filename": "out/" + variables.scenario_name + "/surface_emg", "binary": True, "fixedFormat": False, "combineFiles": True},
            ],
            "face": "0+",
            "StaticBidomainSolver": {
              "numberTimeSteps":        1,
              "timeStepOutputInterval": 50,
              "durationLogKey":         "duration_bidomain",
              "solverName":             "muscularEMGSolver",
              "initialGuessNonzero":    variables.emg_initial_guess_nonzero,
              "PotentialFlow": {
                "FiniteElementMethod" : {  
                  "meshName":           "3Dmesh",
                  "solverName":         "potentialFlowSolver",
                  "prefactor":          1.0,
                  "dirichletBoundaryConditions": variables.potential_flow_dirichlet_bc,
                  "dirichletOutputFilename":     None,                # filename for a vtp file that contains the Dirichlet boundary condition nodes and their values, set to None to disable
                  "neumannBoundaryConditions":   [],
                  "inputMeshIsGlobal":  True,
                },
              },
              "Activation": {
                "FiniteElementMethod" : {  
                  "meshName":           "3Dmesh",
                  "solverName":         "muscularEMGSolver",
                  "prefactor":          1.0,
                  "inputMeshIsGlobal":  True,
                  "dirichletBoundaryConditions": {},
                  "dirichletOutputFilename":     None,                # filename for a vtp file that contains the Dirichlet boundary condition nodes and their values, set to None to disable
                  "neumannBoundaryConditions":   [],
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
              "OutputWriter" : variables.output_writer_emg,
            }
          },
        },
        "Term2": {
          "MultipleInstances": {
            "logKey":                     "duration_fat",
            "ranksAllComputedInstances":  variables.fat_global_rank_nos,
            "nInstances":                 1,
            "instances": 
            [{
              "ranks":                    list(range(len(variables.fat_global_rank_nos))),
              "RepeatedCallStatic": {     # RepeatedCallStatic is only needed to call a static solver within a timestepping scheme
                "numberTimeSteps":        1,    # 1 call to the static solver per timestep
                "timeStepOutputInterval": 10,   # >1: do not call
                "FiniteElementMethod": {
                  "meshName":             "3DFatMesh",
                  "solverName":           "fatEMGSolver",
                  "prefactor":            0.4,
                  "dirichletBoundaryConditions": variables.fat_dirichlet_bc,
                  "dirichletOutputFilename":     None,            # filename for a vtp file that contains the Dirichlet boundary condition nodes and their values, set to None to disable
                  "updatePrescribedValuesFromSolution": True,     # update the prescribed Dirichlet boundary condition values at the beginning of each timestep by the values that got transferred from the intra-muscular domain by the coupling.
                  "neumannBoundaryConditions": [],
                  "inputMeshIsGlobal":  False,
                },
              },
            }],
            "OutputWriter": [
              {"format": "Paraview", "outputInterval": int(1./variables.dt_3D*variables.output_timestep_smaller_files), "filename": "out/" + variables.scenario_name + "/hd_emg_fat", "binary": True, "fixedFormat": False, "combineFiles": True},
            ]  
          },
        },
      },
    }
  }
}

# stop timer and calculate how long parsing lasted
if rank_no == 0:
  t_stop_script = timeit.default_timer()
  print("Python config parsed in {:.1f}s.".format(t_stop_script - t_start_script))
