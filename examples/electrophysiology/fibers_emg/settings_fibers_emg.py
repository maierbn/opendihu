# Multiple 1D variables.fibers (monodomain) with 3D EMG (static bidomain), biceps geometry
# to see all available arguments, execute: ./fibers_emg ../settings_fibers_emg.py -help
##
# if variables.fiber_file=cuboid.bin, it uses a small cuboid test example
#
# You have to set n_subdomains such that it matches the number of processes, e.g. 2x2x1 = 4 processes.
# Decomposition is in x,y,z direction, the fibers are aligned with the z axis.
# E.g. --n_subdomains 2 2 1 which is 2x2x1 means no subdivision per fiber, 
# --n_subdomains 8 8 4 means every fiber will be subdivided to 4 processes and all fibers will be computed by 8x8 processes.
#
# Example with 4 processes and end time 5:
#   mpirun -n 4 ./fibers_emg ../settings_fibers_emg.py --n_subdomains 2 2 1 --variables.end_time=5.0
#
# Three files contribute to the settings:
# A lot of variables are set by the helper.py script, the variables and their values are defined in variables.py and this file
# creates the composite config that is needed by opendihu.
# You should only make changes in the variables.py file to set different parameters or in this file to add functionality.
# This is the most complex example, you can ask me (BM) how it works exactly.


import sys
import timeit
import argparse

sys.path.insert(0, '..')
import variables              # file variables.py, defined default values for all parameters, you can set the parameters there
from create_partitioned_meshes_for_settings import *   # file create_partitioned_meshes_for_settings with helper functions about own subdomain

# define command line arguments
parser = argparse.ArgumentParser(description='fibers_emg')
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
parser.add_argument('--fiber_file',                          help='The filename of the file that contains the fiber data.', default=variables.fiber_file)
parser.add_argument('--cellml_file',                         help='The filename of the file that contains the cellml model.', default=variables.cellml_file)
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
parser.add_argument('--rank_reordering',                     help='Enable rank reordering in the c++ code', action="store_true")
parser.add_argument('--linear_elasticity',                   help='Enable linear elasticity', action="store_true")
 
# parse command line arguments and assign values to variables module
args = parser.parse_args(args=sys.argv[:-2], namespace=variables)

# initialize some dependend variables
if variables.n_subdomains is not None:
  variables.n_subdomains_x = variables.n_subdomains[0]
  variables.n_subdomains_y = variables.n_subdomains[1]
  variables.n_subdomains_z = variables.n_subdomains[2]
  
if variables.linear_elasticity:
  variables.cellml_file = "../../input/shorten.cpp"
  variables.emg_solver_type = "cg"

# parse arguments
rank_no = (int)(sys.argv[-2])
n_ranks = (int)(sys.argv[-1])
  
# output information of run
if rank_no == 0:
  print("scenario_name: {},  n_subdomains: {} {} {},  n_ranks: {},  end_time: {}".format(variables.scenario_name, variables.n_subdomains_x, variables.n_subdomains_y, variables.n_subdomains_z, n_ranks, variables.end_time))
  print("dt_0D:           {:0.0e}, diffusion_solver_type:      {}".format(variables.dt_0D, variables.diffusion_solver_type))
  print("dt_1D:           {:0.0e}, potential_flow_solver_type: {}".format(variables.dt_1D, variables.potential_flow_solver_type))
  print("dt_splitting:    {:0.0e}, emg_solver_type:            {}, emg_initial_guess_nonzero: {}".format(variables.dt_splitting, variables.emg_solver_type, variables.emg_initial_guess_nonzero))
  print("dt_3D:           {:0.0e}, paraview_output: {}".format(variables.dt_3D, variables.paraview_output))
  print("output_timestep: {:0.0e}  stimulation_frequency: {} 1/ms = {} Hz".format(variables.output_timestep, variables.stimulation_frequency, variables.stimulation_frequency*1e3))
  print("fiber_file:              {}".format(variables.fiber_file))
  print("cellml_file:             {}".format(variables.cellml_file))
  print("fiber_distribution_file: {}".format(variables.fiber_distribution_file))
  print("firing_times_file:       {}".format(variables.firing_times_file))
  print("********************************************************************************")
  
  # start timer to measure duration of parsing of this script  
  t_start_script = timeit.default_timer()
    
# initialize all helper variables
from helper import *

variables.n_subdomains_xy = variables.n_subdomains_x * variables.n_subdomains_y
variables.n_fibers_total = variables.n_fibers_x * variables.n_fibers_y

# define the config dict
config = {
  "scenarioName": variables.scenario_name,
  "Meshes": variables.meshes,
  "MappingsBetweenMeshes": {"MeshFiber_{}".format(i) : "3Dmesh" for i in range(variables.n_fibers_total)},
  "Solvers": {
    "implicitSolver": {     # solver for the implicit timestepping scheme of the diffusion time step
      "maxIterations":      1e4,
      "relativeTolerance":  1e-10,
      "solverType":         variables.diffusion_solver_type,
      "preconditionerType": variables.diffusion_preconditioner_type,
      "dumpFilename":       "",   # "out/dump_"
      "dumpFormat":         "matlab",
    },
    "potentialFlowSolver": {# solver for the initial potential flow, that is needed to estimate fiber directions for the bidomain equation
      "relativeTolerance":  1e-10,
      "maxIterations":      1e4,
      "solverType":         variables.potential_flow_solver_type,
      "preconditionerType": variables.potential_flow_preconditioner_type,
      "dumpFilename":       "",
      "dumpFormat":         "matlab",
    },
    "activationSolver": {   # solver for the static Bidomain equation and the EMG
      "relativeTolerance":  1e-100,
      "maxIterations":      1e4,
      "solverType":         variables.emg_solver_type,
      "preconditionerType": variables.emg_preconditioner_type,
      "dumpFilename":       "",
      "dumpFormat":         "matlab",
    },
    "linearElasticitySolver": {   # solver for linear elasticity
      "relativeTolerance":  1e-1,
      "maxIterations":      1e4,
#      "solverType":        variables.emg_solver_type,
#      "preconditionerType": variables.emg_preconditioner_type,
      "dumpFilename":       "",
      "dumpFormat":         "matlab",
    }, 
  },
  "Coupling": {
    "timeStepWidth":          variables.dt_3D,  # 1e-1
    "logTimeStepWidthAsKey":  "dt_3D",
    "durationLogKey":         "duration_total",
    "timeStepOutputInterval": 10,
    "endTime":                variables.end_time,
    "transferSlotName":       "intermediates" if variables.linear_elasticity else "states",
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
            "transferSlotName":       "states",   # which output slot of the cellml adapter ("states" or "intermediates") to use for transfer to diffusion, in this case we need "states", states[0] which is Vm

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
                      
                    "CellML" : {
                      "sourceFilename":                         variables.cellml_file,                          # input C++ source file, can be either generated by OpenCMISS or OpenCOR from cellml model
                      "compilerFlags":                          "-fPIC -O3 -shared ",
                      #"simdSourceFilename" :                   "simdcode.cpp",                                 # transformed C++ source file that gets generated from sourceFilename and is ready for multiple instances
                      #"libraryFilename":                       "cellml_simd_lib.so",                           # compiled library
                      "useGivenLibrary":                        False,
                      #"statesInitialValues":                   [],
                      #"setSpecificParametersFunction":         set_specific_parameters,                        # callback function that sets parameters like stimulation current
                      #"setSpecificParametersCallInterval":     int(1./variables.stimulation_frequency/variables.dt_0D),     # set_parameters should be called every 0.1, 5e-5 * 1e3 = 5e-2 = 0.05
                      "setSpecificStatesFunction":              set_specific_states,                            # callback function that sets states like Vm, activation can be implemented by using this method and directly setting Vm values, or by using setParameters/setSpecificParameters
                      #"setSpecificStatesCallInterval":         2*int(1./variables.stimulation_frequency/variables.dt_0D),     # set_specific_states should be called variables.stimulation_frequency times per ms, the factor 2 is needed because every Heun step includes two calls to rhs
                      "setSpecificStatesCallInterval":          0,                                              # 0 means disabled
                      "setSpecificStatesCallFrequency":         variables.stimulation_frequency,                # set_specific_states should be called variables.stimulation_frequency times per ms, the factor 2 is needed because every Heun step includes two calls to rhs
                      "setSpecificStatesRepeatAfterFirstCall":  0.01,                                           # simulation time span for which the setSpecificStates callback will be called after a call was triggered
                      "additionalArgument":                     fiber_no(subdomain_coordinate_x, subdomain_coordinate_y, fiber_in_subdomain_coordinate_x, fiber_in_subdomain_coordinate_y),
                      
                      "outputIntermediateIndex":                0,                                              # which intermediate value to use in further computation
                      "outputStateIndex":                       0,                                              # Shorten / Hodgkin Huxley: state 0 = Vm, Shorten: rate 28 = gamma, intermediate 0 = gamma
                      "parametersUsedAsIntermediate":           variables.parameters_used_as_intermediate,      #[32],       # list of intermediate value indices, that will be set by parameters. Explicitely defined parameters that will be copied to intermediates, this vector contains the indices of the algebraic array. This is ignored if the input is generated from OpenCMISS generated c code.
                      "parametersUsedAsConstant":               variables.parameters_used_as_constant,          #[65],           # list of constant value indices, that will be set by parameters. This is ignored if the input is generated from OpenCMISS generated c code.
                      "parametersInitialValues":                variables.parameters_initial_values,            #[0.0, 1.0],      # initial values for the parameters: I_Stim, l_hs
                      "meshName":                               "MeshFiber_{}".format(fiber_no(subdomain_coordinate_x, subdomain_coordinate_y, fiber_in_subdomain_coordinate_x, fiber_in_subdomain_coordinate_y)),
                      "prefactor":                              1.0,
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
                  "ranks":                         list(range(variables.n_subdomains_z)),    # these rank nos are local nos to the outer instance of MultipleInstances, i.e. from 0 to number of ranks in z direction
                  "ImplicitEuler" : {
                    "initialValues":               [],
                    #"numberTimeSteps":            1,
                    "timeStepWidth":               variables.dt_1D,  # 1e-5
                    "logTimeStepWidthAsKey":       "dt_1D",
                    "durationLogKey":              "duration_1D",
                    "timeStepOutputInterval":      1e4,
                    "dirichletBoundaryConditions": {0: -75.0036, -1: -75.0036},
                    "inputMeshIsGlobal":           True,
                    "solverName":                  "implicitSolver",
                    "FiniteElementMethod" : {
                      "maxIterations":             1e4,
                      "relativeTolerance":         1e-10,
                      "inputMeshIsGlobal":         True,
                      "meshName":                  "MeshFiber_{}".format(fiber_no(subdomain_coordinate_x, subdomain_coordinate_y, fiber_in_subdomain_coordinate_x, fiber_in_subdomain_coordinate_y)),
                      "prefactor":                 variables.Conductivity/(variables.Am*variables.Cm),
                      "solverName":                "implicitSolver",
                    },
                    "OutputWriter" : [
                      #{"format": "Paraview", "outputInterval": int(1./variables.dt_1D*variables.output_timestep), "filename": "out/fiber_"+str(fiber_no(subdomain_coordinate_x, subdomain_coordinate_y, fiber_in_subdomain_coordinate_x, fiber_in_subdomain_coordinate_y)), "binary": True, "fixedFormat": False, "combineFiles": True},
                      #{"format": "Paraview", "outputInterval": 1./variables.dt_1D*variables.output_timestep, "filename": "out/fiber_"+str(i)+"_txt", "binary": False, "fixedFormat": False},
                      #{"format": "ExFile", "filename": "out/fiber_"+str(i), "outputInterval": 1./variables.dt_1D*variables.output_timestep, "sphereSize": "0.02*0.02*0.02"},
                      #{"format": "PythonFile", "filename": "out/fiber_"+str(i), "outputInterval": 1./variables.dt_1D*variables.output_timestep, "binary":True, "onlyNodalValues":True},
                    ]
                  },
                } for fiber_in_subdomain_coordinate_y in range(n_fibers_in_subdomain_y(subdomain_coordinate_y)) \
                    for fiber_in_subdomain_coordinate_x in range(n_fibers_in_subdomain_x(subdomain_coordinate_x))],
                "OutputWriter" : variables.output_writer_fibers,
              },
            },
          }
        } if (subdomain_coordinate_x,subdomain_coordinate_y) == (variables.own_subdomain_coordinate_x,variables.own_subdomain_coordinate_y) else None
        for subdomain_coordinate_y in range(variables.n_subdomains_y)
            for subdomain_coordinate_x in range(variables.n_subdomains_x)]
      }
    },
    "Term2": {        # Bidomain, EMG
      "StaticBidomainSolver": {       # version for fibers_emg
        "timeStepWidth":          variables.dt_3D,
        "timeStepOutputInterval": 50,
        "durationLogKey":         "duration_bidomain",
        "solverName":             "activationSolver",
        "initialGuessNonzero":    variables.emg_initial_guess_nonzero,
        "PotentialFlow": {
          "FiniteElementMethod" : {  
            "meshName":           "3Dmesh",
            "solverName":         "potentialFlowSolver",
            "prefactor":          1.0,
            "dirichletBoundaryConditions": variables.potential_flow_dirichlet_bc,
            "neumannBoundaryConditions":   [],
            "inputMeshIsGlobal":  True,
          },
        },
        "Activation": {
          "FiniteElementMethod" : {  
            "meshName":           "3Dmesh",
            "solverName":         "activationSolver",
            "prefactor":          1.0,
            "inputMeshIsGlobal":  True,
            "dirichletBoundaryConditions": {},
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
      },
      "OutputSurface": {        # version for fibers_emg_2d_output
        "OutputWriter": [
          {"format": "Paraview", "outputInterval": int(1./variables.dt_3D*variables.output_timestep), "filename": "out/" + variables.scenario_name + "/surface_emg", "binary": True, "fixedFormat": False, "combineFiles": True},
        ],
        "face": "1-",
        "StaticBidomainSolver": {
          "timeStepWidth":          variables.dt_3D,
          "timeStepOutputInterval": 50,
          "durationLogKey":         "duration_bidomain",
          "solverName":             "activationSolver",
          "initialGuessNonzero":    variables.emg_initial_guess_nonzero,
          "PotentialFlow": {
            "FiniteElementMethod" : {  
              "meshName":           "3Dmesh",
              "solverName":         "potentialFlowSolver",
              "prefactor":          1.0,
              "dirichletBoundaryConditions": variables.potential_flow_dirichlet_bc,
              "neumannBoundaryConditions":   [],
              "inputMeshIsGlobal":  True,
            },
          },
          "Activation": {
            "FiniteElementMethod" : {  
              "meshName":           "3Dmesh",
              "solverName":         "activationSolver",
              "prefactor":          1.0,
              "inputMeshIsGlobal":  True,
              "dirichletBoundaryConditions": {},
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
      "QuasiStaticLinearElasticitySolver": {
        "PotentialFlow": {        # potential flow for fiber directions in the 3D mesh
          "FiniteElementMethod" : {  
            "meshName":           "3Dmesh",
            "solverName":         "potentialFlowSolver",
            "prefactor":          1.0,
            "dirichletBoundaryConditions": variables.potential_flow_dirichlet_bc,
            "neumannBoundaryConditions":   [],
            "inputMeshIsGlobal":  True,
          },
        },
        "FiniteElementMethod" : {   # linear elasticity finite element method
          "meshName":             "3Dmesh",
          "solverName":           "linearElasticitySolver",
          "prefactor":            1.0,
          "inputMeshIsGlobal":    True,
          "dirichletBoundaryConditions": variables.linear_elasticity_dirichlet_bc,
          "neumannBoundaryConditions":   variables.linear_elasticity_neumann_bc,
          "bulkModulus":          4.0, #40e3 # https://www.researchgate.net/publication/230248067_Bulk_Modulus
          "shearModulus":         3.9, #39e3 # https://onlinelibrary.wiley.com/doi/full/10.1002/mus.24104
        },
        "maximumActiveStress":      1.0,
        "strainScalingCurveWidth":  1.0,
        "scalingFactor":        1e0, 
        "OutputWriter" : [
          {"format": "Paraview", "outputInterval": int(1./variables.dt_3D*variables.output_timestep), "filename": "out/deformation", "binary": True, "fixedFormat": False, "onlyNodalValues":True, "combineFiles":True},
          #{"format": "PythonFile", "filename": "out/deformation", "outputInterval": 1, "binary":False, "onlyNodalValues":True},
        ]
      }
    }
  }
}

# stop timer and calculate how long parsing lasted
if rank_no == 0:
  t_stop_script = timeit.default_timer()
  print("Python config parsed in {:.1f}s.".format(t_stop_script - t_start_script))
