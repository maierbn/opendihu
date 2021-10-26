# biceps with fibres, with previous contraction to obtain a smaller geometry and a stretching operation on this new geometry

import numpy as np
import scipy.stats
import pickle
import sys,os 
import timeit
import argparse
import importlib
import struct
sys.path.insert(0, "..")

# own MPI rank no and number of MPI ranks
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
  if rank_no == 0:
    print("Warning: There is no variables file, e.g:\n ./spindles_fibers ../settings_spindles_fibers.py fibers.py\n")
  exit(0)
  
# -------------------------------------------------------- begin parameters ---------------------------------------------------------
# -------------------------------------------------------- end parameters ------------------------------------------------------------
parser = argparse.ArgumentParser(description='spindles_fibers')
parser.add_argument('--scenario_name',                       help='The name to identify this run in the log.',   default=variables.scenario_name)
parser.add_argument('--n_subdomains', nargs=3,               help='Number of subdomains in x,y,z direction.',    type=int)
parser.add_argument('--n_subdomains_x', '-x',                help='Number of subdomains in x direction.',        type=int, default=variables.n_subdomains_x)
parser.add_argument('--n_subdomains_y', '-y',                help='Number of subdomains in y direction.',        type=int, default=variables.n_subdomains_y)
parser.add_argument('--n_subdomains_z', '-z',                help='Number of subdomains in z direction.',        type=int, default=variables.n_subdomains_z)
parser.add_argument('--potential_flow_solver_type',          help='The solver for the potential flow (non-spd matrix).',  default=variables.potential_flow_solver_type, choices=["gmres","cg","lu","gamg","richardson","chebyshev","cholesky","jacobi","sor","preonly"])
parser.add_argument('--potential_flow_preconditioner_type',  help='The preconditioner for the potential flow.',           default=variables.potential_flow_preconditioner_type, choices=["jacobi","sor","lu","ilu","gamg","none"])
parser.add_argument('--potential_flow_solver_maxit',         help='Maximum number of iterations for potential flow solver', type=int, default=variables.potential_flow_solver_maxit)
parser.add_argument('--potential_flow_solver_reltol',        help='Relative tolerance for potential flow solver',         type=float, default=variables.potential_flow_solver_reltol)
parser.add_argument('--emg_solver_type',                     help='The solver for the static bidomain.',                  default=variables.emg_solver_type)
#parser.add_argument('--emg_solver_type',                    help='The solver for the static bidomain.',                  default=variables.emg_solver_type, choices=["gmres","cg","lu","gamg","richardson","chebyshev","cholesky","jacobi","sor","preonly"])
parser.add_argument('--emg_preconditioner_type',             help='The preconditioner for the static bidomain.',          default=variables.emg_preconditioner_type, choices=["jacobi","sor","lu","ilu","gamg","none"])
parser.add_argument('--emg_solver_maxit',                    help='Maximum number of iterations for activation solver',   type=int, default=variables.emg_solver_maxit)
parser.add_argument('--emg_solver_reltol',                   help='Relative tolerance for activation solver',             type=float, default=variables.emg_solver_reltol)
parser.add_argument('--emg_solver_abstol',                   help='Absolute tolerance for activation solver',             type=float, default=variables.emg_solver_abstol)
parser.add_argument('--emg_initial_guess_nonzero',           help='If the initial guess for the emg linear system should be set to the previous solution.', default=variables.emg_initial_guess_nonzero, action='store_true')
parser.add_argument('--paraview_output',                     help='Enable the paraview output writer.',          default=variables.paraview_output, action='store_true')
parser.add_argument('--adios_output',                        help='Enable the MegaMol/ADIOS output writer.',          default=variables.adios_output, action='store_true')
parser.add_argument('--fiber_file',                          help='The filename of the file that contains the fiber data.', default=variables.fiber_file)
parser.add_argument('--firing_times_file',                   help='The filename of the file that contains the cellml model.', default=variables.firing_times_file)
parser.add_argument('--end_time', '--tend', '-t',            help='The end simulation time.',                    type=float, default=variables.end_time)
parser.add_argument('--approximate_exponential_function',    help='Approximate the exp function by a Taylor series',      default=variables.approximate_exponential_function, action="store_true")

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
    for j in [1]:
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
  print("dt_0D:           {:0.1e}, diffusion_solver_type:      {}".format(variables.dt_0D, variables.diffusion_solver_type))
  print("dt_1D:           {:0.1e}, potential_flow_solver_type: {}, approx. exp.: {}".format(variables.dt_1D, variables.potential_flow_solver_type, variables.approximate_exponential_function))
  print("dt_splitting:    {:0.1e}, emg_solver_type:            {}, emg_initial_guess_nonzero: {}".format(variables.dt_splitting, variables.emg_solver_type, variables.emg_initial_guess_nonzero))
  print("dt_3D:           {:0.1e}, paraview_output: {}, optimization_type: {}".format(variables.dt_3D, variables.paraview_output, variables.optimization_type))
  print("dt_elasticity:   {:0.0e}    elasticity solver: {}, preconditioner: {}".format(variables.dt_elasticity, variables.elasticity_solver_type, variables.elasticity_preconditioner_type))
  print("fiber_file:              {}".format(variables.fiber_file))
  print("fat_mesh_file:           {}".format(variables.fat_mesh_file))
  print("cellml_file:             {}".format(variables.cellml_file))
  print("fiber_distribution_file: {}".format(variables.fiber_distribution_file))
  print("firing_times_file:       {}".format(variables.firing_times_file))
  print("********************************************************************************")

  # start timer to measure duration of parsing of this script  
  t_start_script = timeit.default_timer()
    
# initialize all helper variables
from helper import *
  
# add neuron meshes
neuron_meshes = {
  "motoneuronMesh": {
    "nElements" :         variables.n_motoneurons-1 if n_ranks == 1 else variables.n_motoneurons*n_ranks,  # the last dof is empty in parallel
    "physicalExtent":     0,
    "physicalOffset":     0,
    "logKey":             "motoneuron",
    "inputMeshIsGlobal":  True,
    "nRanks":             n_ranks
  },
  "muscleSpindleMesh": {
    "nElements" :         variables.n_muscle_spindles-1 if n_ranks == 1 else variables.n_muscle_spindles*n_ranks,
    "physicalExtent":     0,
    "physicalOffset":     0,
    "logKey":             "muscle_spindle",
    "inputMeshIsGlobal":  True,
    "nRanks":             n_ranks
  }
}
variables.meshes.update(neuron_meshes)

print("stimulation_node_nos: {}".format(stimulation_node_nos))
print("muscle_spindle_node_nos: {}".format(muscle_spindle_node_nos))
print("updateSystemMatrixInterval: {}".format(int(variables.dt_elasticity/variables.dt_splitting)))
  
config = {
  "scenarioName":          variables.scenario_name,
  "logFormat":                      "csv",                                # "csv" or "json", format of the lines in the log file, csv gives smaller files
  "solverStructureDiagramFile":     "out/"+variables.scenario_name+"/solver_structure.txt",               # output file of a diagram that shows data connection between solvers
  "mappingsBetweenMeshesLogFile":   "out/"+variables.scenario_name+"/mappings_between_meshes_log.txt",    # log file for mappings 
  "meta": {                                                               # additional fields that will appear in the log
    "partitioning":         [variables.n_subdomains_x, variables.n_subdomains_y, variables.n_subdomains_z]
  },
  "Meshes":                variables.meshes,
  "MappingsBetweenMeshes": {
    "3Dmesh": [
       {"name": "3Dmesh_elasticity_quadratic+3DFatMesh_elasticity_quadratic", "xiTolerance": 0.35, "enableWarnings": True, "compositeUseOnlyInitializedMappings": True, "fixUnmappedDofs": True, "defaultValue": 0},
    ],
    "3DFatMesh":  [
       {"name": "3Dmesh_elasticity_quadratic+3DFatMesh_elasticity_quadratic", "xiTolerance": 0.2, "enableWarnings": True, "compositeUseOnlyInitializedMappings": True, "fixUnmappedDofs": True, "defaultValue": 0},
    ],
  },
  
  # connections of the slots, identified by slot name
  "connectedSlots": [
    ("mn",     "mn_out"),
    ("mn_out", "mn"),
    ("ms0",    "ms_in0"),
    ("ms1",    "ms_in1"),
    ("ms2",    "ms_in2"),
    ("ms3",    "ms_in3"),
    ("ms4",    "ms_in4"),
 #   ("stress", "m_g_in")     # connection of activation from subcellular model to muscle contraction
  ],
  
  "Solvers": {
    "precontractionMechanicsSolver": {   # solver for the preprocessing contraction simulation (static mechanics problem)
      "relativeTolerance":   1e-5,           # 1e-10 relative tolerance of the linear solver
      "absoluteTolerance":   1e-10,           # 1e-10 absolute tolerance of the residual of the linear solver
      "solverType":          "lu",            # type of the linear solver
      "preconditionerType":  "none",          # type of the preconditioner
      "maxIterations":       1e4,                                         # maximum number of iterations in the linear solver
      "snesMaxFunctionEvaluations": 1e8,                                  # maximum number of function iterations
      "snesMaxIterations":   34,              # maximum number of iterations in the nonlinear solver
      "snesRelativeTolerance": 1e-2,         # relative tolerance of the nonlinear solver
      "snesAbsoluteTolerance": 1e-2,         # absolute tolerance of the nonlinear solver
      "snesLineSearchType": "l2",                                         # type of linesearch, possible values: "bt" "nleqerr" "basic" "l2" "cp" "ncglinear"
      "snesRebuildJacobianFrequency": 1,    # how often the jacobian should be recomputed, -1 indicates NEVER rebuild, 1 means rebuild every time the Jacobian is computed within a single nonlinear solve, 2 means every second time the Jacobian is built etc. -2 means rebuild at next chance but then never again 
      "hypreOptions":        "",                                          # additional options for the hypre solvers could be given here
      "dumpFilename":        "",                                          # dump system matrix and right hand side after every solve
      "dumpFormat":          "matlab",                                    # default, ascii, matlab
    },
    "prestretchMechanicsSolver": {   # solver for the prestretch simulation (static mechanics problem)
      "relativeTolerance":   1e-5,          # 1e-10 relative tolerance of the linear solver
      "absoluteTolerance":   1e-10,         # 1e-10 absolute tolerance of the residual of the linear solver
      "solverType":          "lu",          # type of the linear solver
      "preconditionerType":  "none",        # type of the preconditioner
      "maxIterations":       1e4,           # maximum number of iterations in the linear solver
      "snesMaxFunctionEvaluations": 1e8,    # maximum number of function iterations
      "snesMaxIterations":   34,            # maximum number of iterations in the nonlinear solver
      "snesRelativeTolerance": 1e-10,       # relative tolerance of the nonlinear solver
      "snesAbsoluteTolerance": 1e-10,       # absolute tolerance of the nonlinear solver
      "snesLineSearchType": "l2",           # type of linesearch, possible values: "bt" "nleqerr" "basic" "l2" "cp" "ncglinear"
      "snesRebuildJacobianFrequency": 1,    # how often the jacobian should be recomputed, -1 indicates NEVER rebuild, 1 means rebuild every time the Jacobian is computed within a single nonlinear solve, 2 means every second time the Jacobian is built etc. -2 means rebuild at next chance but then never again 
      "hypreOptions":        "",            # additional options for the hypre solvers could be given here
      "dumpFilename":        "",            # dump system matrix and right hand side after every solve
      "dumpFormat":          "matlab",      # default, ascii, matlab
    },
    "implicitSolver": {     # solver for the implicit timestepping scheme of the diffusion time step
      "relativeTolerance":  variables.diffusion_solver_reltol,
      "absoluteTolerance":  1e-10,         # 1e-10 absolute tolerance of the residual    
      "maxIterations":      variables.diffusion_solver_maxit,
      "solverType":         variables.diffusion_solver_type,
      "preconditionerType": variables.diffusion_preconditioner_type,
      "dumpFilename":       "",   # "out/dump_"
      "dumpFormat":         "matlab",
    },
    "potentialFlowSolver": {
      "relativeTolerance":  1e-10,
      "absoluteTolerance":  1e-10,         # 1e-10 absolute tolerance of the residual          
      "maxIterations":      1e5,
      "solverType":         "gmres",
      "preconditionerType": "none",
      "dumpFormat":         "default",
      "dumpFilename":       "",
    },
    "activationSolver": {   # solver for the static Bidomain equation and the EMG
      "relativeTolerance":  variables.emg_solver_reltol,
      "absoluteTolerance":  variables.emg_solver_abstol, 
      "maxIterations":      variables.emg_solver_maxit,
      "solverType":         variables.emg_solver_type,
      "preconditionerType": variables.emg_preconditioner_type,
      "dumpFilename":       "",
      "dumpFormat":         "ascii",
    },
    "mechanicsSolver": {   # solver for the dynamic mechanics problem
      "relativeTolerance":   variables.linear_relative_tolerance,           # 1e-10 relative tolerance of the linear solver
      "absoluteTolerance":   variables.linear_absolute_tolerance,           # 1e-10 absolute tolerance of the residual of the linear solver
      "solverType":          variables.elasticity_solver_type,            # type of the linear solver
      "preconditionerType":  variables.elasticity_preconditioner_type,    # type of the preconditioner
      "maxIterations":       1e4,                                         # maximum number of iterations in the linear solver
      "snesMaxFunctionEvaluations": 1e8,                                  # maximum number of function iterations
      "snesMaxIterations":   variables.snes_max_iterations,               # maximum number of iterations in the nonlinear solver
      "snesRelativeTolerance": variables.snes_relative_tolerance,         # relative tolerance of the nonlinear solver
      "snesAbsoluteTolerance": variables.snes_absolute_tolerance,         # absolute tolerance of the nonlinear solver
      "snesLineSearchType": "l2",                                         # type of linesearch, possible values: "bt" "nleqerr" "basic" "l2" "cp" "ncglinear"
      "snesRebuildJacobianFrequency": variables.snes_rebuild_jacobian_frequency,    # how often the jacobian should be recomputed, -1 indicates NEVER rebuild, 1 means rebuild every time the Jacobian is computed within a single nonlinear solve, 2 means every second time the Jacobian is built etc. -2 means rebuild at next chance but then never again 
      "hypreOptions":        "",                                          # additional options for the hypre solvers could be given here
      "dumpFilename":        "",                                          # dump system matrix and right hand side after every solve
      "dumpFormat":          "matlab",                                    # default, ascii, matlab
    }
  },
  
  # preprocessing (prescribe activation, precontraction, prestretch) and main solver
  # "MultipleCoupling": {
  #   "timeStepWidth":          variables.end_time,
  #   "logTimeStepWidthAsKey":  "dt_artifical_contraction",
  #   "durationLogKey":         "duration_artifical_contraction",
  #   "timeStepOutputInterval": 1,
  #   "endTime":                variables.end_time,
  #   "deferInitialization":    True,       # initialize nested solvers only right before computation of first timestep
  #   "connectedSlotsTerm1To2": None,       # connect lambda to slot 0 and gamma to slot 2
  #   "connectedSlotsTerm2To1": None,       # transfer nothing back
  # 
  #   # prescribed activation
  #   "Term1": {
  #     "PrescribedValues": {
  #       "meshName":               ["3Dmesh", "3DFatMesh"],      # reference to the fibers mesh
  #       "numberTimeSteps":        1,             # number of timesteps to call the callback functions subsequently, this is usually 1 for prescribed values, because it is enough to set the reaction term only once per time step
  #       "timeStepOutputInterval": 20,            # if the time step should be written to console, a value > 10 produces no output
  #       "slotNames":              ["lambda", "gamma"],
  # 
  #       # a list of field variables that will get values assigned in every timestep, by the provided callback function
  #       "fieldVariables1": [
  #         {"name": "lambda", "callback": variables.set_lambda_values},
  #         {"name": "gamma", "callback": variables.set_gamma_values},
  #       ],
  #       "fieldVariables2":     [],
  #       "additionalArgument":  None,         # a custom argument to the fieldVariables callback functions, this will be passed on as the last argument
  # 
  #       "OutputWriter" : [
  #         {"format": "Paraview", "outputInterval": 1, "filename": "out/" + variables.scenario_name + "/0_prescribed_stress", "binary": True, "fixedFormat": False, "combineFiles": True, "fileNumbering": "incremental"}
  #       ]
  #     }
  #   },
  #   # artifical contraction with constant gamma
  #   "Term2": {  
  #     "MuscleContractionSolver": {
  #       "numberTimeSteps":              1,                         # only use 1 timestep per interval
  #       "timeStepOutputInterval":       100,                       # do not output time steps
  #       "Pmax":                         variables.pmax,            # maximum PK2 active stress
  #       "enableForceLengthRelation":    False,                     # if the factor f_l(λ_f) modeling the force-length relation (as in Heidlauf2013) should be multiplied. Set to false if this relation is already considered in the CellML model.
  #       "lambdaDotScalingFactor":       1.0,                       # scaling factor for the output of the lambda dot slot, i.e. the contraction velocity. Use this to scale the unit-less quantity to, e.g., micrometers per millisecond for the subcellular model.
  #       "slotNames":                    ["lambda", "ldot", "gamma", "T"],
  #       "OutputWriter" : [
  #         {"format": "Paraview", "outputInterval": int(1./variables.dt_elasticity*variables.output_timestep_elasticity), "filename": "out/"+variables.scenario_name+"/1_precontraction_mechanics_3D", "binary": True, "fixedFormat": False, "onlyNodalValues":True, "combineFiles": True, "fileNumbering": "incremental"},
  #       ],
  #       "mapGeometryToMeshes":          ["3Dmesh_elasticity_quadratic+3DFatMesh_elasticity_quadratic"],    # the mesh names of the meshes that will get the geometry transferred
  #       "reverseMappingOrder":          True,                      # if the mapping target->own mesh should be used instead of own->target mesh. This gives better results in some cases.
  #       "dynamic":                      False,                     # if the dynamic solid mechanics solver should be used, else it computes the quasi-static problem
  # 
  #       # the actual solid mechanics solver, this is either "DynamicHyperelasticitySolver" or "HyperelasticitySolver", depending on the value of "dynamic"
  # 
  #       "HyperelasticitySolver": {
  #         "durationLogKey":             "duration_precontraction",               # key to find duration of this solver in the log file
  # 
  #         "materialParameters":         variables.material_parameters,  # material parameters of the Mooney-Rivlin material
  #         "displacementsScalingFactor": 1.0,                       # scaling factor for displacements, only set to sth. other than 1 only to increase visual appearance for very small displacements
  #         "residualNormLogFilename":    "out/"+variables.scenario_name+"/1_log_precontraction_residual_norm.txt",   # log file where residual norm values of the nonlinear solver will be written
  #         "useAnalyticJacobian":        True,                      # whether to use the analytically computed jacobian matrix in the nonlinear solver (fast)
  #         "useNumericJacobian":         False,                     # whether to use the numerically computed jacobian matrix in the nonlinear solver (slow), only works with non-nested matrices, if both numeric and analytic are enable, it uses the analytic for the preconditioner and the numeric as normal jacobian
  # 
  #         "dumpDenseMatlabVariables":   False,                     # whether to have extra output of matlab vectors, x,r, jacobian matrix (very slow)
  #         # if useAnalyticJacobian,useNumericJacobian and dumpDenseMatlabVariables all all three true, the analytic and numeric jacobian matrices will get compared to see if there are programming errors for the analytic jacobian
  # 
  #         # mesh
  #         "inputMeshIsGlobal":          True,                     # boundary conditions and initial values are given as global numbers (every process has all information)
  #         "meshName":                   ["3Dmesh_elasticity_quadratic_precontraction", "3DFatMesh_elasticity_quadratic_precontraction"],       # name of the 3D mesh, it is defined under "Meshes" at the beginning of this config
  #         "fiberMeshNames":             [],                       # fiber meshes that will be used to determine the fiber direction, there are no fibers in multidomain, so this is empty
  #         "fiberDirection":             [0,0,1],                  # if fiberMeshNames is empty, directly set the constant fiber direction, in element coordinate system
  # 
  #         # solving
  #         "solverName":                 "precontractionMechanicsSolver",         # name of the nonlinear solver configuration, it is defined under "Solvers" at the beginning of this config
  #         #"loadFactors":                [0.5, 1.0],                # load factors for every timestep
  #         "loadFactors":                list(np.logspace(-1,0,3)),        # load factors
  #         #"loadFactors":                [],                        # no load factors
  #         "loadFactorGiveUpThreshold":  1e-5,                      # a threshold for the load factor, when to abort the solve of the current time step. The load factors are adjusted automatically if the nonlinear solver diverged. If the load factors get too small, it aborts the solve.
  #         "scaleInitialGuess":          False,                      # when load stepping is used, scale initial guess between load steps a and b by sqrt(a*b)/a. This potentially reduces the number of iterations per load step (but not always).
  #         "nNonlinearSolveCalls":       1,                         # how often the nonlinear solve should be repeated
  # 
  #         # boundary and initial conditions
  #         "dirichletBoundaryConditions": variables.precontraction_elasticity_dirichlet_bc,   # the initial Dirichlet boundary conditions that define values for displacements u and velocity v
  #         "neumannBoundaryConditions":   variables.precontraction_elasticity_neumann_bc,     # Neumann boundary conditions that define traction forces on surfaces of elements
  #         "divideNeumannBoundaryConditionValuesByTotalArea": True,            # if the given Neumann boundary condition values under "neumannBoundaryConditions" are total forces instead of surface loads and therefore should be scaled by the surface area of all elements where Neumann BC are applied
  #         "updateDirichletBoundaryConditionsFunction": None,                  # function that updates the dirichlet BCs while the simulation is running
  #         "updateDirichletBoundaryConditionsFunctionCallInterval": 1,         # every which step the update function should be called, 1 means every time step
  # 
  #         "initialValuesDisplacements":  [[0.0,0.0,0.0] for _ in range(variables.n_points_global_composite_mesh)],     # the initial values for the displacements, vector of values for every node [[node1-x,y,z], [node2-x,y,z], ...]
  #         "initialValuesVelocities":     [[0.0,0.0,0.0] for _ in range(variables.n_points_global_composite_mesh)],     # the initial values for the velocities, vector of values for every node [[node1-x,y,z], [node2-x,y,z], ...]
  #         "extrapolateInitialGuess":     True,                                # if the initial values for the dynamic nonlinear problem should be computed by extrapolating the previous displacements and velocities
  #         "constantBodyForce":           variables.precontraction_constant_body_force,       # a constant force that acts on the whole body, e.g. for gravity
  # 
  #         "dirichletOutputFilename":    "out/"+variables.scenario_name+"/1_precontraction_dirichlet_boundary_conditions",     # output filename for the dirichlet boundary conditions, set to "" to have no output
  # 
  #         # define which file formats should be written
  #         # 1. main output writer that writes output files using the quadratic elements function space. Writes displacements, velocities and PK2 stresses.
  #         "OutputWriter" : [
  # 
  #           # Paraview files
  #           {"format": "Paraview", "outputInterval": 1, "filename": "out/"+variables.scenario_name+"/1_precontraction_u", "binary": True, "fixedFormat": False, "onlyNodalValues":True, "combineFiles":True, "fileNumbering": "incremental"},
  #         ],
  #         # 2. additional output writer that writes also the hydrostatic pressure
  #         "pressure": {   # output files for pressure function space (linear elements), contains pressure values, as well as displacements and velocities
  #           "OutputWriter" : [
  #             {"format": "Paraview", "outputInterval": 1, "filename": "out/"+variables.scenario_name+"/1_precontraction_p", "binary": True, "fixedFormat": False, "onlyNodalValues":True, "combineFiles":True, "fileNumbering": "incremental"},
  #           ]
  #         },
  #         # 4. output writer for debugging, outputs files after each load increment, the geometry is not changed but u and v are written
  #         "LoadIncrements": {   
  #           "OutputWriter" : [
  #             {"format": "Paraview", "outputInterval": 1, "filename": "out/"+variables.scenario_name+"/1_precontraction_load_increments", "binary": True, "fixedFormat": False, "onlyNodalValues":True, "combineFiles":True, "fileNumbering": "incremental"},
  #           ]
  #         }
  #       }
  #     }
  #   },
  #   # prestretch simulation
  #   "Term3": {
  #     "HyperelasticitySolver": {
  #       "durationLogKey":             "duration_prestretch",               # key to find duration of this solver in the log file
  # 
  #       "materialParameters":         variables.material_parameters,  # material parameters of the Mooney-Rivlin material
  #       "displacementsScalingFactor": 1.0,                       # scaling factor for displacements, only set to sth. other than 1 only to increase visual appearance for very small displacements
  #       "residualNormLogFilename":    "out/"+variables.scenario_name+"/2_log_prestretch_residual_norm.txt",   # log file where residual norm values of the nonlinear solver will be written
  #       "useAnalyticJacobian":        True,                      # whether to use the analytically computed jacobian matrix in the nonlinear solver (fast)
  #       "useNumericJacobian":         False,                     # whether to use the numerically computed jacobian matrix in the nonlinear solver (slow), only works with non-nested matrices, if both numeric and analytic are enable, it uses the analytic for the preconditioner and the numeric as normal jacobian
  #       "slotNames":                  ["m_ux", "m_uy", "m_uz"],  # slot names of the data connector slots, there are three slots, namely the displacement components ux, uy, uz
  # 
  #       "dumpDenseMatlabVariables":   False,                     # whether to have extra output of matlab vectors, x,r, jacobian matrix (very slow)
  #       # if useAnalyticJacobian,useNumericJacobian and dumpDenseMatlabVariables all all three true, the analytic and numeric jacobian matrices will get compared to see if there are programming errors for the analytic jacobian
  # 
  #       # mesh
  #       "inputMeshIsGlobal":          True,                     # boundary conditions and initial values are given as global numbers (every process has all information)
  #       "meshName":                   ["3Dmesh_elasticity_quadratic", "3DFatMesh_elasticity_quadratic"],       # name of the 3D mesh, it is defined under "Meshes" at the beginning of this config
  #       "fiberMeshNames":             [],                       # fiber meshes that will be used to determine the fiber direction, there are no fibers in multidomain, so this is empty
  #       "fiberDirection":             [0,0,1],                  # if fiberMeshNames is empty, directly set the constant fiber direction, in element coordinate system
  # 
  #       # solving
  #       "solverName":                 "prestretchMechanicsSolver",         # name of the nonlinear solver configuration, it is defined under "Solvers" at the beginning of this config
  #       #"loadFactors":                [0.5, 1.0],                # load factors for every timestep
  #       "loadFactors":                list(np.logspace(-1,0,3)),        # load factors
  #       #"loadFactors":                [],                        # no load factors
  #       "loadFactorGiveUpThreshold":  1e-5,                      # a threshold for the load factor, when to abort the solve of the current time step. The load factors are adjusted automatically if the nonlinear solver diverged. If the load factors get too small, it aborts the solve.
  #       "scaleInitialGuess":          False,                      # when load stepping is used, scale initial guess between load steps a and b by sqrt(a*b)/a. This potentially reduces the number of iterations per load step (but not always).
  #       "nNonlinearSolveCalls":       1,                         # how often the nonlinear solve should be repeated
  # 
  #       # boundary and initial conditions
  #       "dirichletBoundaryConditions": variables.prestretch_elasticity_dirichlet_bc,   # the initial Dirichlet boundary conditions that define values for displacements u and velocity v
  #       "neumannBoundaryConditions":   variables.prestretch_elasticity_neumann_bc,     # Neumann boundary conditions that define traction forces on surfaces of elements
  #       "divideNeumannBoundaryConditionValuesByTotalArea": True,            # if the given Neumann boundary condition values under "neumannBoundaryConditions" are total forces instead of surface loads and therefore should be scaled by the surface area of all elements where Neumann BC are applied
  #       "updateDirichletBoundaryConditionsFunction": None,                  # function that updates the dirichlet BCs while the simulation is running
  #       "updateDirichletBoundaryConditionsFunctionCallInterval": 1,         # every which step the update function should be called, 1 means every time step
  # 
  #       "initialValuesDisplacements":  [[0.0,0.0,0.0] for _ in range(variables.n_points_global_composite_mesh)],     # the initial values for the displacements, vector of values for every node [[node1-x,y,z], [node2-x,y,z], ...]
  #       "initialValuesVelocities":     [[0.0,0.0,0.0] for _ in range(variables.n_points_global_composite_mesh)],     # the initial values for the velocities, vector of values for every node [[node1-x,y,z], [node2-x,y,z], ...]
  #       "extrapolateInitialGuess":     True,                                # if the initial values for the dynamic nonlinear problem should be computed by extrapolating the previous displacements and velocities
  #       "constantBodyForce":           variables.prestretch_constant_body_force,       # a constant force that acts on the whole body, e.g. for gravity
  # 
  #       "dirichletOutputFilename":    "out/"+variables.scenario_name+"/2_prestretch_dirichlet_boundary_conditions",     # output filename for the dirichlet boundary conditions, set to "" to have no output
  # 
  #       # define which file formats should be written
  #       # 1. main output writer that writes output files using the quadratic elements function space. Writes displacements, velocities and PK2 stresses.
  #       "OutputWriter" : [
  # 
  #         # Paraview files
  #         {"format": "Paraview", "outputInterval": 1, "filename": "out/"+variables.scenario_name+"/2_prestretch_u", "binary": True, "fixedFormat": False, "onlyNodalValues":True, "combineFiles":True, "fileNumbering": "incremental"},
  #       ],
  #       # 2. additional output writer that writes also the hydrostatic pressure
  #       "pressure": {   # output files for pressure function space (linear elements), contains pressure values, as well as displacements and velocities
  #         "OutputWriter" : [
  #           {"format": "Paraview", "outputInterval": 1, "filename": "out/"+variables.scenario_name+"/2_prestretch_p", "binary": True, "fixedFormat": False, "onlyNodalValues":True, "combineFiles":True, "fileNumbering": "incremental"},
  #         ]
  #       },
  #       # 4. output writer for debugging, outputs files after each load increment, the geometry is not changed but u and v are written
  #       "LoadIncrements": {   
  #         "OutputWriter" : [
  #           {"format": "Paraview", "outputInterval": 1, "filename": "out/"+variables.scenario_name+"/2_prestretch_load_increments", "binary": True, "fixedFormat": False, "onlyNodalValues":True, "combineFiles":True, "fileNumbering": "incremental"},
  #         ]
  #       },
  #     }
  #   },
  # 
  #   # neuron solvers and fibers simulation
  #   "Term4": {
  # 
  #     # coupling scheme of neuron solvers
  #     "MultipleCoupling": {
  #       "endTime":                variables.end_time,                         # end time of the simulation
  #       "timeStepWidth":          variables.dt_neuron_transfer,               # time step width of the data transfer between the sub solvers
  #       "logTimeStepWidthAsKey":  "dt_neuron_transfer",                       # string under which the timestep width will be stored in the log file
  #       "durationLogKey":         "duration_total",                           # string under which the total duration will be stored in the log file
  #       "timeStepOutputInterval": 500,                                        # how often to print the current time step
  #       "connectedSlotsTerm1To2": None,                                       # this would be the connected slots for a normal Coupling scheme, but here it should be set to None, use global option "connectedSlots" instead
  #       "connectedSlotsTerm2To1": None,                                       # this would be the connected slots for a normal Coupling scheme, but here it should be set to None, use global option "connectedSlots" instead
  # 
  #       # muscle spindles
  #       "Term1": {
  #         # mapping muscle spindles output -> motor neuron input
  #         "MapDofs": {
  #           "description":                "muscle_spindles_to_motoneurons",   # description that will be shown in solver structure visualization
  #           "nAdditionalFieldVariables":  1,                                  # number of additional field variables that are defined by this object. They have 1 component, use the templated function space and mesh given by meshName.
  #           "additionalSlotNames":        ["mn_in"],
  #           "meshName":                   "motoneuronMesh",                   # the mesh on which the additional field variables will be defined
  #           "beforeComputation": None,
  #           "afterComputation": [                                         # transfer/mapping of dofs that will be performed after the computation of the nested solver
  #             {                                                 
  #               "fromConnectorSlot":                "ms_out",
  #               "toConnectorSlots":                 ["mn_in"],
  #               "fromSlotConnectorArrayIndex":      0,                    # which fiber/compartment
  #               "toSlotConnectorArrayIndex":        0,
  #               "mode":                             "callback",           # "copyLocal", "copyLocalIfPositive", "localSetIfAboveThreshold" or "communicate"
  #               "fromDofNosNumbering":              "local",
  #               "toDofNosNumbering":                "local",
  #               "dofsMapping":                      None,
  #               "inputDofs":                        list(range(variables.n_muscle_spindles)),   # [0,1,...,n_muscle_spindles], this is for mesh "muscleSpindleMesh"
  #               "outputDofs":                       [list(range(variables.n_motoneurons))],   # [0,1,...,n_motor_neurons], this is for mesh "motoneuronMesh"
  #               "callback":                         variables.callback_muscle_spindles_to_motoneurons,
  #               #"thresholdValue":                   20,                  # if mode is "localSetIfAboveThreshold", this is the threshold, if the value is above it, set the value `valueToSet`
  #               #"valueToSet":                       20,                  # if mode is "localSetIfAboveThreshold", this is the value to set the target dof to, if the source dof is above thresholdValue.
  #             }
  #           ],
  # 
  #           # muscle spindle model solver
  #           "Heun" : {
  #             "description":                  "muscle spindle",
  #             "timeStepWidth":                variables.dt_muscle_spindles,
  #             "logTimeStepWidthAsKey":        "dt_muscle_spindles",
  #             "durationLogKey":               "duration_muscle_spindles",
  #             "initialValues":                [],
  #             "timeStepOutputInterval":       1e4,
  #             "inputMeshIsGlobal":            True,
  #             "dirichletBoundaryConditions":  {},
  #             "dirichletOutputFilename":      None,                                 # filename for a vtp file that contains the Dirichlet boundary condition nodes and their values, set to None to disable
  #             "checkForNanInf":               True,             # check if the solution vector contains nan or +/-inf values, if yes, an error is printed. This is a time-consuming check.
  #             "nAdditionalFieldVariables":    0,
  #             "additionalSlotNames":          [],
  # 
  #             # cellml model of muscle spindle
  #             "CellML" : {
  #               "modelFilename":                          variables.muscle_spindle_cellml_file,           # input C++ source file or cellml XML file
  #               "initializeStatesToEquilibrium":          False,                                          # if the equilibrium values of the states should be computed before the simulation starts
  #               "initializeStatesToEquilibriumTimestepWidth": 1e-4,                                       # if initializeStatesToEquilibrium is enable, the timestep width to use to solve the equilibrium equation
  # 
  #               # optimization parameters
  #               "optimizationType":                       "vc",                                           # "vc", "simd", "openmp" type of generated optimizated source file
  #               "approximateExponentialFunction":         True,                                           # if optimizationType is "vc", whether the exponential function exp(x) should be approximate by (1+x/n)^n with n=1024
  #               "compilerFlags":                          "-fPIC -O3 -march=native -shared ",             # compiler flags used to compile the optimized model code
  #               "maximumNumberOfThreads":                 0,                                              # if optimizationType is "openmp", the maximum number of threads to use. Default value 0 means no restriction.
  # 
  #               # stimulation callbacks, motor neuron is not stimulated by a callback function, but has a constant stimulation current
  #               "setSpecificStatesFunction":              None,                                           # callback function that sets states like Vm, activation can be implemented by using this method and directly setting Vm values, or by using setParameters/setSpecificParameters
  #               #"setSpecificStatesCallInterval":         0,                                              # set_specific_states should be called variables.stimulation_frequency times per ms, the factor 2 is needed because every Heun step includes two calls to rhs
  #               "setSpecificStatesCallInterval":          0,                                              # 0 means disabled
  #               "setSpecificStatesCallFrequency":         0,                                              # set_specific_states should be called variables.stimulation_frequency times per ms
  #               "setSpecificStatesFrequencyJitter":       0,                                              # random value to add or substract to setSpecificStatesCallFrequency every stimulation, this is to add random jitter to the frequency
  #               "setSpecificStatesRepeatAfterFirstCall":  0.01,                                           # [ms] simulation time span for which the setSpecificStates callback will be called after a call was triggered
  #               "setSpecificStatesCallEnableBegin":       0,                                              # [ms] first time when to call setSpecificStates
  #               "additionalArgument":                     None,
  # 
  #               "handleResultFunction":                   None, #handle_result,              # callback function that gets all current values and can do something with them
  #               "handleResultCallInterval":               1,                          # interval in which handle_result will be called
  #               "handleResultFunctionAdditionalParameter": None,                      # additional last argument for handle_result
  # 
  #               "mappings":                               variables.muscle_spindle_mappings,              # mappings between parameters and algebraics/constants and between outputConnectorSlots and states, algebraics or parameters, they are defined in helper.py
  #               "parametersInitialValues":                variables.muscle_spindle_parameters_initial_values,  # initial values for the parameters, either once for all instances, or for all instances in array of struct ordering with nParameters_ parameters per instance: [inst0p0, inst0p1, ... inst0pn, inst1p0, inst1p1, ...]
  # 
  #               "meshName":                               "muscleSpindleMesh",                                       
  #               "stimulationLogFilename":                 "out/stimulation.log",
  # 
  #               # output writer for states, algebraics and parameters                
  #               "OutputWriter" : [
  #                 {"format": "PythonFile", "outputInterval": (int)(2./variables.dt_muscle_spindles*variables.output_timestep_neurons), "filename": "out/muscle_spindles", "binary": True, "fixedFormat": False, "combineFiles": True, "onlyNodalValues": True, "fileNumbering": "incremental"}
  #               ]
  #             }
  #           }
  #         }
  #       },
  # 
  #       "Term2": {
  # 
                              "MapDofs": {
                                "description":                "motoneurons->stimulated nodes",  # description that will be shown in solver structure visualization
                                "nAdditionalFieldVariables":  0,                              # number of additional field variables that are defined by this object. They have 1 component, use the templated function space and mesh given by meshName.
                                "additionalSlotNames":        [],
                                "meshName":                   "motoneuronMesh",               # the mesh on which the additional field variables will be defined
                                "beforeComputation":          None,

                                # mapping from motoneuronMesh which contains on every rank as many nodes as there are motoneurons to the 3D domain
                                # map from motoneuronMesh (algebraics) to the fiber meshes (solution)
                                "afterComputation": [                                        # transfer/mapping of dofs that will be performed before the computation of the nested solver
                                  {
                                    "fromConnectorSlot":                "mn_out",                # source slot of the dofs mapping
                                    "toConnectorSlots":                 "mn_out",                # target slot of the dofs mapping
                                    "fromSlotConnectorArrayIndex":      0,
                                    "toSlotConnectorArrayIndex":        0,      # which fiber in this motor unit
                                    "mode":                             "localSetIfAboveThreshold",          # "copyLocal", "copyLocalIfPositive", "localSetIfAboveThreshold" or "communicate"
                                    "fromDofNosNumbering":              "local",
                                    "toDofNosNumbering":                "global",
                                    "dofsMapping":                      {i:i for i in range(variables.n_motoneurons)},   # map from the motor unit to the stimulated node of the fiber mesh
                                    "inputDofs":                        None,                # this option is only needed in mode "callback"
                                    "outputDofs":                       None,                # this option is only needed in mode "callback"
                                    "callback":                         None,                # this option is only needed in mode "callback"
                                    "thresholdValue":                   100000,                  # if mode is "localSetIfAboveThreshold", this is the threshold, if the value is above it, set the value `valueToSet`
                                    "valueToSet":                       42,       # if mode is "localSetIfAboveThreshold", this is the value to set the target dof to, if the source dof is above thresholdValue.
                        }],
                     
          # motoneuron solver
          "Heun" : {
            "endTime": 2e-3,
            "description":                  "motoneurons",
            "timeStepWidth":                1e-3,#variables.dt_motoneuron,
            "logTimeStepWidthAsKey":        "dt_motoneuron",
            "durationLogKey":               "duration_motoneuron",
            "initialValues":                [],
            "timeStepOutputInterval":       1,
            "inputMeshIsGlobal":            True,
            "dirichletBoundaryConditions":  {},
            "dirichletOutputFilename":      None,                                 # filename for a vtp file that contains the Dirichlet boundary condition nodes and their values, set to None to disable
            "checkForNanInf":               True,             # check if the solution vector contains nan or +/-inf values, if yes, an error is printed. This is a time-consuming check.
            "nAdditionalFieldVariables":    0,
            "additionalSlotNames":          [],
                
            # cellml model of motorneuron
            "CellML" : {
              "modelFilename":                          variables.motoneuron_cellml_file,               # input C++ source file or cellml XML file
              "initializeStatesToEquilibrium":          False,                                          # if the equilibrium values of the states should be computed before the simulation starts
              "initializeStatesToEquilibriumTimestepWidth": 1e-4,                                       # if initializeStatesToEquilibrium is enable, the timestep width to use to solve the equilibrium equation
              
              # optimization parameters
              "optimizationType":                       "vc",                                           # "vc", "simd", "openmp" type of generated optimizated source file
              "approximateExponentialFunction":         True,                                           # if optimizationType is "vc", whether the exponential function exp(x) should be approximate by (1+x/n)^n with n=1024
              "compilerFlags":                          "-fPIC -O3 -march=native -shared ",             # compiler flags used to compile the optimized model code
              "maximumNumberOfThreads":                 0,                                              # if optimizationType is "openmp", the maximum number of threads to use. Default value 0 means no restriction.
              
              # stimulation callbacks, motor neuron is not stimulated by a callback function, but has a constant stimulation current
              "setSpecificStatesFunction":              None,                                           # callback function that sets states like Vm, activation can be implemented by using this method and directly setting Vm values, or by using setParameters/setSpecificParameters
              #"setSpecificStatesCallInterval":         0,                                              # set_specific_states should be called variables.stimulation_frequency times per ms, the factor 2 is needed because every Heun step includes two calls to rhs
              "setSpecificStatesCallInterval":          0,                                              # 0 means disabled
              "setSpecificStatesCallFrequency":         0,                                              # set_specific_states should be called variables.stimulation_frequency times per ms
              "setSpecificStatesFrequencyJitter":       0,                                              # random value to add or substract to setSpecificStatesCallFrequency every stimulation, this is to add random jitter to the frequency
              "setSpecificStatesRepeatAfterFirstCall":  0.01,                                           # [ms] simulation time span for which the setSpecificStates callback will be called after a call was triggered
              "setSpecificStatesCallEnableBegin":       0,                                              # [ms] first time when to call setSpecificStates
              "additionalArgument":                     None,
              
              "mappings":                               variables.motoneuron_mappings,                  # mappings between parameters and algebraics/constants and between outputConnectorSlots and states, algebraics or parameters, they are defined in helper.py
              "parametersInitialValues":                variables.motoneuron_parameters_initial_values, # initial values for the parameters, either once for all instances, or for all instances in array of struct ordering with nParameters_ parameters per instance: [inst0p0, inst0p1, ... inst0pn, inst1p0, inst1p1, ...]
              
              "meshName":                               "motoneuronMesh",                               # use the linear mesh, it was partitioned by the helper.py script which called opendihu/scripts/create_partitioned_meshes_for_settings.py
              "stimulationLogFilename":                 "out/stimulation.log",

              # output writer for states, algebraics and parameters                
              "OutputWriter" : [
                {"format": "Paraview",   "outputInterval": 1+0*int(2./variables.dt_motoneuron*variables.output_timestep_motoneuron), "filename": "out/motoneurons", "binary": True, "fixedFormat": False, "combineFiles": True, "onlyNodalValues": True, "fileNumbering": "incremental"},
                # {"format": "PythonFile", "outputInterval": (int)(2./variables.dt_motoneuron*variables.output_timestep_motoneuron), "filename": "out/motoneurons", "binary": True, "fixedFormat": False, "combineFiles": True, "onlyNodalValues": True, "fileNumbering": "incremental"}
              ]
            }
          }
          }
  #       },
  # 
  #       # electrophysiology and mechanics solver
  #       "Term3": {
  #         # map from λ in the 3D mesh to muscle spindles input
  #         "MapDofs": {
  #           "description":                "muscle_spindles_input",        # description that will be shown in solver structure visualization
  #           "nAdditionalFieldVariables":  5,                              # number of additional field variables that are defined by this object. They have 1 component, use the templated function space and mesh given by meshName.
  #           "additionalSlotNames":        ["ms0","ms1","ms2","ms3","ms4"],
  #           "meshName":                   "muscleSpindleMesh",            # the mesh on which the additional field variables will be defined
  #           "beforeComputation": [       # transfer/mapping of dofs that will be performed before the computation of the nested solver
  #             # read spindle stretch (slot m_lda) and communicate to all processes
  #             {                                                 
  #               "fromConnectorSlot":                "m_lda",
  #               "toConnectorSlots":                 "ms0",
  #               "fromSlotConnectorArrayIndex":      0,                    # which fiber/compartment
  #               "toSlotConnectorArrayIndex":        0,
  #               "mode":                             "communicate",        # "copyLocal", "copyLocalIfPositive", "localSetIfAboveThreshold" or "communicate"
  #               "fromDofNosNumbering":              "global",
  #               "toDofNosNumbering":                "global",
  #               "dofsMapping":                      
  #                 {muscle_spindle_dof : [rank_no*variables.n_muscle_spindles + i for rank_no in range(n_ranks)]
  #                  for i,muscle_spindle_dof in enumerate(muscle_spindle_node_nos)},
  #               "inputDofs":                        None,                
  #               "outputDofs":                       None,                
  #               "callback":                         None,
  #               #"thresholdValue":                   20,                  # if mode is "localSetIfAboveThreshold", this is the threshold, if the value is above it, set the value `valueToSet`
  #               #"valueToSet":                       20,                  # if mode is "localSetIfAboveThreshold", this is the value to set the target dof to, if the source dof is above thresholdValue.
  #             }, 
  #             # call callback_muscle_spindles_input 
  #             {                                                 
  #               "fromConnectorSlot":                "ms0",
  #               "toConnectorSlots":                 ["ms0","ms1","ms2","ms3","ms4"],
  #               "fromSlotConnectorArrayIndex":      0,                    # which fiber/compartment
  #               "toSlotConnectorArrayIndex":        0,
  #               "mode":                             "callback",           # "copyLocal", "copyLocalIfPositive", "localSetIfAboveThreshold" or "communicate"
  #               "fromDofNosNumbering":              "local",
  #               "toDofNosNumbering":                "local",
  #               "dofsMapping":                      None,
  #               "inputDofs":                        list(range(variables.n_muscle_spindles)),
  #               "outputDofs":                       [list(range(variables.n_muscle_spindles)) for _ in range(5)],   # [0,1,...,n_muscle_spindles]
  #               "callback":                         variables.callback_muscle_spindles_input,
  #               #"thresholdValue":                   20,                  # if mode is "localSetIfAboveThreshold", this is the threshold, if the value is above it, set the value `valueToSet`
  #               #"valueToSet":                       20,                  # if mode is "localSetIfAboveThreshold", this is the value to set the target dof to, if the source dof is above thresholdValue.
  #             }
  #           ],
  #           "afterComputation":  None,
  # 
  #           # map from motoneuronMesh to stimulated nodes
  #           "MapDofs": {
  #             "description":                "motoneurons->stimulated nodes",  # description that will be shown in solver structure visualization
  #             "nAdditionalFieldVariables":  1,                              # number of additional field variables that are defined by this object. They have 1 component, use the templated function space and mesh given by meshName.
  #             "additionalSlotNames":        ["mn"],
  #             "meshName":                   "motoneuronMesh",               # the mesh on which the additional field variables will be defined
  #             "beforeComputation":          None, 
  # 
  #             # mapping from motoneuronMesh which contains on every rank as many nodes as there are motoneurons to the 3D domain
  #             # map from motoneuronMesh (algebraics) to the fiber meshes (solution)
  #             "afterComputation": [                                        # transfer/mapping of dofs that will be performed before the computation of the nested solver
  #               {
  #                 "fromConnectorSlot":                "mn",                # source slot of the dofs mapping
  #                 "toConnectorSlots":                 "vm",                # target slot of the dofs mapping
  #                 "fromSlotConnectorArrayIndex":      0,                  
  #                 "toSlotConnectorArrayIndex":        get_fiber_index_in_motor_unit(fiber_index, motor_unit_no),      # which fiber in this motor unit
  #                 "mode":                             "localSetIfAboveThreshold",          # "copyLocal", "copyLocalIfPositive", "localSetIfAboveThreshold" or "communicate"
  #                 "fromDofNosNumbering":              "local",
  #                 "toDofNosNumbering":                "global",
  #                 "dofsMapping":                      {motor_unit_no: stimulation_node_nos},   # map from the motor unit to the stimulated node of the fiber mesh
  #                 "inputDofs":                        None,                # this option is only needed in mode "callback"
  #                 "outputDofs":                       None,                # this option is only needed in mode "callback"
  #                 "callback":                         None,                # this option is only needed in mode "callback"
  #                 "thresholdValue":                   20,                  # if mode is "localSetIfAboveThreshold", this is the threshold, if the value is above it, set the value `valueToSet`
  #                 "valueToSet":                       variables.vm_value_stimulated,       # if mode is "localSetIfAboveThreshold", this is the value to set the target dof to, if the source dof is above thresholdValue.
  #               }
  #             for motor_unit_no in range(variables.n_motor_units)
  #               for fiber_index in range(get_n_fibers_in_motor_unit(motor_unit_no))],     # iterate over all motor units and all fibers in every motor unit
  # 
  # 
  #             # actual fibers_emg simulation plus electro-mechanics solver
  #             "Coupling": {
  #               "description":            "fibers+mechanics simulation",    # description that will be shown in solver structure visualization
  #               "timeStepWidth":          variables.dt_elasticity,
  #               "logTimeStepWidthAsKey":  "dt_elasticity",
  #               "durationLogKey":         "duration_total",
  #               "timeStepOutputInterval": 1,
  #               "endTime":                variables.end_time,
  #               "connectedSlotsTerm1To2": {1:2},  # {4:2} #None,       # data transfer is configured using global option "connectedSlots"
  #               "connectedSlotsTerm2To1": None,       # transfer nothing back
  #               "Term1": {        # fibers_emg
  # 
  # 
  #                 # fibers_emg
  #                 "Coupling": {
  #                   "description":            "fibers_emg",    # description that will be shown in solver structure visualization
  # 
  #                   "timeStepWidth":          variables.dt_3D,  # 1e-1
  #                   "logTimeStepWidthAsKey":  "dt_3D",
  #                   "durationLogKey":         "duration_total",
  #                   "timeStepOutputInterval": 10,
  #                   "endTime":                variables.end_time,
  #                   "connectedSlotsTerm1To2": {0:0}, 
  #                   "connectedSlotsTerm2To1": [None],  
  #                   "Term1": {        # monodomain, fibers
  #                     "MultipleInstances": {
  #                       "logKey":                     "duration_subdomains_xy",
  #                       "ranksAllComputedInstances":  list(range(n_ranks)),
  #                       "nInstances":                 variables.n_subdomains_xy,
  #                       "instances": 
  #                       [{
  #                         "ranks": list(range(subdomain_coordinate_y*variables.n_subdomains_x + subdomain_coordinate_x, n_ranks, variables.n_subdomains_x*variables.n_subdomains_y)),
  #                         "StrangSplitting": {
  #                           #"numberTimeSteps": 1,
  #                           "timeStepWidth":          variables.dt_splitting,  # 1e-1
  #                           "logTimeStepWidthAsKey":  "dt_splitting",
  #                           "durationLogKey":         "duration_monodomain",
  #                           "timeStepOutputInterval": 100,
  #                           "endTime":                variables.dt_splitting,
  #                           "connectedSlotsTerm1To2": {0:0},   # transfer slot 0 = state Vm from Term1 (CellML) to Term2 (Diffusion)
  #                           "connectedSlotsTerm2To1": {0:0},   # transfer the same back
  # 
  #                           "Term1": {      # CellML, i.e. reaction term of Monodomain equation
  #                             "MultipleInstances": {
  #                               "logKey":             "duration_subdomains_z",
  #                               "nInstances":         n_fibers_in_subdomain_x(subdomain_coordinate_x)*n_fibers_in_subdomain_y(subdomain_coordinate_y),
  #                               "instances": 
  #                               [{
  #                                 "ranks":                          list(range(variables.n_subdomains_z)),    # these rank nos are local nos to the outer instance of MultipleInstances, i.e. from 0 to number of ranks in z direction
  #                                 "Heun" : {
  #                                   "timeStepWidth":                variables.dt_0D,  # 5e-5
  #                                   "logTimeStepWidthAsKey":        "dt_0D",
  #                                   "durationLogKey":               "duration_0D",
  #                                   "initialValues":                [],
  #                                   "timeStepOutputInterval":       1e4,
  #                                   "inputMeshIsGlobal":            True,
  #                                   "checkForNanInf":               False,
  #                                   "dirichletBoundaryConditions":  {},
  #                                   "dirichletOutputFilename":      None,                # filename for a vtp file that contains the Dirichlet boundary condition nodes and their values, set to None to disable
  #                                   "nAdditionalFieldVariables":    0,
  #                                   "additionalSlotNames":          [],
  # 
  #                                   "CellML" : {
  #                                     "modelFilename":                          variables.cellml_file,                          # input C++ source file or cellml XML file
  #                                     "initializeStatesToEquilibrium":          False,                                          # if the equilibrium values of the states should be computed before the simulation starts
  #                                     "initializeStatesToEquilibriumTimestepWidth": 1e-4,                                       # if initializeStatesToEquilibrium is enable, the timestep width to use to solve the equilibrium equation
  # 
  #                                     # optimization parameters
  #                                     "optimizationType":                       variables.optimization_type,                    # "vc", "simd", "openmp" type of generated optimizated source file
  #                                     "approximateExponentialFunction":         variables.approximate_exponential_function,     # if optimizationType is "vc", whether the exponential function exp(x) should be approximate by (1+x/n)^n with n=1024
  #                                     "compilerFlags":                          "-fPIC -O3 -march=native -shared ",             # compiler flags used to compile the optimized model code
  #                                     "maximumNumberOfThreads":                 0,                                              # if optimizationType is "openmp", the maximum number of threads to use. Default value 0 means no restriction.
  #                                     #"libraryFilename":                        "lib/lib.so",                                   # manually compiled library
  # 
  #                                     # stimulation callbacks
  #                                     "setSpecificStatesFunction":              None,                                           # callback function that sets states like Vm, activation can be implemented by using this method and directly setting Vm values, or by using setParameters/setSpecificParameters
  #                                     "setSpecificStatesCallInterval":          0,                                              # 0 means disabled
  #                                     "setSpecificStatesCallFrequency":         0,                                              # set_specific_states should be called variables.stimulation_frequency times per ms
  #                                     "setSpecificStatesFrequencyJitter":       None,                                           # random value to add or substract to setSpecificStatesCallFrequency every stimulation, this is to add random jitter to the frequency
  #                                     "setSpecificStatesRepeatAfterFirstCall":  0.01,                                           # [ms] simulation time span for which the setSpecificStates callback will be called after a call was triggered
  #                                     "setSpecificStatesCallEnableBegin":       0,                                              # [ms] first time when to call setSpecificStates
  #                                     "additionalArgument":                     fiber_no,
  # 
  #                                     "mappings":                               variables.mappings,                             # mappings between parameters and algebraics/constants and between outputConnectorSlots and states, algebraics or parameters, they are defined in helper.py
  #                                     "parametersInitialValues":                variables.parameters_initial_values,            #[0.0, 1.0],      # initial values for the parameters: I_Stim, l_hs
  # 
  #                                     "meshName":                               "MeshFiber_{}".format(fiber_no),
  #                                     "stimulationLogFilename":                 "out/" + variables.scenario_name + "/stimulation.log",
  # 
  #                                     # output writer for states, algebraics and parameters                
  #                                     "OutputWriter" : [
  #                                       {"format": "Paraview", "outputInterval": (int)(1./variables.dt_0D*variables.output_timestep_0D_states), "filename": "out/" + variables.scenario_name + "/3_0D_all", "binary": True, "fixedFormat": False, "combineFiles": True, "fileNumbering": "incremental"}
  #                                     ] if variables.states_output else []
  #                                   }, 
  #                                   "OutputWriter" : [
  #                                     {"format": "Paraview", "outputInterval": 1, "filename": "out/" + variables.scenario_name + "/0D_states({},{})".format(fiber_in_subdomain_coordinate_x,fiber_in_subdomain_coordinate_y), "binary": True, "fixedFormat": False, "combineFiles": True, "fileNumbering": "incremental"}
  #                                   ] if variables.states_output else []
  # 
  #                                 },
  #                               } for fiber_in_subdomain_coordinate_y in range(n_fibers_in_subdomain_y(subdomain_coordinate_y)) \
  #                                   for fiber_in_subdomain_coordinate_x in range(n_fibers_in_subdomain_x(subdomain_coordinate_x)) \
  #                                     for fiber_no in [get_fiber_no(subdomain_coordinate_x, subdomain_coordinate_y, fiber_in_subdomain_coordinate_x, fiber_in_subdomain_coordinate_y)] \
  #                                       for motor_unit_no in [get_motor_unit_no(fiber_no)]],
  #                             }
  #                           },
  #                           "Term2": {     # Diffusion
  #                             "MultipleInstances": {
  #                               "nInstances": n_fibers_in_subdomain_x(subdomain_coordinate_x)*n_fibers_in_subdomain_y(subdomain_coordinate_y),
  #                               "instances": 
  #                               [{
  #                                 "ranks":                         list(range(variables.n_subdomains_z)),    # these rank nos are local nos to the outer instance of MultipleInstances, i.e. from 0 to number of ranks in z direction
  #                                 "ImplicitEuler" : {     # include both CrankNicolson and ImplicitEuler in the settings such that both variants in C++ file are possible
  #                                   "initialValues":               [],
  #                                   #"numberTimeSteps":            1,
  #                                   "timeStepWidth":               variables.dt_1D,  # 1e-5
  #                                   "timeStepWidthRelativeTolerance": 1e-10,
  #                                   "logTimeStepWidthAsKey":       "dt_1D",
  #                                   "durationLogKey":              "duration_1D",
  #                                   "timeStepOutputInterval":      1e4,
  #                                   "dirichletBoundaryConditions": {},            # old Dirichlet BC that are not used in FastMonodomainSolver: {0: -75.0036, -1: -75.0036},
  #                                   "dirichletOutputFilename":     None,          # filename for a vtp file that contains the Dirichlet boundary condition nodes and their values, set to None to disable
  #                                   "inputMeshIsGlobal":           True,
  #                                   "solverName":                  "implicitSolver",
  #                                   "checkForNanInf":              False,
  #                                   "nAdditionalFieldVariables":   3,
  #                                   "additionalSlotNames":         ["stress", "alpha"],
  #                                   "FiniteElementMethod" : {
  #                                     "inputMeshIsGlobal":         True,
  #                                     "meshName":                  "MeshFiber_{}".format(fiber_no),
  #                                     "prefactor":                 get_diffusion_prefactor(fiber_no, motor_unit_no),  # resolves to Conductivity / (Am * Cm)
  #                                     "solverName":                "implicitSolver",
  #                                     "slotName":                  "",
  #                                   },
  #                                   "OutputWriter" : [
  #                                     #{"format": "Paraview", "outputInterval": int(1./variables.dt_1D*variables.output_timestep), "filename": "out/fiber_"+str(fiber_no), "binary": True, "fixedFormat": False, "combineFiles": True},
  #                                     #{"format": "Paraview", "outputInterval": 1./variables.dt_1D*variables.output_timestep, "filename": "out/fiber_"+str(i)+"_txt", "binary": False, "fixedFormat": False},
  #                                     #{"format": "ExFile", "filename": "out/fiber_"+str(i), "outputInterval": 1./variables.dt_1D*variables.output_timestep, "sphereSize": "0.02*0.02*0.02"},
  #                                     #{"format": "PythonFile", "filename": "out/fiber_"+str(i), "outputInterval": 1./variables.dt_1D*variables.output_timestep, "binary":True, "onlyNodalValues":True},
  #                                   ]
  #                                 },
  #                                 "CrankNicolson": {    # include both CrankNicolson and ImplicitEuler in the settings such that both variants in C++ file are possible
  #                                   "initialValues":               [], 
  #                                   #"numberTimeSteps":            1,
  #                                   "timeStepWidth":               variables.dt_1D,  # 1e-5
  #                                   "timeStepWidthRelativeTolerance": 1e-10,
  #                                   "logTimeStepWidthAsKey":       "dt_1D",
  #                                   "durationLogKey":              "duration_1D",
  #                                   "timeStepOutputInterval":      1e4,
  #                                   "dirichletBoundaryConditions": {},            # old Dirichlet BC that are not used in FastMonodomainSolver: {0: -75.0036, -1: -75.0036},
  #                                   "dirichletOutputFilename":     None,          # filename for a vtp file that contains the Dirichlet boundary condition nodes and their values, set to None to disable
  #                                   "inputMeshIsGlobal":           True,
  #                                   "solverName":                  "implicitSolver",
  #                                   "checkForNanInf":              False,
  #                                   "nAdditionalFieldVariables":   2,
  #                                   "additionalSlotNames":         ["stress", "alpha"], 
  #                                   "FiniteElementMethod" : { 
  #                                     "inputMeshIsGlobal":         True,
  #                                     "meshName":                  "MeshFiber_{}".format(fiber_no),
  #                                     "prefactor":                 get_diffusion_prefactor(fiber_no, motor_unit_no),  # resolves to Conductivity / (Am * Cm)
  #                                     "solverName":                "implicitSolver",
  #                                     "slotName":                  "", 
  #                                   },  
  #                                   "OutputWriter" : [ 
  #                                   ]   
  #                                 },   
  #                               } for fiber_in_subdomain_coordinate_y in range(n_fibers_in_subdomain_y(subdomain_coordinate_y)) \
  #                                   for fiber_in_subdomain_coordinate_x in range(n_fibers_in_subdomain_x(subdomain_coordinate_x)) \
  #                                     for fiber_no in [get_fiber_no(subdomain_coordinate_x, subdomain_coordinate_y, fiber_in_subdomain_coordinate_x, fiber_in_subdomain_coordinate_y)] \
  #                                       for motor_unit_no in [get_motor_unit_no(fiber_no)]],
  #                               "OutputWriter" : [
  #                                 {"format": "Paraview", "outputInterval": int(1./variables.dt_splitting*variables.output_timestep_fibers), "filename": "out/" + variables.scenario_name + "/4_fibers", "binary": True, "fixedFormat": False, "combineFiles": True, "fileNumbering": "incremental"}
  #                               ]
  #                             },
  #                           },
  #                         }
  #                       } if (subdomain_coordinate_x,subdomain_coordinate_y) == (variables.own_subdomain_coordinate_x,variables.own_subdomain_coordinate_y) else None
  #                       for subdomain_coordinate_y in range(variables.n_subdomains_y)
  #                           for subdomain_coordinate_x in range(variables.n_subdomains_x)]
  #                     },
  #                     "fiberDistributionFile":    variables.fiber_distribution_file,   # for FastMonodomainSolver, e.g. MU_fibre_distribution_3780.txt
  #                     "firingTimesFile":          variables.firing_times_file,         # for FastMonodomainSolver, e.g. MU_firing_times_real.txt
  #                     "onlyComputeIfHasBeenStimulated": False,                          # only compute fibers after they have been stimulated for the first time
  #                     "disableComputationWhenStatesAreCloseToEquilibrium": True,       # optimization where states that are close to their equilibrium will not be computed again
  #                     "valueForStimulatedPoint":  variables.vm_value_stimulated,       # to which value of Vm the stimulated node should be set
  #                     "neuromuscularJunctionRelativeSize": 0.1,                        # range where the neuromuscular junction is located around the center, relative to fiber length. The actual position is draws randomly from the interval [0.5-s/2, 0.5+s/2) with s being this option. 0 means sharply at the center, 0.1 means located approximately at the center, but it can vary 10% in total between all fibers.
  #                     "generateGPUSource":        True,                                # (set to True) only effective if optimizationType=="gpu", whether the source code for the GPU should be generated. If False, an existing source code file (which has to have the correct name) is used and compiled, i.e. the code generator is bypassed. This is useful for debugging, such that you can adjust the source code yourself. (You can also add "-g -save-temps " to compilerFlags under CellMLAdapter)
  #                     "useSinglePrecision":       False,                               # only effective if optimizationType=="gpu", whether single precision computation should be used on the GPU. Some GPUs have poor double precision performance. Note, this drastically increases the error and, in consequence, the timestep widths should be reduced.
  #                     #"preCompileCommand":        "bash -c 'module load argon-tesla/gcc/11-20210110-openmp; module list; gcc --version",     # only effective if optimizationType=="gpu", system command to be executed right before the compilation
  #                     #"postCompileCommand":       "'",   # only effective if optimizationType=="gpu", system command to be executed right after the compilation
  #                   },
  #                   "Term2": {        # Bidomain, EMG
  #                     "OutputSurface": {
  #                       "OutputWriter": [
  #                         {"format": "Paraview", "outputInterval": int(1./variables.dt_3D*variables.output_timestep), "filename": "out/" + variables.scenario_name + "/4_surface_emg", "binary": True, "fixedFormat": False, "combineFiles": True, "fileNumbering": "incremental",},
  #                       ] if variables.enable_surface_emg else [],
  #                       "face":                     ["1+"],              # which faces of the 3D mesh should be written into the 2D mesh
  #                       "samplingPoints":           None,                # the electrode positions, they are created in the helper.py script
  #                       "enableCsvFile":            False,               # if the values at the sampling points should be written to csv files
  #                       "enableVtpFile":            False,               # if the values at the sampling points should be written to vtp files
  #                       "enableGeometryInCsvFile":  False,               # if the csv output file should contain geometry of the electrodes in every time step. This increases the file size and only makes sense if the geometry changed throughout time, i.e. when computing with contraction
  # 
  #                       "StaticBidomainSolver": {
  #                         "timeStepWidth":          variables.dt_3D,
  #                         "timeStepOutputInterval": 50,
  #                         "durationLogKey":         "duration_bidomain",
  #                         "solverName":             "activationSolver",
  #                         "initialGuessNonzero":    variables.emg_initial_guess_nonzero,
  #                         "slotNames":             [],
  #                         "nAdditionalFieldVariables": 0,                # number of additional field variable that can be used for additional output
  # 
  #                         "PotentialFlow": {
  #                           "FiniteElementMethod" : {
  #                             "meshName":           "3Dmesh",
  #                             "solverName":         "potentialFlowSolver",
  #                             "prefactor":          1.0,
  #                             "dirichletBoundaryConditions": variables.potential_flow_dirichlet_bc,
  #                             "dirichletOutputFilename":     None,                # filename for a vtp file that contains the Dirichlet boundary condition nodes and their values, set to None to disable
  #                             "neumannBoundaryConditions":   [],
  #                             "inputMeshIsGlobal":  True,
  #                             "slotName":           "",
  #                           },
  #                         },
  #                         "Activation": {
  #                           "FiniteElementMethod" : {
  #                             "meshName":           "3Dmesh",
  #                             "solverName":         "activationSolver",
  #                             "prefactor":          1.0,
  #                             "inputMeshIsGlobal":  True,
  #                             "dirichletBoundaryConditions": {},
  #                             "dirichletOutputFilename":     None,                # filename for a vtp file that contains the Dirichlet boundary condition nodes and their values, set to None to disable
  #                             "neumannBoundaryConditions":   [],
  #                             "slotName":           "",
  #                             "diffusionTensor": [[      # sigma_i, fiber direction is (1,0,0), one list item = same tensor for all elements, multiple list items = a different tensor for each element
  #                               8.93, 0, 0,
  #                               0, 0.893, 0,
  #                               0, 0, 0.893
  #                             ]],
  #                             "extracellularDiffusionTensor": [[      # sigma_e, one list item = same tensor for all elements, multiple list items = a different tensor for each element
  #                               6.7, 0, 0,
  #                               0, 6.7, 0,
  #                               0, 0, 6.7,
  #                             ]],
  #                           },
  #                         },
  #                         "OutputWriter" : [
  #                           {"format": "Paraview", "outputInterval": int(1./variables.dt_3D*variables.output_timestep_3D_emg), "filename": "out/" + variables.scenario_name + "/4_hd_emg", "binary": True, "fixedFormat": False, "combineFiles": True, "fileNumbering": "incremental"},
  #                         ]
  #                       }
  #                     },
  #                   }
  #                 }    
  #               },
  #               "Term2": {        # solid mechanics
  #                 "MuscleContractionSolver": {
  #                   "numberTimeSteps":              1,                         # only use 1 timestep per interval
  #                   "timeStepOutputInterval":       100,                       # do not output time steps
  #                   "Pmax":                         variables.pmax,            # maximum PK2 active stress
  #                   "enableForceLengthRelation":    True,                      # if the factor f_l(λ_f) modeling the force-length relation (as in Heidlauf2013) should be multiplied. Set to false if this relation is already considered in the CellML model.
  #                   "lambdaDotScalingFactor":       1.0,                       # scaling factor for the output of the lambda dot slot, i.e. the contraction velocity. Use this to scale the unit-less quantity to, e.g., micrometers per millisecond for the subcellular model.
  #                   "slotNames":                    ["m_lda", "m_ldot", "m_g_in", "m_T", "m_ux", "m_uy", "m_uz"],  # slot names of the data connector slots: lambda, lambdaDot, gamma, traction
  #                   "OutputWriter" : [
  #                     {"format": "Paraview", "outputInterval": int(1./variables.dt_elasticity*variables.output_timestep_elasticity), "filename": "out/" + variables.scenario_name + "/4_mechanics_3D", "binary": True, "fixedFormat": False, "onlyNodalValues":True, "combineFiles": True, "fileNumbering": "incremental"},
  #                   ],
  #                   "mapGeometryToMeshes":          ["3Dmesh","3DFatMesh"],    # the mesh names of the meshes that will get the geometry transferred
  #                   "reverseMappingOrder":          True,                      # if the mapping target->own mesh should be used instead of own->target mesh. This gives better results in some cases.
  #                   "dynamic":                      variables.dynamic,                      # if the dynamic solid mechanics solver should be used, else it computes the quasi-static problem
  # 
  #                   # the actual solid mechanics solver, this is either "DynamicHyperelasticitySolver" or "HyperelasticitySolver", depending on the value of "dynamic"
  #                   "DynamicHyperelasticitySolver": {
  #                     "timeStepWidth":              variables.dt_elasticity,           # time step width 
  #                     "durationLogKey":             "duration_mechanics",               # key to find duration of this solver in the log file
  #                     "timeStepOutputInterval":     1,                         # how often the current time step should be printed to console
  # 
  #                     "materialParameters":         variables.material_parameters,  # material parameters of the Mooney-Rivlin material
  #                     "density":                    variables.rho,             # density of the material
  #                     "dampingFactor":              variables.damping_factor,  # factor for velocity dependent damping
  #                     "displacementsScalingFactor": 1.0,                       # scaling factor for displacements, only set to sth. other than 1 only to increase visual appearance for very small displacements
  #                     "residualNormLogFilename":    "out/"+variables.scenario_name+"/4_main_log_residual_norm.txt",   # log file where residual norm values of the nonlinear solver will be written
  #                     "useAnalyticJacobian":        True,                      # whether to use the analytically computed jacobian matrix in the nonlinear solver (fast)
  #                     "useNumericJacobian":         False,                     # whether to use the numerically computed jacobian matrix in the nonlinear solver (slow), only works with non-nested matrices, if both numeric and analytic are enable, it uses the analytic for the preconditioner and the numeric as normal jacobian
  # 
  #                     "dumpDenseMatlabVariables":   False,                     # whether to have extra output of matlab vectors, x,r, jacobian matrix (very slow)
  #                     # if useAnalyticJacobian,useNumericJacobian and dumpDenseMatlabVariables all all three true, the analytic and numeric jacobian matrices will get compared to see if there are programming errors for the analytic jacobian
  # 
  #                     # mesh
  #                     "inputMeshIsGlobal":          True,                     # boundary conditions and initial values are given as global numbers (every process has all information)
  #                     "meshName":                   ["3Dmesh_elasticity_quadratic", "3DFatMesh_elasticity_quadratic"],       # name of the 3D mesh, it is defined under "Meshes" at the beginning of this config
  #                     "fiberMeshNames":             [],                       # fiber meshes that will be used to determine the fiber direction
  #                     "fiberDirection":             [0,0,1],                  # if fiberMeshNames is empty, directly set the constant fiber direction, in element coordinate system
  # 
  #                     # solving
  #                     "solverName":                 "mechanicsSolver",         # name of the nonlinear solver configuration, it is defined under "Solvers" at the beginning of this config
  #                     #"loadFactors":                [0.5, 1.0],                # load factors for every timestep
  #                     "loadFactors":                [],                        # no load factors, solve problem directly
  #                     "loadFactorGiveUpThreshold":  0.25,                       # a threshold for the load factor, when to abort the solve of the current time step. The load factors are adjusted automatically if the nonlinear solver diverged. If the load factors get too small, it aborts the solve.
  #                     "scaleInitialGuess":          False,                     # when load stepping is used, scale initial guess between load steps a and b by sqrt(a*b)/a. This potentially reduces the number of iterations per load step (but not always).
  #                     "nNonlinearSolveCalls":       1,                         # how often the nonlinear solve should be repeated
  # 
  #                     # boundary and initial conditions
  #                     "dirichletBoundaryConditions": variables.main_elasticity_dirichlet_bc,   # the initial Dirichlet boundary conditions that define values for displacements u and velocity v
  #                     "neumannBoundaryConditions":   variables.main_elasticity_neumann_bc,     # Neumann boundary conditions that define traction forces on surfaces of elements
  #                     "divideNeumannBoundaryConditionValuesByTotalArea": True,            # if the given Neumann boundary condition values under "neumannBoundaryConditions" are total forces instead of surface loads and therefore should be scaled by the surface area of all elements where Neumann BC are applied
  #                     "updateDirichletBoundaryConditionsFunction": update_dirichlet_boundary_conditions_helper,                  # function that updates the dirichlet BCs while the simulation is running
  #                     "updateDirichletBoundaryConditionsFunctionCallInterval": 1,         # every which step the update function should be called, 1 means every time step
  #                     "updateNeumannBoundaryConditionsFunction":   update_neumann_boundary_conditions_helper,                    # function that updates the Neumann BCs while the simulation is running
  #                     "updateNeumannBoundaryConditionsFunctionCallInterval": 1,           # every which step the update function should be called, 1 means every time step
  # 
  # 
  #                     "initialValuesDisplacements":  [[0.0,0.0,0.0] for _ in range(variables.n_points_global_composite_mesh)],     # the initial values for the displacements, vector of values for every node [[node1-x,y,z], [node2-x,y,z], ...]
  #                     "initialValuesVelocities":     [[0.0,0.0,0.0] for _ in range(variables.n_points_global_composite_mesh)],     # the initial values for the velocities, vector of values for every node [[node1-x,y,z], [node2-x,y,z], ...]
  #                     "extrapolateInitialGuess":     True,                                # if the initial values for the dynamic nonlinear problem should be computed by extrapolating the previous displacements and velocities
  #                     "constantBodyForce":           variables.main_constant_body_force,       # a constant force that acts on the whole body, e.g. for gravity
  # 
  #                     "dirichletOutputFilename":     "out/"+variables.scenario_name+"/4_main_dirichlet_boundary_conditions",     # output filename for the dirichlet boundary conditions, set to "" to have no output
  #                     "totalForceLogFilename":       "out/"+variables.scenario_name+"/4_main_tendon_force.csv",              # filename of a log file that will contain the total (bearing) forces and moments at the top and bottom of the volume
  #                     "totalForceLogOutputInterval":       10,                                  # output interval when to write the totalForceLog file
  #                     "totalForceBottomElementNosGlobal":  [j*nx + i for j in range(ny) for i in range(nx)],                  # global element nos of the bottom elements used to compute the total forces in the log file totalForceLogFilename
  #                     "totalForceTopElementNosGlobal":     [(nz-1)*ny*nx + j*nx + i for j in range(ny) for i in range(nx)],   # global element nos of the top elements used to compute the total forces in the log file totalForceTopElementsGlobal
  # 
  #                     # define which file formats should be written
  #                     # 1. main output writer that writes output files using the quadratic elements function space. Writes displacements, velocities and PK2 stresses.
  #                     "OutputWriter" : [
  # 
  #                       # Paraview files
  #                       {"format": "Paraview", "outputInterval": int(1./variables.dt_elasticity*variables.output_timestep_elasticity), "filename": "out/"+variables.scenario_name+"/4_displacements", "binary": True, "fixedFormat": False, "onlyNodalValues":True, "combineFiles":True, "fileNumbering": "incremental"},
  # 
  #                       # Python callback function "postprocess"
  #                       #{"format": "PythonCallback", "outputInterval": 1, "callback": postprocess, "onlyNodalValues":True, "filename": ""},
  #                     ],
  #                     # 2. additional output writer that writes also the hydrostatic pressure
  #                     "pressure": {   # output files for pressure function space (linear elements), contains pressure values, as well as displacements and velocities
  #                       "OutputWriter" : [
  #                         #{"format": "Paraview", "outputInterval": int(1./variables.dt_elasticity*variables.output_timestep_elasticity), "filename": "out/"+variables.scenario_name+"/4_pressure", "binary": True, "fixedFormat": False, "onlyNodalValues":True, "combineFiles":True, "fileNumbering": "incremental"},
  #                       ]
  #                     },
  #                     # 3. additional output writer that writes virtual work terms
  #                     "dynamic": {    # output of the dynamic solver, has additional virtual work values 
  #                       "OutputWriter" : [   # output files for displacements function space (quadratic elements)
  #                         #{"format": "Paraview", "outputInterval": 1, "filename": "out/"+variables.scenario_name+"/dynamic", "binary": True, "fixedFormat": False, "onlyNodalValues":True, "combineFiles":True, "fileNumbering": "incremental"},
  #                         #{"format": "Paraview", "outputInterval": int(1./variables.dt_elasticity*variables.output_timestep_elasticity), "filename": "out/"+variables.scenario_name+"/4_virtual_work", "binary": True, "fixedFormat": False, "onlyNodalValues":True, "combineFiles":True, "fileNumbering": "incremental"},
  #                       ],
  #                     },
  #                     # 4. output writer for debugging, outputs files after each load increment, the geometry is not changed but u and v are written
  #                     "LoadIncrements": {   
  #                       "OutputWriter" : [
  #                         #{"format": "Paraview", "outputInterval": int(1./variables.dt_elasticity*variables.output_timestep_elasticity), "filename": "out/"+variables.scenario_name+"/4_load_increments", "binary": False, "fixedFormat": False, "onlyNodalValues":True, "combineFiles":True, "fileNumbering": "incremental"},
  #                       ]
  #                     },
  #                   },
  # 
  #                   "HyperelasticitySolver": {
  #                     "durationLogKey":             "duration_mechanics",               # key to find duration of this solver in the log file
  #                     "timeStepOutputInterval":     1,                         # how often the current time step should be printed to console
  # 
  #                     "materialParameters":         variables.material_parameters,  # material parameters of the Mooney-Rivlin material
  #                     "displacementsScalingFactor": 1.0,                       # scaling factor for displacements, only set to sth. other than 1 only to increase visual appearance for very small displacements
  #                     "residualNormLogFilename":    "out/"+variables.scenario_name+"/4_main_log_residual_norm.txt",   # log file where residual norm values of the nonlinear solver will be written
  #                     "useAnalyticJacobian":        True,                      # whether to use the analytically computed jacobian matrix in the nonlinear solver (fast)
  #                     "useNumericJacobian":         False,                     # whether to use the numerically computed jacobian matrix in the nonlinear solver (slow), only works with non-nested matrices, if both numeric and analytic are enable, it uses the analytic for the preconditioner and the numeric as normal jacobian
  # 
  #                     "dumpDenseMatlabVariables":   False,                     # whether to have extra output of matlab vectors, x,r, jacobian matrix (very slow)
  #                     # if useAnalyticJacobian,useNumericJacobian and dumpDenseMatlabVariables all all three true, the analytic and numeric jacobian matrices will get compared to see if there are programming errors for the analytic jacobian
  # 
  #                     # mesh
  #                     "inputMeshIsGlobal":          True,                     # boundary conditions and initial values are given as global numbers (every process has all information)
  #                     "meshName":                   ["3Dmesh_elasticity_quadratic", "3DFatMesh_elasticity_quadratic"],       # name of the 3D mesh, it is defined under "Meshes" at the beginning of this config
  #                     "fiberMeshNames":             [],                       # fiber meshes that will be used to determine the fiber direction
  #                     "fiberDirection":             [0,0,1],                  # if fiberMeshNames is empty, directly set the constant fiber direction, in element coordinate system
  # 
  #                     # solving
  #                     "solverName":                 "mechanicsSolver",         # name of the nonlinear solver configuration, it is defined under "Solvers" at the beginning of this config
  #                     #"loadFactors":                [0.5, 1.0],                # load factors for every timestep
  #                     "loadFactors":                [],                        # no load factors, solve problem directly
  #                     "loadFactorGiveUpThreshold":  0.25,                       # a threshold for the load factor, when to abort the solve of the current time step. The load factors are adjusted automatically if the nonlinear solver diverged. If the load factors get too small, it aborts the solve.
  #                     "scaleInitialGuess":          False,                     # when load stepping is used, scale initial guess between load steps a and b by sqrt(a*b)/a. This potentially reduces the number of iterations per load step (but not always).
  #                     "nNonlinearSolveCalls":       1,                         # how often the nonlinear solve should be repeated
  # 
  #                     # boundary and initial conditions
  #                     "dirichletBoundaryConditions": variables.main_elasticity_dirichlet_bc,   # the initial Dirichlet boundary conditions that define values for displacements u and velocity v
  #                     "neumannBoundaryConditions":   variables.main_elasticity_neumann_bc,     # Neumann boundary conditions that define traction forces on surfaces of elements
  #                     "divideNeumannBoundaryConditionValuesByTotalArea": True,            # if the given Neumann boundary condition values under "neumannBoundaryConditions" are total forces instead of surface loads and therefore should be scaled by the surface area of all elements where Neumann BC are applied
  #                     "updateDirichletBoundaryConditionsFunction": None, #update_dirichlet_boundary_conditions_helper,                  # function that updates the dirichlet BCs while the simulation is running
  #                     "updateDirichletBoundaryConditionsFunctionCallInterval": 1,         # every which step the update function should be called, 1 means every time step
  #                     "updateNeumannBoundaryConditionsFunction":   update_neumann_boundary_conditions_helper,                    # function that updates the Neumann BCs while the simulation is running
  #                     "updateNeumannBoundaryConditionsFunctionCallInterval": 1,           # every which step the update function should be called, 1 means every time step
  # 
  # 
  #                     "initialValuesDisplacements":  [[0.0,0.0,0.0] for _ in range(variables.n_points_global_composite_mesh)],     # the initial values for the displacements, vector of values for every node [[node1-x,y,z], [node2-x,y,z], ...]
  #                     "extrapolateInitialGuess":     True,                                # if the initial values for the dynamic nonlinear problem should be computed by extrapolating the previous displacements and velocities
  #                     "constantBodyForce":           variables.main_constant_body_force,       # a constant force that acts on the whole body, e.g. for gravity
  # 
  #                     "dirichletOutputFilename":     "out/"+variables.scenario_name+"/4_main_dirichlet_boundary_conditions",     # output filename for the dirichlet boundary conditions, set to "" to have no output
  #                     "totalForceLogFilename":       "out/"+variables.scenario_name+"/4_main_tendon_force.csv",              # filename of a log file that will contain the total (bearing) forces and moments at the top and bottom of the volume
  #                     "totalForceLogOutputInterval":       10,                                  # output interval when to write the totalForceLog file
  #                     "totalForceBottomElementNosGlobal":  [j*nx + i for j in range(ny) for i in range(nx)],                  # global element nos of the bottom elements used to compute the total forces in the log file totalForceLogFilename
  #                     "totalForceTopElementNosGlobal":     [(nz-1)*ny*nx + j*nx + i for j in range(ny) for i in range(nx)],   # global element nos of the top elements used to compute the total forces in the log file totalForceTopElementsGlobal
  # 
  #                     # define which file formats should be written
  #                     # 1. main output writer that writes output files using the quadratic elements function space. Writes displacements, velocities and PK2 stresses.
  #                     "OutputWriter" : [
  # 
  #                       # Paraview files
  #                       {"format": "Paraview", "outputInterval": int(1./variables.dt_elasticity*variables.output_timestep_elasticity), "filename": "out/"+variables.scenario_name+"/4_displacements", "binary": True, "fixedFormat": False, "onlyNodalValues":True, "combineFiles":True, "fileNumbering": "incremental"},
  # 
  #                       # Python callback function "postprocess"
  #                       #{"format": "PythonCallback", "outputInterval": 1, "callback": postprocess, "onlyNodalValues":True, "filename": ""},
  #                     ],
  #                     # 2. additional output writer that writes also the hydrostatic pressure
  #                     "pressure": {   # output files for pressure function space (linear elements), contains pressure values, as well as displacements and velocities
  #                       "OutputWriter" : [
  #                         #{"format": "Paraview", "outputInterval": int(1./variables.dt_elasticity*variables.output_timestep_elasticity), "filename": "out/"+variables.scenario_name+"/4_pressure", "binary": True, "fixedFormat": False, "onlyNodalValues":True, "combineFiles":True, "fileNumbering": "incremental"},
  #                       ]
  #                     },
  #                     # 3. additional output writer that writes virtual work terms
  #                     "dynamic": {    # output of the dynamic solver, has additional virtual work values 
  #                       "OutputWriter" : [   # output files for displacements function space (quadratic elements)
  #                         #{"format": "Paraview", "outputInterval": 1, "filename": "out/"+variables.scenario_name+"/dynamic", "binary": True, "fixedFormat": False, "onlyNodalValues":True, "combineFiles":True, "fileNumbering": "incremental"},
  #                         #{"format": "Paraview", "outputInterval": int(1./variables.dt_elasticity*variables.output_timestep_elasticity), "filename": "out/"+variables.scenario_name+"/4_virtual_work", "binary": True, "fixedFormat": False, "onlyNodalValues":True, "combineFiles":True, "fileNumbering": "incremental"},
  #                       ],
  #                     },
  #                     # 4. output writer for debugging, outputs files after each load increment, the geometry is not changed but u and v are written
  #                     "LoadIncrements": {   
  #                       "OutputWriter" : [
  #                         #{"format": "Paraview", "outputInterval": int(1./variables.dt_elasticity*variables.output_timestep_elasticity), "filename": "out/"+variables.scenario_name+"/4_load_increments", "binary": False, "fixedFormat": False, "onlyNodalValues":True, "combineFiles":True, "fileNumbering": "incremental"},
  #                       ]
  #                     }
  #                   }
  #                 }
  #               }
  #             }
  #           }
  #         }
  #       }
  #     }
  #   }
  # }
}

print('##S CONFIG##')
print(config)
print('##E CONFIG##')

# stop timer and calculate how long parsing lasted
if rank_no == 0:
  t_stop_script = timeit.default_timer()
  print("Python config parsed in {:.1f}s.".format(t_stop_script - t_start_script))
