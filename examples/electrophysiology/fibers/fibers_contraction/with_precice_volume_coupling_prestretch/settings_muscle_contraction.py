# Multiple 1D fibers (monodomain) with 3D static Mooney Rivlin with active contraction term, on biceps geometry.
# This example uses precice to couple the Monodomain eq. on fibers with the 3D contraction, both inside opendihu.
# You need to run both executables at the same time.

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
  print("Warning: There is no variables file, e.g:\n ./muscle_contraction ../settings_muscle_contraction.py ramp.py\n")
  exit(0)

# -------------- begin user parameters ----------------
# -------------- end user parameters ----------------

# define command line arguments
parser = argparse.ArgumentParser(description='static_biceps_emg')
parser.add_argument('--scenario_name',                       help='The name to identify this run in the log.',   default=variables.scenario_name)
parser.add_argument('--n_subdomains', nargs=3,               help='Number of subdomains in x,y,z direction.',    type=int)
parser.add_argument('--n_subdomains_x', '-x',                help='Number of subdomains in x direction.',        type=int, default=variables.n_subdomains_x)
parser.add_argument('--n_subdomains_y', '-y',                help='Number of subdomains in y direction.',        type=int, default=variables.n_subdomains_y)
parser.add_argument('--n_subdomains_z', '-z',                help='Number of subdomains in z direction.',        type=int, default=variables.n_subdomains_z)
parser.add_argument('--paraview_output',                     help='Enable the paraview output writer.',          default=variables.paraview_output, action='store_true')
parser.add_argument('--adios_output',                        help='Enable the MegaMol/ADIOS output writer.',          default=variables.adios_output, action='store_true')
parser.add_argument('--fiber_file',                          help='The filename of the file that contains the fiber data.', default=variables.fiber_file)
parser.add_argument('--end_time', '--tend', '-t',            help='The end simulation time.',                    type=float, default=variables.end_time)
parser.add_argument('--output_timestep',                     help='The timestep for writing outputs.',           type=float, default=variables.output_timestep)
parser.add_argument('--dt_3D',                               help='The timestep for the 3D model, i.e. dynamic solid mechanics.', type=float, default=variables.dt_3D)
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
  print("dt_3D:           {:0.0e}, paraview_output: {}".format(variables.dt_3D, variables.paraview_output))
  print("output_timestep: {:0.0e}".format(variables.output_timestep))
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
  "scenarioName":          variables.scenario_name,
  "logFormat":             "csv",
  "solverStructureDiagramFile":     "out/" + variables.scenario_name + "/solver_structure_muscle_contraction.txt",     # output file of a diagram that shows data connection between solvers
  "mappingsBetweenMeshesLogFile":   "out/" + variables.scenario_name + "/mappings_between_meshes_muscle_contraction.txt",
  "regularization":        None,                         # None to disable or [tol,eps]. if the absolute determinant |J|=|det(F)| within an element is below tol, add eps*I to the matrix, to regularize the inversion of the nearly singular matrix
  "Meshes":                variables.meshes,
  "MappingsBetweenMeshes": {
    "3Dmesh": [
       {"name": "3Dmesh_quadratic", "xiTolerance": 0.35, "enableWarnings": False, "compositeUseOnlyInitializedMappings": True, "fixUnmappedDofs": True, "defaultValue": 0},
    ],
  },
  "connectedSlots": [
    ("p_lda", "lambda"),
    ("p_ldot","ldot"),
    ("p_gam", "gamma"),
    ("p_T",   "T"),
  ],
  "Solvers": {
    "precontractionMechanicsSolver": {   # solver for the preprocessing contraction simulation (static mechanics problem)
      "relativeTolerance":   1e-10,           # 1e-10 relative tolerance of the linear solver
      "absoluteTolerance":   1e-10,           # 1e-10 absolute tolerance of the residual of the linear solver
      "solverType":          "lu",            # type of the linear solver
      "preconditionerType":  "none",          # type of the preconditioner
      "maxIterations":       1e4,                                         # maximum number of iterations in the linear solver
      "snesMaxFunctionEvaluations": 1e8,                                  # maximum number of function iterations
      "snesMaxIterations":   34,              # maximum number of iterations in the nonlinear solver
      "snesRelativeTolerance": 2e-2,         # relative tolerance of the nonlinear solver
      "snesAbsoluteTolerance": 2e-2,         # absolute tolerance of the nonlinear solver
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
    "mechanicsSolver": {   # solver for the dynamic mechanics problem
      "relativeTolerance":  1e-10,          # 1e-10 relative tolerance of the linear solver
      "absoluteTolerance":  1e-10,          # 1e-10 absolute tolerance of the residual of the linear solver       
      "solverType":         "preonly",      # type of the linear solver: cg groppcg pipecg pipecgrr cgne nash stcg gltr richardson chebyshev gmres tcqmr fcg pipefcg bcgs ibcgs fbcgs fbcgsr bcgsl cgs tfqmr cr pipecr lsqr preonly qcg bicg fgmres pipefgmres minres symmlq lgmres lcd gcr pipegcr pgmres dgmres tsirm cgls
      "preconditionerType": "lu",           # type of the preconditioner
      "maxIterations":       1e4,           # maximum number of iterations in the linear solver
      "snesMaxFunctionEvaluations": 1e8,    # maximum number of function iterations
      "snesMaxIterations":   40,            # maximum number of iterations in the nonlinear solver
      "snesRelativeTolerance": 1e-5,       # relative tolerance of the nonlinear solver
      "snesAbsoluteTolerance": 1e-5,        # absolute tolerance of the nonlinear solver
      "snesLineSearchType": "l2",        # type of linesearch, possible values: "bt" "nleqerr" "basic" "l2" "cp" "ncglinear"
      "snesRebuildJacobianFrequency": 1,
      "dumpFilename":        "",
      "dumpFormat":          "matlab",
    }
  },
  
  # preprocessing (prescribe activation, precontraction, prestretch) and main solver
  "MultipleCoupling": {
    "timeStepWidth":          variables.end_time,
    "logTimeStepWidthAsKey":  "dt_artifical_contraction",
    "durationLogKey":         "duration_artifical_contraction",
    "timeStepOutputInterval": 1,
    "endTime":                variables.end_time,
    "deferInitialization":    True,       # initialize nested solvers only right before computation of first timestep
    "connectedSlotsTerm1To2": None,       # connect lambda to slot 0 and gamma to slot 2
    "connectedSlotsTerm2To1": None,       # transfer nothing back
    
    # prescribed activation
    "Term1": {
      "PrescribedValues": {
        "meshName":               "3Dmesh",      # reference to the fibers mesh
        "numberTimeSteps":        1,             # number of timesteps to call the callback functions subsequently, this is usually 1 for prescribed values, because it is enough to set the reaction term only once per time step
        "timeStepOutputInterval": 20,            # if the time step should be written to console, a value > 10 produces no output
        "slotNames":              ["p_lda", "p_gam"],
        
        # a list of field variables that will get values assigned in every timestep, by the provided callback function
        "fieldVariables1": [
          {"name": "lambda", "callback": variables.set_lambda_values},
          {"name": "gamma", "callback": variables.set_gamma_values},
        ],
        "fieldVariables2":     [],
        "additionalArgument":  None,         # a custom argument to the fieldVariables callback functions, this will be passed on as the last argument
        
        "OutputWriter" : [
          {"format": "Paraview", "outputInterval": 1, "filename": "out/" + variables.scenario_name + "/0_prescribed_stress", "binary": True, "fixedFormat": False, "combineFiles": True, "fileNumbering": "incremental"}
        ]
      }
    },
    # artifical contraction with constant gamma
    "Term2": {  
      "MuscleContractionSolver": {
        "numberTimeSteps":              1,                         # only use 1 timestep per interval
        "timeStepOutputInterval":       100,                       # do not output time steps
        "Pmax":                         variables.pmax,            # maximum PK2 active stress
        "enableForceLengthRelation":    False,                     # if the factor f_l(λ_f) modeling the force-length relation (as in Heidlauf2013) should be multiplied. Set to false if this relation is already considered in the CellML model.
        "lambdaDotScalingFactor":       1.0,                       # scaling factor for the output of the lambda dot slot, i.e. the contraction velocity. Use this to scale the unit-less quantity to, e.g., micrometers per millisecond for the subcellular model.
        "slotNames":                    ["p_lda", "p_ldot", "p_gam", "p_T"],
        "OutputWriter" : [
          {"format": "Paraview", "outputInterval": int(1./variables.dt_elasticity*variables.output_timestep_elasticity), "filename": "out/"+variables.scenario_name+"/1_precontraction_mechanics_3D", "binary": True, "fixedFormat": False, "onlyNodalValues":True, "combineFiles": True, "fileNumbering": "incremental"},
        ],
        "mapGeometryToMeshes":          "3Dmesh_quadratic",    # the mesh names of the meshes that will get the geometry transferred
        "reverseMappingOrder":          True,                      # if the mapping target->own mesh should be used instead of own->target mesh. This gives better results in some cases.
        "dynamic":                      False,                     # if the dynamic solid mechanics solver should be used, else it computes the quasi-static problem
        
        # the actual solid mechanics solver, this is either "DynamicHyperelasticitySolver" or "HyperelasticitySolver", depending on the value of "dynamic"
                                 
        "HyperelasticitySolver": {
          "durationLogKey":             "duration_precontraction",               # key to find duration of this solver in the log file
          
          "materialParameters":         variables.material_parameters,  # material parameters of the Mooney-Rivlin material
          "displacementsScalingFactor": 1.0,                       # scaling factor for displacements, only set to sth. other than 1 only to increase visual appearance for very small displacements
          "residualNormLogFilename":    "out/"+variables.scenario_name+"/1_log_precontraction_residual_norm.txt",   # log file where residual norm values of the nonlinear solver will be written
          "useAnalyticJacobian":        True,                      # whether to use the analytically computed jacobian matrix in the nonlinear solver (fast)
          "useNumericJacobian":         False,                     # whether to use the numerically computed jacobian matrix in the nonlinear solver (slow), only works with non-nested matrices, if both numeric and analytic are enable, it uses the analytic for the preconditioner and the numeric as normal jacobian
            
          "dumpDenseMatlabVariables":   False,                     # whether to have extra output of matlab vectors, x,r, jacobian matrix (very slow)
          # if useAnalyticJacobian,useNumericJacobian and dumpDenseMatlabVariables all all three true, the analytic and numeric jacobian matrices will get compared to see if there are programming errors for the analytic jacobian
          
          # mesh
          "inputMeshIsGlobal":          True,                     # boundary conditions and initial values are given as global numbers (every process has all information)
          "meshName":                   "3Dmesh_quadratic_precontraction",       # name of the 3D mesh, it is defined under "Meshes" at the beginning of this config
          "fiberMeshNames":             [],                       # fiber meshes that will be used to determine the fiber direction, there are no fibers in multidomain, so this is empty
          "fiberDirection":             [0,0,1],                  # if fiberMeshNames is empty, directly set the constant fiber direction, in element coordinate system

          # solving
          "solverName":                 "precontractionMechanicsSolver",         # name of the nonlinear solver configuration, it is defined under "Solvers" at the beginning of this config
          #"loadFactors":                [0.5, 1.0],                # load factors for every timestep
          "loadFactors":                list(np.logspace(-1,0,3)),        # load factors
          #"loadFactors":                [],                        # no load factors
          "loadFactorGiveUpThreshold":  1e-5,                      # a threshold for the load factor, when to abort the solve of the current time step. The load factors are adjusted automatically if the nonlinear solver diverged. If the load factors get too small, it aborts the solve.
          "scaleInitialGuess":          False,                      # when load stepping is used, scale initial guess between load steps a and b by sqrt(a*b)/a. This potentially reduces the number of iterations per load step (but not always).
          "nNonlinearSolveCalls":       1,                         # how often the nonlinear solve should be repeated
          
          # boundary and initial conditions
          "dirichletBoundaryConditions": variables.precontraction_elasticity_dirichlet_bc,   # the initial Dirichlet boundary conditions that define values for displacements u and velocity v
          "neumannBoundaryConditions":   variables.precontraction_elasticity_neumann_bc,     # Neumann boundary conditions that define traction forces on surfaces of elements
          "divideNeumannBoundaryConditionValuesByTotalArea": True,            # if the given Neumann boundary condition values under "neumannBoundaryConditions" are total forces instead of surface loads and therefore should be scaled by the surface area of all elements where Neumann BC are applied
          "updateDirichletBoundaryConditionsFunction": None,                  # function that updates the dirichlet BCs while the simulation is running
          "updateDirichletBoundaryConditionsFunctionCallInterval": 1,         # every which step the update function should be called, 1 means every time step
          
          "initialValuesDisplacements":  [[0.0,0.0,0.0] for _ in range(variables.n_points_global_mesh)],     # the initial values for the displacements, vector of values for every node [[node1-x,y,z], [node2-x,y,z], ...]
          "initialValuesVelocities":     [[0.0,0.0,0.0] for _ in range(variables.n_points_global_mesh)],     # the initial values for the velocities, vector of values for every node [[node1-x,y,z], [node2-x,y,z], ...]
          "extrapolateInitialGuess":     True,                                # if the initial values for the dynamic nonlinear problem should be computed by extrapolating the previous displacements and velocities
          "constantBodyForce":           variables.precontraction_constant_body_force,       # a constant force that acts on the whole body, e.g. for gravity
          
          "dirichletOutputFilename":    "out/"+variables.scenario_name+"/1_precontraction_dirichlet_boundary_conditions",     # output filename for the dirichlet boundary conditions, set to "" to have no output
          
          # define which file formats should be written
          # 1. main output writer that writes output files using the quadratic elements function space. Writes displacements, velocities and PK2 stresses.
          "OutputWriter" : [
            
            # Paraview files
            {"format": "Paraview", "outputInterval": 1, "filename": "out/"+variables.scenario_name+"/1_precontraction_u", "binary": True, "fixedFormat": False, "onlyNodalValues":True, "combineFiles":True, "fileNumbering": "incremental"},
          ],
          # 2. additional output writer that writes also the hydrostatic pressure
          "pressure": {   # output files for pressure function space (linear elements), contains pressure values, as well as displacements and velocities
            "OutputWriter" : [
              {"format": "Paraview", "outputInterval": 1, "filename": "out/"+variables.scenario_name+"/1_precontraction_p", "binary": True, "fixedFormat": False, "onlyNodalValues":True, "combineFiles":True, "fileNumbering": "incremental"},
            ]
          },
          # 4. output writer for debugging, outputs files after each load increment, the geometry is not changed but u and v are written
          "LoadIncrements": {   
            "OutputWriter" : [
              {"format": "Paraview", "outputInterval": 1, "filename": "out/"+variables.scenario_name+"/1_precontraction_load_increments", "binary": True, "fixedFormat": False, "onlyNodalValues":True, "combineFiles":True, "fileNumbering": "incremental"},
            ]
          }
        }
      }
    },
    # prestretch simulation
    "Term3": {
      "HyperelasticitySolver": {
        "durationLogKey":             "duration_prestretch",               # key to find duration of this solver in the log file
        
        "materialParameters":         variables.material_parameters,  # material parameters of the Mooney-Rivlin material
        "displacementsScalingFactor": 1.0,                       # scaling factor for displacements, only set to sth. other than 1 only to increase visual appearance for very small displacements
        "residualNormLogFilename":    "out/"+variables.scenario_name+"/2_log_prestretch_residual_norm.txt",   # log file where residual norm values of the nonlinear solver will be written
        "useAnalyticJacobian":        True,                      # whether to use the analytically computed jacobian matrix in the nonlinear solver (fast)
        "useNumericJacobian":         False,                     # whether to use the numerically computed jacobian matrix in the nonlinear solver (slow), only works with non-nested matrices, if both numeric and analytic are enable, it uses the analytic for the preconditioner and the numeric as normal jacobian
        "slotNames":                  ["m_ux", "m_uy", "m_uz"],  # slot names of the data connector slots, there are three slots, namely the displacement components ux, uy, uz
          
        "dumpDenseMatlabVariables":   False,                     # whether to have extra output of matlab vectors, x,r, jacobian matrix (very slow)
        # if useAnalyticJacobian,useNumericJacobian and dumpDenseMatlabVariables all all three true, the analytic and numeric jacobian matrices will get compared to see if there are programming errors for the analytic jacobian
        
        # mesh
        "inputMeshIsGlobal":          True,                     # boundary conditions and initial values are given as global numbers (every process has all information)
        "meshName":                   "3Dmesh_quadratic",       # name of the 3D mesh, it is defined under "Meshes" at the beginning of this config
        "fiberMeshNames":             [],                       # fiber meshes that will be used to determine the fiber direction, there are no fibers in multidomain, so this is empty
        "fiberDirection":             [0,0,1],                  # if fiberMeshNames is empty, directly set the constant fiber direction, in element coordinate system

        # solving
        "solverName":                 "prestretchMechanicsSolver",         # name of the nonlinear solver configuration, it is defined under "Solvers" at the beginning of this config
        #"loadFactors":                [0.5, 1.0],                # load factors for every timestep
        "loadFactors":                list(np.logspace(-1,0,3)),        # load factors
        #"loadFactors":                [],                        # no load factors
        "loadFactorGiveUpThreshold":  1e-5,                      # a threshold for the load factor, when to abort the solve of the current time step. The load factors are adjusted automatically if the nonlinear solver diverged. If the load factors get too small, it aborts the solve.
        "scaleInitialGuess":          False,                      # when load stepping is used, scale initial guess between load steps a and b by sqrt(a*b)/a. This potentially reduces the number of iterations per load step (but not always).
        "nNonlinearSolveCalls":       1,                         # how often the nonlinear solve should be repeated
        
        # boundary and initial conditions
        "dirichletBoundaryConditions": variables.prestretch_elasticity_dirichlet_bc,   # the initial Dirichlet boundary conditions that define values for displacements u and velocity v
        "neumannBoundaryConditions":   variables.prestretch_elasticity_neumann_bc,     # Neumann boundary conditions that define traction forces on surfaces of elements
        "divideNeumannBoundaryConditionValuesByTotalArea": True,            # if the given Neumann boundary condition values under "neumannBoundaryConditions" are total forces instead of surface loads and therefore should be scaled by the surface area of all elements where Neumann BC are applied
        "updateDirichletBoundaryConditionsFunction": None,                  # function that updates the dirichlet BCs while the simulation is running
        "updateDirichletBoundaryConditionsFunctionCallInterval": 1,         # every which step the update function should be called, 1 means every time step
        
        "initialValuesDisplacements":  [[0.0,0.0,0.0] for _ in range(variables.n_points_global_mesh)],     # the initial values for the displacements, vector of values for every node [[node1-x,y,z], [node2-x,y,z], ...]
        "initialValuesVelocities":     [[0.0,0.0,0.0] for _ in range(variables.n_points_global_mesh)],     # the initial values for the velocities, vector of values for every node [[node1-x,y,z], [node2-x,y,z], ...]
        "extrapolateInitialGuess":     True,                                # if the initial values for the dynamic nonlinear problem should be computed by extrapolating the previous displacements and velocities
        "constantBodyForce":           variables.prestretch_constant_body_force,       # a constant force that acts on the whole body, e.g. for gravity
        
        "dirichletOutputFilename":    "out/"+variables.scenario_name+"/2_prestretch_dirichlet_boundary_conditions",     # output filename for the dirichlet boundary conditions, set to "" to have no output
        
        # define which file formats should be written
        # 1. main output writer that writes output files using the quadratic elements function space. Writes displacements, velocities and PK2 stresses.
        "OutputWriter" : [
          
          # Paraview files
          {"format": "Paraview", "outputInterval": 1, "filename": "out/"+variables.scenario_name+"/2_prestretch_u", "binary": True, "fixedFormat": False, "onlyNodalValues":True, "combineFiles":True, "fileNumbering": "incremental"},
        ],
        # 2. additional output writer that writes also the hydrostatic pressure
        "pressure": {   # output files for pressure function space (linear elements), contains pressure values, as well as displacements and velocities
          "OutputWriter" : [
            {"format": "Paraview", "outputInterval": 1, "filename": "out/"+variables.scenario_name+"/2_prestretch_p", "binary": True, "fixedFormat": False, "onlyNodalValues":True, "combineFiles":True, "fileNumbering": "incremental"},
          ]
        },
        # 4. output writer for debugging, outputs files after each load increment, the geometry is not changed but u and v are written
        "LoadIncrements": {   
          "OutputWriter" : [
            {"format": "Paraview", "outputInterval": 1, "filename": "out/"+variables.scenario_name+"/2_prestretch_load_increments", "binary": True, "fixedFormat": False, "onlyNodalValues":True, "combineFiles":True, "fileNumbering": "incremental"},
          ]
        },
      }
    },
    
    # actual simulation
    "Term4": {
      "PreciceAdapterVolumeCoupling": {
        "timeStepOutputInterval":   100,                        # interval in which to display current timestep and time in console
        "timestepWidth":            1,                          # coupling time step width, must match the value in the precice config
        "couplingEnabled":          True,                       # if the precice coupling is enabled, if not, it simply calls the nested solver, for debugging
        "preciceConfigFilename":    "../precice_config.xml",    # the preCICE configuration file
        "preciceParticipantName":   "MuscleContraction",        # name of the own precice participant, has to match the name given in the precice xml config file
        "scalingFactor":            1,                          # a factor to scale the exchanged data, prior to communication
        "outputOnlyConvergedTimeSteps": True,                   # if the output writers should be called only after a time window of precice is complete, this means the timestep has converged
        
        "preciceData": [
          {
            "mode":                 "write",                    # mode is one of "read" or "write"
            "preciceDataName":      "Geometry",                 # name of the vector or scalar to transfer, as given in the precice xml settings file
            "preciceMeshName":      "MuscleContractionMesh",    # name of the precice coupling mesh, as given in the precice xml settings file
            "opendihuMeshName":     None,                       # extra specification of the opendihu mesh that is used for the initialization of the precice mapping. If None or "", the mesh of the field variable is used.
            "slotName":             None,                       # name of the existing slot of the opendihu data connector to which this variable is associated to (only relevant if not isGeometryField)
            "isGeometryField":      True,                       # if this is the geometry field of the mesh
          },
          {
            "mode":                 "read",                     # mode is one of "read" or "write"
            "preciceDataName":      "Gamma",                    # name of the vector or scalar to transfer, as given in the precice xml settings file
            "preciceMeshName":      "MuscleContractionMesh",    # name of the precice coupling mesh, as given in the precice xml settings file
            "opendihuMeshName":     None,                       # extra specification of the opendihu mesh that is used for the initialization of the precice mapping. If None or "", the mesh of the field variable is used.
            "slotName":             "gamma",                    # name of the existing slot of the opendihu data connector to which this variable is associated to (only relevant if not isGeometryField)
            "isGeometryField":      False,                      # if this is the geometry field of the mesh
          },
        ],
        
        "MuscleContractionSolver": {
          "mapGeometryToMeshes":          [],
          "numberTimeSteps":              1,                         # only use 1 timestep per interval
          "timeStepOutputInterval":       100,                       # do not output time steps
          "endTime":                      variables.end_time,        # end time of the simulation (if precice coupling is disabled)
          "Pmax":                         variables.pmax,            # maximum PK2 active stress
          "enableForceLengthRelation":    True,                      # if the factor f_l(λ_f) modeling the force-length relation (as in Heidlauf2013) should be multiplied. Set to false if this relation is already considered in the CellML model.
          "lambdaDotScalingFactor":       1.0,                       # scaling factor for the output of the lambda dot slot, i.e. the contraction velocity. Use this to scale the unit-less quantity to, e.g., micrometers per millisecond for the subcellular model.
          "slotNames":                    ["lambda", "ldot", "gamma", "T", "m_ux", "m_uy", "m_uz"],   # slot names of the data connector slots: lambda, lambdaDot, gamma, traction
          "dynamic":                      False,                     # if the dynamic formulation with velocity or the quasi-static formulation is computed
          
          "OutputWriter" : [
            {"format": "Paraview", "outputInterval": int(1./variables.dt_3D*variables.output_timestep), "filename": "out/" + variables.scenario_name + "/mechanics", "binary": True, "fixedFormat": False, "onlyNodalValues":True, "combineFiles":True, "fileNumbering": "incremental"},
          ],
          
          # the mechanics solver for the dynamic problem (if "dynamic": True)
          "DynamicHyperelasticitySolver": {
            "timeStepWidth":              variables.dt_3D,           # time step width 
            "durationLogKey":             "nonlinear",               # key to find duration of this solver in the log file
            "timeStepOutputInterval":     1,                         # how often the current time step should be printed to console
            
            "materialParameters":         variables.material_parameters,  # material parameters of the Mooney-Rivlin material
            "density":                    variables.rho,             # density of the material
            "displacementsScalingFactor": 1.0,                       # scaling factor for displacements, only set to sth. other than 1 only to increase visual appearance for very small displacements
            "residualNormLogFilename":    "log_residual_norm.txt",   # log file where residual norm values of the nonlinear solver will be written
            "useAnalyticJacobian":        True,                      # whether to use the analytically computed jacobian matrix in the nonlinear solver (fast)
            "useNumericJacobian":         False,                     # whether to use the numerically computed jacobian matrix in the nonlinear solver (slow), only works with non-nested matrices, if both numeric and analytic are enable, it uses the analytic for the preconditioner and the numeric as normal jacobian
              
            "dumpDenseMatlabVariables":   False,                     # whether to have extra output of matlab vectors, x,r, jacobian matrix (very slow)
            # if useAnalyticJacobian,useNumericJacobian and dumpDenseMatlabVariables all all three true, the analytic and numeric jacobian matrices will get compared to see if there are programming errors for the analytic jacobian
            
            "loadFactors":                [],                        # no load factors, solve problem directly
            "loadFactorGiveUpThreshold":  4e-2,                      # a threshold for the load factor, when to abort the solve of the current time step. The load factors are adjusted automatically if the nonlinear solver diverged. If the progression between two subsequent load factors gets smaller than this value, the solution is aborted.
            "scaleInitialGuess":          False,                     # when load stepping is used, scale initial guess between load steps a and b by sqrt(a*b)/a. This potentially reduces the number of iterations per load step (but not always).
            "nNonlinearSolveCalls":       1,                         # how often the nonlinear solve should be repeated
        
            # mesh
            "inputMeshIsGlobal":          True,                      # the boundary conditions for the mesh are given globally
            "meshName":                   "3Dmesh_quadratic",        # name of the 3D mesh, it is defined under "Meshes" at the beginning of this config
            "fiberDirection":             None,                      # if fiberMeshNames is empty, directly set the constant fiber direction, in element coordinate system
            "fiberDirectionInElement":    [0,0,1],                   # if fiberMeshNames is empty, directly set the constant fiber direction, in element coordinate system
          
            # solver
            "solverName":                 "mechanicsSolver",         # name of the nonlinear solver configuration, it is defined under "Solvers" at the beginning of this config
            
            # boundary and initial conditions
            "dirichletBoundaryConditions": variables.elasticity_dirichlet_bc,   # the initial Dirichlet boundary conditions that define values for displacements u and velocity v
            "neumannBoundaryConditions":   variables.elasticity_neumann_bc,     # Neumann boundary conditions that define traction forces on surfaces of elements
            "updateDirichletBoundaryConditionsFunction": None,                  # function that updates the dirichlet BCs while the simulation is running
            "updateDirichletBoundaryConditionsFunctionCallInterval": 1,         # every which step the update function should be called, 1 means every time step
            "divideNeumannBoundaryConditionValuesByTotalArea": True,            # if the given Neumann boundary condition values under "neumannBoundaryConditions" are total forces instead of surface loads and therefore should be scaled by the surface area of all elements where Neumann BC are applied
            
            "initialValuesDisplacements":  [[0.0,0.0,0.0] for _ in range(mx*my*mz)],     # the initial values for the displacements, vector of values for every node [[node1-x,y,z], [node2-x,y,z], ...]
            "initialValuesVelocities":     [[0.0,0.0,0.0] for _ in range(mx*my*mz)],     # the initial values for the velocities, vector of values for every node [[node1-x,y,z], [node2-x,y,z], ...]
            "constantBodyForce":           variables.constant_body_force,                 # a constant force that acts on the whole body, e.g. for gravity
            
            "dirichletOutputFilename":     None,                                    # filename for a vtp file that contains the Dirichlet boundary condition nodes and their values, set to None to disable
            
            # define which file formats should be written
            # 1. main output writer that writes output files using the quadratic elements function space. Writes displacements, velocities and PK2 stresses.
            "OutputWriter" : [
              
              # Paraview files
              #{"format": "Paraview", "outputInterval": 1, "filename": "out/u", "binary": False, "fixedFormat": False, "onlyNodalValues":True, "combineFiles":True},
              
              # Python callback function "postprocess"
              #{"format": "PythonCallback", "outputInterval": 1, "callback": postprocess, "onlyNodalValues":True, "filename": ""},
            ],
            # 2. additional output writer that writes also the hydrostatic pressure
            "pressure": {   # output files for pressure function space (linear elements), contains pressure values, as well as displacements and velocities
              "OutputWriter" : [
                #{"format": "Paraview", "outputInterval": 1, "filename": "out/p", "binary": False, "fixedFormat": False, "onlyNodalValues":True, "combineFiles":True},
              ]
            },
            # 3. additional output writer that writes virtual work terms
            "dynamic": {    # output of the dynamic solver, has additional virtual work values 
              "OutputWriter" : [   # output files for displacements function space (quadratic elements)
                #{"format": "Paraview", "outputInterval": int(output_interval/dt), "filename": "out/dynamic", "binary": False, "fixedFormat": False, "onlyNodalValues":True, "combineFiles":True},
                #{"format": "Paraview", "outputInterval": 1, "filename": "out/dynamic", "binary": False, "fixedFormat": False, "onlyNodalValues":True, "combineFiles":True},
              ],
            }
          },
          
          # the mechanics solver for the static problem (if "dynamic": False)
          "HyperelasticitySolver": {
            "durationLogKey":             "hyperelasticity",         # key to find duration of this solver in the log file
            "timeStepOutputInterval":     1,                         # how often the current time step should be printed to console
            
            "materialParameters":         variables.material_parameters,  # material parameters of the Mooney-Rivlin material
            "density":                    variables.rho,             # density of the material
            "displacementsScalingFactor": 1.0,                       # scaling factor for displacements, only set to sth. other than 1 only to increase visual appearance for very small displacements
            "residualNormLogFilename":    "log_residual_norm.txt",   # log file where residual norm values of the nonlinear solver will be written
            "useAnalyticJacobian":        True,                      # whether to use the analytically computed jacobian matrix in the nonlinear solver (fast)
            "useNumericJacobian":         False,                     # whether to use the numerically computed jacobian matrix in the nonlinear solver (slow), only works with non-nested matrices, if both numeric and analytic are enable, it uses the analytic for the preconditioner and the numeric as normal jacobian
              
            "dumpDenseMatlabVariables":   False,                     # whether to have extra output of matlab vectors, x,r, jacobian matrix (very slow)
            # if useAnalyticJacobian,useNumericJacobian and dumpDenseMatlabVariables all all three true, the analytic and numeric jacobian matrices will get compared to see if there are programming errors for the analytic jacobian
            
            #"loadFactors":  [0.1, 0.2, 0.35, 0.5, 1.0],             # load factors for every timestep
            "loadFactors":                [],                        # no load factors, solve problem directly
            #"loadFactorGiveUpThreshold":  4e-2,                      # a threshold for the load factor, when to abort the solve of the current time step. The load factors are adjusted automatically if the nonlinear solver diverged. If the progression between two subsequent load factors gets smaller than this value, the solution is aborted.
            "loadFactorGiveUpThreshold":  0.6,                       # a threshold for the load factor, when to abort the solve of the current time step. The load factors are adjusted automatically if the nonlinear solver diverged. If the progression between two subsequent load factors gets smaller than this value, the solution is aborted.
            "scaleInitialGuess":          False,                     # when load stepping is used, scale initial guess between load steps a and b by sqrt(a*b)/a. This potentially reduces the number of iterations per load step (but not always).
            "nNonlinearSolveCalls": 1,                               # how often the nonlinear solve should be repeated
        
            # mesh
            "inputMeshIsGlobal":          True,                      # the boundary conditions for the mesh are given globally
            "meshName":                   "3Dmesh_quadratic",        # name of the 3D mesh, it is defined under "Meshes" at the beginning of this config
            "fiberMeshNames":             [],                        # fiber meshes that will be used to determine the fiber direction, for multidomain there are no fibers so this would be empty list
            "fiberDirection":             None,                      # if fiberMeshNames is empty, directly set the constant fiber direction, in element coordinate system
            "fiberDirectionInElement":    [0,0,1],                   # if fiberMeshNames is empty, directly set the constant fiber direction, in element coordinate system
          
            # solver
            "solverName":                 "mechanicsSolver",         # name of the nonlinear solver configuration, it is defined under "Solvers" at the beginning of this config
            
            # boundary and initial conditions
            "dirichletBoundaryConditions": variables.main_elasticity_dirichlet_bc,   # the initial Dirichlet boundary conditions that define values for displacements u and velocity v
            "neumannBoundaryConditions":   variables.main_elasticity_neumann_bc,     # Neumann boundary conditions that define traction forces on surfaces of elements
            "updateDirichletBoundaryConditionsFunction": None,                  # function that updates the dirichlet BCs while the simulation is running
            "updateDirichletBoundaryConditionsFunctionCallInterval": 1,         # every which step the update function should be called, 1 means every time step
            "divideNeumannBoundaryConditionValuesByTotalArea": True,            # if the given Neumann boundary condition values under "neumannBoundaryConditions" are total forces instead of surface loads and therefore should be scaled by the surface area of all elements where Neumann BC are applied
            
            "initialValuesDisplacements":  [[0.0,0.0,0.0] for _ in range(mx*my*mz)],     # the initial values for the displacements, vector of values for every node [[node1-x,y,z], [node2-x,y,z], ...]
            "constantBodyForce":           variables.constant_body_force,                 # a constant force that acts on the whole body, e.g. for gravity
            
            "dirichletOutputFilename":     "out/"+variables.scenario_name+"/dirichlet_boundary_conditions",    # filename for a vtp file that contains the Dirichlet boundary condition nodes and their values, set to None to disable
            
            # define which file formats should be written
            # 1. main output writer that writes output files using the quadratic elements function space. Writes displacements, velocities and PK2 stresses.
            "OutputWriter" : [
              
              # Paraview files
              #{"format": "Paraview", "outputInterval": 1, "filename": "out/u", "binary": False, "fixedFormat": False, "onlyNodalValues":True, "combineFiles":True},
              
              # Python callback function "postprocess"
              #{"format": "PythonCallback", "outputInterval": 1, "callback": postprocess, "onlyNodalValues":True, "filename": ""},
            ],
            # 2. additional output writer that writes also the hydrostatic pressure
            "pressure": {   # output files for pressure function space (linear elements), contains pressure values, as well as displacements and velocities
              "OutputWriter" : [
                #{"format": "Paraview", "outputInterval": 1, "filename": "out/p", "binary": False, "fixedFormat": False, "onlyNodalValues":True, "combineFiles":True},
              ]
            },
            # 3. additional output writer that writes virtual work terms
            "dynamic": {    # output of the dynamic solver, has additional virtual work values 
              "OutputWriter" : [   # output files for displacements function space (quadratic elements)
                #{"format": "Paraview", "outputInterval": int(output_interval/dt), "filename": "out/dynamic", "binary": False, "fixedFormat": False, "onlyNodalValues":True, "combineFiles":True},
                #{"format": "Paraview", "outputInterval": 1, "filename": "out/dynamic", "binary": False, "fixedFormat": False, "onlyNodalValues":True, "combineFiles":True},
              ],
            },
            # output writer for debugging, outputs files after each load increment, the geometry is not changed but u and v are written
            "LoadIncrements": {   
              "OutputWriter" : [
                #{"format": "Paraview", "outputInterval": 1, "filename": "out/load_increment_", "binary": False, "fixedFormat": False, "onlyNodalValues":True, "combineFiles":True},
              ]
            },
          }
        }
      }
    }
  }
}

# stop timer and calculate how long parsing lasted
if rank_no == 0:
  t_stop_script = timeit.default_timer()
  print("Python config parsed in {:.1f}s.".format(t_stop_script - t_start_script))
