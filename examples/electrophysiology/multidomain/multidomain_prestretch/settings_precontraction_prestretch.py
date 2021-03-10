# biceps, contraction by constant activation
#

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
    print("Warning: There is no variables file, e.g:\n ./multidomain_prestretch ../settings_prestretch.py shortened.py\n")
  exit(0)
  
# -------------------------------------------------------- begin parameters ---------------------------------------------------------
# -------------------------------------------------------- end parameters ------------------------------------------------------------
parser = argparse.ArgumentParser(description='multidomain_contraction')
parser.add_argument('--scenario_name',                       help='The name to identify this run in the log.',   default=variables.scenario_name)
parser.add_argument('--n_subdomains', nargs=3,               help='Number of subdomains in x,y,z direction.',    type=int)
parser.add_argument('--n_subdomains_x', '-x',                help='Number of subdomains in x direction.',        type=int, default=variables.n_subdomains_x)
parser.add_argument('--n_subdomains_y', '-y',                help='Number of subdomains in y direction.',        type=int, default=variables.n_subdomains_y)
parser.add_argument('--n_subdomains_z', '-z',                help='Number of subdomains in z direction.',        type=int, default=variables.n_subdomains_z)
parser.add_argument('--paraview_output',                     help='Enable the paraview output writer.',          default=variables.paraview_output, action='store_true')
parser.add_argument('--adios_output',                        help='Enable the MegaMol/ADIOS output writer.',          default=variables.adios_output, action='store_true')
parser.add_argument('--fiber_file',                          help='The filename of the file that contains the fiber data.', default=variables.fiber_file)

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
  print("dt_0D:           {:0.0e}    multidomain solver:         {}, lumped mass matrix: {}".format(variables.dt_0D, variables.multidomain_solver_type, variables.use_lumped_mass_matrix))
  print("dt_multidomain:  {:0.0e}    multidomain preconditioner: {}, symmetric precond.: {}".format(variables.dt_multidomain, variables.multidomain_preconditioner_type, variables.use_symmetric_preconditioner_matrix))
  print("dt_splitting:    {:0.0e}    theta: {}, solver tolerances, abs: {}, rel: {}".format(variables.dt_splitting, variables.theta, variables.multidomain_absolute_tolerance, variables.multidomain_relative_tolerance))
  print("dt_elasticity:   {:0.0e}    elasticity solver: {}, preconditioner: {}".format(variables.dt_elasticity, variables.elasticity_solver_type, variables.elasticity_preconditioner_type))
  print("fiber_file:              {}".format(variables.fiber_file))
  print("fat_mesh_file:           {}".format(variables.fat_mesh_file))
  print("cellml_file:             {}".format(variables.cellml_file))
  print("firing_times_file:       {}".format(variables.firing_times_file))
  print("********************************************************************************")
  
  # start timer to measure duration of parsing of this script  
  t_start_script = timeit.default_timer()
    
# initialize all helper variables
from helper import *
  
config = {
  "scenarioName":          variables.scenario_name,
  "logFormat":                      "csv",                                # "csv" or "json", format of the lines in the log file, csv gives smaller files
  "solverStructureDiagramFile":     "out/"+variables.scenario_name+"/precontraction_solver_structure.txt",               # output file of a diagram that shows data connection between solvers
  "mappingsBetweenMeshesLogFile":   "out/"+variables.scenario_name+"/precontraction_mappings_between_meshes_log.txt",    # log file for mappings 
  "meta": {                                                               # additional fields that will appear in the log
    "partitioning":         [variables.n_subdomains_x, variables.n_subdomains_y, variables.n_subdomains_z]
  },
  "connectedSlots": [
    ("stress", "g_mu"),
    ("g_tot", "g_in")
  ],
  "Meshes":                variables.meshes,
  "MappingsBetweenMeshes": {
    "3Dmesh": [
       {"name": "3Dmesh_elasticity_quadratic",                                 "xiTolerance": 0.5, "enableWarnings": True, "compositeUseOnlyInitializedMappings": True, "fixUnmappedDofs": True, "defaultValue": 0},
    ],
    "3DFatMesh":  [
       {"name": "3DFatMesh_elasticity_quadratic",                              "xiTolerance": 0.5, "enableWarnings": True, "compositeUseOnlyInitializedMappings": True, "fixUnmappedDofs": True, "defaultValue": 0},
    ],
  },
  "Solvers": {
    "potentialFlowSolver": {
      "relativeTolerance":  1e-10,
      "absoluteTolerance":  1e-10,         # 1e-10 absolute tolerance of the residual          
      "maxIterations":      1e5,
      "solverType":         variables.potential_flow_solver_type,
      "preconditionerType": variables.potential_flow_preconditioner_type,
      "dumpFormat":         "default",
      "dumpFilename":       "",
    },
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
    }
  },
  "MultipleCoupling": {
    "timeStepWidth":          1,
    "logTimeStepWidthAsKey":  "dt_artifical_contraction",
    "durationLogKey":         "duration_artifical_contraction",
    "timeStepOutputInterval": 1,
    "endTime":                1,
    "connectedSlotsTerm1To2": None,      # connect lambda to slot 0 and gamma to slot 2
    "connectedSlotsTerm2To1": None,       # transfer nothing back
    
    # prescribed activation
    "Term1": {
      "PrescribedValues": {
        "meshName":               ["3Dmesh", "3DFatMesh"],      # reference to the multidomain mesh
        "numberTimeSteps":        1,             # number of timesteps to call the callback functions subsequently, this is usually 1 for prescribed values, because it is enough to set the reaction term only once per time step
        "timeStepOutputInterval": 20,            # if the time step should be written to console, a value > 10 produces no output
        "slotNames":              ["lambda", "gamma"],
        
        # a list of field variables that will get values assigned in every timestep, by the provided callback function
        "fieldVariables1": [
          {"name": "lambda", "callback": variables.set_lambda_values},
          {"name": "gamma", "callback": variables.set_gamma_values},
        ],
        "fieldVariables2":     [],
        "additionalArgument":  None,         # a custom argument to the fieldVariables callback functions, this will be passed on as the last argument
        
        "OutputWriter" : [
          {"format": "Paraview", "outputInterval": 1, "filename": "out/" + variables.scenario_name + "/prescribed_stress", "binary": True, "fixedFormat": False, "combineFiles": True, "fileNumbering": "incremental"}
        ]
      }
    },
    "Term2": {  
      "MuscleContractionSolver": {
        "numberTimeSteps":              1,                         # only use 1 timestep per interval
        "timeStepOutputInterval":       100,                       # do not output time steps
        "Pmax":                         variables.pmax,            # maximum PK2 active stress
        "enableForceLengthRelation":    False,                     # if the factor f_l(λ_f) modeling the force-length relation (as in Heidlauf2013) should be multiplied. Set to false if this relation is already considered in the CellML model.
        "lambdaDotScalingFactor":       1.0,                       # scaling factor for the output of the lambda dot slot, i.e. the contraction velocity. Use this to scale the unit-less quantity to, e.g., micrometers per millisecond for the subcellular model.
        "slotNames":                    ["lambda", "ldot", "gamma", "T"],
        "OutputWriter" : [
          {"format": "Paraview", "outputInterval": int(1./variables.dt_elasticity*variables.output_timestep_elasticity), "filename": "out/"+variables.scenario_name+"/precontraction_mechanics_3D", "binary": True, "fixedFormat": False, "onlyNodalValues":True, "combineFiles": True, "fileNumbering": "incremental"},
        ],
        "mapGeometryToMeshes":          ["3Dmesh_elasticity_quadratic+3DFatMesh_elasticity_quadratic"],    # the mesh names of the meshes that will get the geometry transferred
        "reverseMappingOrder":          True,                      # if the mapping target->own mesh should be used instead of own->target mesh. This gives better results in some cases.
        "dynamic":                      False,                     # if the dynamic solid mechanics solver should be used, else it computes the quasi-static problem
        
        # the actual solid mechanics solver, this is either "DynamicHyperelasticitySolver" or "HyperelasticitySolver", depending on the value of "dynamic"
                                 
        "HyperelasticitySolver": {
          "durationLogKey":             "duration_precontraction",               # key to find duration of this solver in the log file
          
          "materialParameters":         variables.material_parameters,  # material parameters of the Mooney-Rivlin material
          "displacementsScalingFactor": 1.0,                       # scaling factor for displacements, only set to sth. other than 1 only to increase visual appearance for very small displacements
          "residualNormLogFilename":    "out/"+variables.scenario_name+"/log_precontraction_residual_norm.txt",   # log file where residual norm values of the nonlinear solver will be written
          "useAnalyticJacobian":        True,                      # whether to use the analytically computed jacobian matrix in the nonlinear solver (fast)
          "useNumericJacobian":         False,                     # whether to use the numerically computed jacobian matrix in the nonlinear solver (slow), only works with non-nested matrices, if both numeric and analytic are enable, it uses the analytic for the preconditioner and the numeric as normal jacobian
            
          "dumpDenseMatlabVariables":   False,                     # whether to have extra output of matlab vectors, x,r, jacobian matrix (very slow)
          # if useAnalyticJacobian,useNumericJacobian and dumpDenseMatlabVariables all all three true, the analytic and numeric jacobian matrices will get compared to see if there are programming errors for the analytic jacobian
          
          # mesh
          "inputMeshIsGlobal":          True,                     # boundary conditions and initial values are given as global numbers (every process has all information)
          "meshName":                   ["3Dmesh_elasticity_quadratic_precontraction", "3DFatMesh_elasticity_quadratic_precontraction"],       # name of the 3D mesh, it is defined under "Meshes" at the beginning of this config
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
          
          "initialValuesDisplacements":  [[0.0,0.0,0.0] for _ in range(variables.n_points_global_composite_mesh)],     # the initial values for the displacements, vector of values for every node [[node1-x,y,z], [node2-x,y,z], ...]
          "initialValuesVelocities":     [[0.0,0.0,0.0] for _ in range(variables.n_points_global_composite_mesh)],     # the initial values for the velocities, vector of values for every node [[node1-x,y,z], [node2-x,y,z], ...]
          "extrapolateInitialGuess":     True,                                # if the initial values for the dynamic nonlinear problem should be computed by extrapolating the previous displacements and velocities
          "constantBodyForce":           variables.precontraction_constant_body_force,       # a constant force that acts on the whole body, e.g. for gravity
          
          "dirichletOutputFilename":    "out/"+variables.scenario_name+"/precontraction_dirichlet_boundary_conditions",     # output filename for the dirichlet boundary conditions, set to "" to have no output
          
          # define which file formats should be written
          # 1. main output writer that writes output files using the quadratic elements function space. Writes displacements, velocities and PK2 stresses.
          "OutputWriter" : [
            
            # Paraview files
            {"format": "Paraview", "outputInterval": 1, "filename": "out/"+variables.scenario_name+"/precontraction_u", "binary": True, "fixedFormat": False, "onlyNodalValues":True, "combineFiles":True, "fileNumbering": "incremental"},
          ],
          # 2. additional output writer that writes also the hydrostatic pressure
          "pressure": {   # output files for pressure function space (linear elements), contains pressure values, as well as displacements and velocities
            "OutputWriter" : [
              {"format": "Paraview", "outputInterval": 1, "filename": "out/"+variables.scenario_name+"/precontraction_p", "binary": True, "fixedFormat": False, "onlyNodalValues":True, "combineFiles":True, "fileNumbering": "incremental"},
            ]
          },
          # 4. output writer for debugging, outputs files after each load increment, the geometry is not changed but u and v are written
          "LoadIncrements": {   
            "OutputWriter" : [
              {"format": "Paraview", "outputInterval": 1, "filename": "out/"+variables.scenario_name+"/precontraction_load_increments", "binary": True, "fixedFormat": False, "onlyNodalValues":True, "combineFiles":True, "fileNumbering": "incremental"},
            ]
          }
        }
      }
    },
    "Term3": {
      "HyperelasticitySolver": {
        "durationLogKey":             "duration_prestretch",               # key to find duration of this solver in the log file
        
        "materialParameters":         variables.material_parameters,  # material parameters of the Mooney-Rivlin material
        "displacementsScalingFactor": 1.0,                       # scaling factor for displacements, only set to sth. other than 1 only to increase visual appearance for very small displacements
        "residualNormLogFilename":    "out/"+variables.scenario_name+"/log_prestretch_residual_norm.txt",   # log file where residual norm values of the nonlinear solver will be written
        "useAnalyticJacobian":        True,                      # whether to use the analytically computed jacobian matrix in the nonlinear solver (fast)
        "useNumericJacobian":         False,                     # whether to use the numerically computed jacobian matrix in the nonlinear solver (slow), only works with non-nested matrices, if both numeric and analytic are enable, it uses the analytic for the preconditioner and the numeric as normal jacobian
          
        "dumpDenseMatlabVariables":   False,                     # whether to have extra output of matlab vectors, x,r, jacobian matrix (very slow)
        # if useAnalyticJacobian,useNumericJacobian and dumpDenseMatlabVariables all all three true, the analytic and numeric jacobian matrices will get compared to see if there are programming errors for the analytic jacobian
        
        # mesh
        "inputMeshIsGlobal":          True,                     # boundary conditions and initial values are given as global numbers (every process has all information)
        "meshName":                   ["3Dmesh_elasticity_quadratic", "3DFatMesh_elasticity_quadratic"],       # name of the 3D mesh, it is defined under "Meshes" at the beginning of this config
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
        
        "initialValuesDisplacements":  [[0.0,0.0,0.0] for _ in range(variables.n_points_global_composite_mesh)],     # the initial values for the displacements, vector of values for every node [[node1-x,y,z], [node2-x,y,z], ...]
        "initialValuesVelocities":     [[0.0,0.0,0.0] for _ in range(variables.n_points_global_composite_mesh)],     # the initial values for the velocities, vector of values for every node [[node1-x,y,z], [node2-x,y,z], ...]
        "extrapolateInitialGuess":     True,                                # if the initial values for the dynamic nonlinear problem should be computed by extrapolating the previous displacements and velocities
        "constantBodyForce":           variables.prestretch_constant_body_force,       # a constant force that acts on the whole body, e.g. for gravity
        
        "dirichletOutputFilename":    "out/"+variables.scenario_name+"/prestretch_dirichlet_boundary_conditions",     # output filename for the dirichlet boundary conditions, set to "" to have no output
        
        # define which file formats should be written
        # 1. main output writer that writes output files using the quadratic elements function space. Writes displacements, velocities and PK2 stresses.
        "OutputWriter" : [
          
          # Paraview files
          {"format": "Paraview", "outputInterval": 1, "filename": "out/"+variables.scenario_name+"/prestretch_u", "binary": True, "fixedFormat": False, "onlyNodalValues":True, "combineFiles":True, "fileNumbering": "incremental"},
        ],
        # 2. additional output writer that writes also the hydrostatic pressure
        "pressure": {   # output files for pressure function space (linear elements), contains pressure values, as well as displacements and velocities
          "OutputWriter" : [
            {"format": "Paraview", "outputInterval": 1, "filename": "out/"+variables.scenario_name+"/prestretch_p", "binary": True, "fixedFormat": False, "onlyNodalValues":True, "combineFiles":True, "fileNumbering": "incremental"},
          ]
        },
        # 4. output writer for debugging, outputs files after each load increment, the geometry is not changed but u and v are written
        "LoadIncrements": {   
          "OutputWriter" : [
            {"format": "Paraview", "outputInterval": 1, "filename": "out/"+variables.scenario_name+"/prestretch_load_increments", "binary": True, "fixedFormat": False, "onlyNodalValues":True, "combineFiles":True, "fileNumbering": "incremental"},
          ]
        },
      }
    }
  }
}

# stop timer and calculate how long parsing lasted
if rank_no == 0:
  t_stop_script = timeit.default_timer()
  print("Python config parsed in {:.1f}s.".format(t_stop_script - t_start_script))

