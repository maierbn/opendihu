# sensory organs, interneurons, motoneurons, multidomain mockup (no electrophyisiology), mechanics
#
# [muscle spindle model] -> [motor neuron model] <- mapping by callback <- [interneuron model] <- mapping by callback <- [golgi tendon organ model]
#                            +-> electro-mechanics model

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
    print("Warning: There is no variables file, e.g:\n ./only_neurons ../settings_only_neurons.py spindles.py\n")
  exit(0)
  
# -------------------------------------------------------- begin parameters ---------------------------------------------------------
# -------------------------------------------------------- end parameters ------------------------------------------------------------
parser = argparse.ArgumentParser(description='only_neurons')
parser.add_argument('--scenario_name',                       help='The name to identify this run in the log.',   default=variables.scenario_name)
parser.add_argument('--n_subdomains', nargs=3,               help='Number of subdomains in x,y,z direction.',    type=int)
parser.add_argument('--n_subdomains_x', '-x',                help='Number of subdomains in x direction.',        type=int, default=variables.n_subdomains_x)
parser.add_argument('--n_subdomains_y', '-y',                help='Number of subdomains in y direction.',        type=int, default=variables.n_subdomains_y)
parser.add_argument('--n_subdomains_z', '-z',                help='Number of subdomains in z direction.',        type=int, default=variables.n_subdomains_z)
parser.add_argument('--potential_flow_solver_type',          help='The solver for the potential flow (non-spd matrix).', default=variables.potential_flow_solver_type, choices=["gmres","cg","lu","gamg","richardson","chebyshev","cholesky","jacobi","sor","preonly"])
parser.add_argument('--potential_flow_preconditioner_type',  help='The preconditioner for the potential flow.',  default=variables.potential_flow_preconditioner_type, choices=["jacobi","sor","lu","ilu","gamg","none"])
parser.add_argument('--paraview_output',                     help='Enable the paraview output writer.',          default=variables.paraview_output, action='store_true')
parser.add_argument('--adios_output',                        help='Enable the MegaMol/ADIOS output writer.',          default=variables.adios_output, action='store_true')
parser.add_argument('--fiber_file',                          help='The filename of the file that contains the fiber data.', default=variables.fiber_file)
parser.add_argument('--firing_times_file',                   help='The filename of the file that contains the cellml model.', default=variables.firing_times_file)
parser.add_argument('--end_time', '--tend', '-t',            help='The end simulation time.',                    type=float, default=variables.end_time)
parser.add_argument('--output_timestep_multidomain',         help='The timestep for writing outputs.',           type=float, default=variables.output_timestep_multidomain)
parser.add_argument('--multidomain_solver_type',             help='Solver for multidomain system.',              default=variables.multidomain_solver_type)
parser.add_argument('--multidomain_preconditioner_type',     help='Preconditioner for multidomain system.',      default=variables.multidomain_preconditioner_type)
parser.add_argument('--multidomain_relative_tolerance',      help='relative residual tol. for multidomain solver.',           type=float, default=variables.multidomain_relative_tolerance)
parser.add_argument('--multidomain_absolute_tolerance',      help='absolute residual tol. for multidomain solver.',           type=float, default=variables.multidomain_absolute_tolerance)
parser.add_argument('--theta',                               help='parameter of Crank-Nicholson in multidomain system.',      type=float, default=variables.theta)
parser.add_argument('--use_lumped_mass_matrix',              help='If the formulation with lumbed mass matrix should be used.',           default=variables.use_lumped_mass_matrix, action='store_true')
parser.add_argument('--use_symmetric_preconditioner_matrix', help='If the preconditioner matrix should be symmetric (only diagnoal part of system matrix).',  default=variables.use_symmetric_preconditioner_matrix, action='store_true')
parser.add_argument('--initial_guess_nonzero',               help='If the initial guess to the linear solver should be the last solution.',  default=variables.initial_guess_nonzero, action='store_true')

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
#variables.meshes = {}

# callback function that receives the whole result values and produces plots while the simulation is running
def handle_result(n_instances, time_step_no, current_time, states, algebraics, name_information, additional_argument):
    
  print("handle result of muscle spindle")
  #print(name_information)
    
  # asign some states to variables
  Vm = states[name_information["stateNames"].index("membrane/V")]
  print("Vm: {}".format(Vm))
    

if False:
  # settings for the multidomain solver
  multidomain_solver = {
    "timeStepWidth":                    variables.dt_multidomain,             # time step width of the subcellular problem
    "endTime":                          variables.end_time,                   # end time, this is not relevant because it will be overridden by the splitting scheme
    "timeStepOutputInterval":           1,                                  # how often the output timestep should be printed
    "durationLogKey":                   "duration_multidomain",               # key for duration in log.csv file
    
    # material parameters for the compartments
    "nCompartments":                    variables.n_compartments,             # number of compartments
    "compartmentRelativeFactors":       variables.relative_factors.tolist(),  # list of lists of the factors for every dof, because "inputIsGlobal": True, this contains the global dofs
    "inputIsGlobal":                    True,                                 # if values and dofs correspond to the global numbering
    "am":                               [variables.get_am(mu_no) for mu_no in range(variables.n_compartments)],   # Am parameter for every motor unit (ration of surface to volume of fibers)
    "cm":                               [variables.get_cm(mu_no) for mu_no in range(variables.n_compartments)],   # Cm parameter for every motor unit (capacitance of the cellular membrane)
    
    # solver options
    "solverName":                       "multidomainLinearSolver",            # reference to the solver used for the global linear system of the multidomain eq.
    "theta":                            variables.theta,                      # weighting factor of implicit term in Crank-Nicolson scheme, 0.5 gives the classic, 2nd-order Crank-Nicolson scheme, 1.0 gives implicit euler
    "useLumpedMassMatrix":              variables.use_lumped_mass_matrix,     # which formulation to use, the formulation with lumped mass matrix (True) is more stable but approximative, the other formulation (False) is exact but needs more iterations
    "useSymmetricPreconditionerMatrix": variables.use_symmetric_preconditioner_matrix,    # if the diagonal blocks of the system matrix should be used as preconditioner matrix
    "initialGuessNonzero":              variables.initial_guess_nonzero,      # if the initial guess for the 3D system should be set as the solution of the previous timestep, this only makes sense for iterative solvers
    "enableFatComputation":             True,                                 # disabling the computation of the fat layer is only for debugging and speeds up computation. If set to False, the respective matrix is set to the identity
    "showLinearSolverOutput":           variables.show_linear_solver_output,  # if convergence information of the linear solver in every timestep should be printed, this is a lot of output for fast computations
    "updateSystemMatrixEveryTimestep":  False,                                # if this multidomain solver will update the system matrix in every first timestep, us this only if the geometry changed, e.g. by contraction
    
    "PotentialFlow": {
      "FiniteElementMethod" : {  
        "meshName":                     "3Dmesh",
        "solverName":                   "potentialFlowSolver",
        "prefactor":                    1.0,
        "dirichletBoundaryConditions":  variables.potential_flow_dirichlet_bc,
        "dirichletOutputFilename":      None,                                 # filename for a vtp file that contains the Dirichlet boundary condition nodes and their values, set to None to disable
        "neumannBoundaryConditions":    [],
        "inputMeshIsGlobal":            True,
      },
    },
    "Activation": {
      "FiniteElementMethod" : {  
        "meshName":                     "3Dmesh",
        "solverName":                   "multidomainLinearSolver",
        "prefactor":                    1.0,
        "inputMeshIsGlobal":            True,
        "dirichletBoundaryConditions":  {},
        "dirichletOutputFilename":      None,                                # filename for a vtp file that contains the Dirichlet boundary condition nodes and their values, set to None to disable
        "neumannBoundaryConditions":    [],
        "diffusionTensor": [[      # sigma_i           # fiber direction is (1,0,0)
          8.93, 0, 0,
          0, 0.0, 0,
          0, 0, 0.0
        ]], 
        "extracellularDiffusionTensor": [[      # sigma_e
          6.7, 0, 0,
          0, 6.7, 0,
          0, 0, 6.7,
        ]],
      },
    },
    "Fat": {
      "FiniteElementMethod" : {  
        "meshName":                     "3DFatMesh",
        "solverName":                   "multidomainLinearSolver",
        "prefactor":                    0.4,
        "inputMeshIsGlobal":            True,
        "dirichletBoundaryConditions":  {},
        "dirichletOutputFilename":      None,                                # filename for a vtp file that contains the Dirichlet boundary condition nodes and their values, set to None to disable
        "neumannBoundaryConditions":    [],
      },
    },
    
    "OutputWriter" : [
      {"format": "Paraview", "outputInterval": (int)(1./variables.dt_multidomain*variables.output_timestep_multidomain), "filename": "out/"+variables.scenario_name+"/multidomain", "binary": True, "fixedFormat": False, "combineFiles": True, "fileNumbering": "incremental"},
      #{"format": "ExFile", "filename": "out/fiber_"+str(i), "outputInterval": 1./dt_1D*output_timestep, "sphereSize": "0.02*0.02*0.02", "fileNumbering": "incremental"},
      #{"format": "PythonFile", "filename": "out/fiber_"+str(i), "outputInterval": int(1./dt_1D*output_timestep), "binary":True, "onlyNodalValues":True, "fileNumbering": "incremental"},
    ]
  }
  
# add neuron meshes
variables.meshes.update(
{
  "motoneuronMesh": {
    "nElements" :         variables.n_motoneurons-1 if n_ranks == 1 else variables.n_motoneurons*n_ranks,
    "physicalExtent":     0,
    "physicalOffset":     0,
    "logKey":             "motoneuron",
    "inputMeshIsGlobal":  True,
    "nRanks":             n_ranks
  },
  "interneuronMesh": {
    "nElements" :         variables.n_interneurons-1 if n_ranks == 1 else variables.n_interneurons*n_ranks,
    "physicalExtent":     0,
    "physicalOffset":     0,
    "logKey":             "interneuron",
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
  },
  "muscleSpindleAndInterneuronMesh": {
    "nElements" :         variables.n_muscle_spindles+variables.n_interneurons-1 if n_ranks == 1 else (variables.n_muscle_spindles+variables.n_interneurons)*n_ranks,
    "physicalExtent":     0,
    "physicalOffset":     0,
    "logKey":             "muscle_spindle_interneuron",
    "inputMeshIsGlobal":  True,
    "nRanks":             n_ranks
  },
  "golgiTendonOrganMesh": {
    "nElements" :         variables.n_golgi_tendon_organs-1 if n_ranks == 1 else variables.n_golgi_tendon_organs*n_ranks,
    "physicalExtent":     0,
    "physicalOffset":     0,
    "logKey":             "golgi_tendon_organ",
    "inputMeshIsGlobal":  True,
    "nRanks":             n_ranks
  },
})
  
config = {
  "scenarioName":          variables.scenario_name,
  "logFormat":                      "csv",                                # "csv" or "json", format of the lines in the log file, csv gives smaller files
  "solverStructureDiagramFile":     "solver_structure.txt",               # output file of a diagram that shows data connection between solvers
  "mappingsBetweenMeshesLogFile":   "mappings_between_meshes_log.txt",    # log file for mappings 
  "meta": {                                                               # additional fields that will appear in the log
    "partitioning":         [variables.n_subdomains_x, variables.n_subdomains_y, variables.n_subdomains_z]
  },
  "Meshes":                variables.meshes,
  "MappingsBetweenMeshes": {
    "3Dmesh": [
       {"name": "3Dmesh_elasticity_quadratic",                                 "xiTolerance": 1.5, "enableWarnings": True, "compositeUseOnlyInitializedMappings": True, "fixUnmappedDofs": True},
       {"name": "3Dmesh_elasticity_quadratic+3DFatMesh_elasticity_quadratic",  "xiTolerance": 1.5, "enableWarnings": False, "compositeUseOnlyInitializedMappings": True, "fixUnmappedDofs": True},    # mapping uses mappings of submeshes (i.e. 3Dmesh_elasticity_quadratic->3Dmesh)
    ],
    "3DFatMesh":  [
       {"name": "3DFatMesh_elasticity_quadratic",                              "xiTolerance": 1.5, "enableWarnings": True, "compositeUseOnlyInitializedMappings": True, "fixUnmappedDofs": True},
       {"name": "3Dmesh_elasticity_quadratic+3DFatMesh_elasticity_quadratic",  "xiTolerance": 1.5, "enableWarnings": False, "compositeUseOnlyInitializedMappings": True, "fixUnmappedDofs": True},    # mapping uses mappings of submeshes (i.e. 3Dmesh_elasticity_quadratic->3Dmesh)
    ]
  },
  "Solvers": {
    # the next two solvers are not used because of the mockup
    #"potentialFlowSolver": {
    #  "relativeTolerance":  1e-10,
    #  "absoluteTolerance":  1e-10,         # 1e-10 absolute tolerance of the residual          
    #  "maxIterations":      1e5,
    #  "solverType":         "gmres",
    #  "preconditionerType": "none",
    #  "dumpFormat":         "default",
    #  "dumpFilename":       "",
    #},
    #"multidomainLinearSolver": {
    #  "relativeTolerance":  variables.multidomain_relative_tolerance,
    #  "absoluteTolerance":  variables.multidomain_absolute_tolerance,     # absolute tolerance of the residual          
    #  "maxIterations":      1e4,
    #  "solverType":         variables.multidomain_solver_type,
    #  "preconditionerType": variables.multidomain_preconditioner_type,
    #  "hypreOptions":       "",                                           # additional options if a hypre preconditioner is selected
    #  "dumpFormat":         "matlab",
    #  "dumpFilename":       "",
    #},
    "mechanicsSolver": {   # solver for the dynamic mechanics problem
      "relativeTolerance":   variables.snes_relative_tolerance,           # 1e-10 relative tolerance of the linear solver
      "absoluteTolerance":   variables.snes_absolute_tolerance,           # 1e-10 absolute tolerance of the residual of the linear solver
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
  "connectedSlots": [
    ("mn_out", "mn"),
    ("in_g",   "in_in"),
    ("msin_s", "msin_i"),
    ("msin_i", "msin_m"),
    ("gt",     "gt_in"),
    ("ms",     "ms_in"),
  ],
  
  # total solver
  "Coupling": {
    "description":            "overall solver",    # description that will be shown in solver structure visualization
    "endTime":                variables.end_time,
    "timeStepWidth":          variables.dt_stimulation_check,
    "logTimeStepWidthAsKey":  "dt_neurons",
    "durationLogKey":         "duration_total",
    "timeStepOutputInterval": 100,
    "connectedSlotsTerm1To2": None,
    "connectedSlotsTerm2To1": None,

    # muscle spindles, golgi tendon organs and interneurons
    "Term1": {
      # muscle spindles, golgi tendon organs and interneurons = motoneuron input solver
      "Coupling": {
        "description":            "sensory organs and interneurons",    # description that will be shown in solver structure visualization
        "endTime":                variables.end_time,
        "timeStepWidth":          variables.dt_neuron_transfer,
        "logTimeStepWidthAsKey":  "dt_neurons",
        "durationLogKey":         "duration_total",
        "timeStepOutputInterval": 100,
        "connectedSlotsTerm1To2": None,
        "connectedSlotsTerm2To1": None,
        
        "Term1": {
          # mapping muscle spindles output -> motor neuron signals
          "MapDofs": {
            "description":                "muscle_spindles_to_motoneurons",   # description that will be shown in solver structure visualization
            "nAdditionalFieldVariables":  1,                                  # number of additional field variables that are defined by this object. They have 1 component, use the templated function space and mesh given by meshName.
            "additionalSlotNames":        ["msin_s"],
            "meshName":                   "muscleSpindleAndInterneuronMesh",  # the mesh on which the additional field variables will be defined
            "beforeComputation": None,
            "afterComputation": [                                         # transfer/mapping of dofs that will be performed after the computation of the nested solver
              {                                                 
                "fromConnectorSlot":                0,
                "toConnectorSlots":                 2,
                "fromSlotConnectorArrayIndex":      0,                    # which fiber/compartment
                "toSlotConnectorArrayIndex":        0,
                "mode":                             "callback",           # "copyLocal", "copyLocalIfPositive", "localSetIfAboveThreshold" or "communicate"
                "fromDofNosNumbering":              "local",
                "toDofNosNumbering":                "local",
                "dofsMapping":                      None,
                "inputDofs":                        list(range(variables.n_muscle_spindles)),   # [0,1,...,n_muscle_spindles], this is for mesh "muscleSpindleMesh"
                "outputDofs":                       [list(range(variables.n_muscle_spindles))],   # [0,1,...,n_muscle_spindles], this is for mesh "muscleSpindleAndInterneuronMesh"
                "callback":                         variables.callback_muscle_spindles_to_motoneurons,
                #"thresholdValue":                   20,                  # if mode is "localSetIfAboveThreshold", this is the threshold, if the value is above it, set the value `valueToSet`
                #"valueToSet":                       20,                  # if mode is "localSetIfAboveThreshold", this is the value to set the target dof to, if the source dof is above thresholdValue.
              }
            ],
                   
            # muscle spindle model solver
            "Heun" : {
              "description":                  "muscle spindle",
              "timeStepWidth":                variables.dt_muscle_spindles,
              "logTimeStepWidthAsKey":        "dt_muscle_spindles",
              "durationLogKey":               "duration_muscle_spindles",
              "initialValues":                [],
              "timeStepOutputInterval":       1e4,
              "inputMeshIsGlobal":            True,
              "dirichletBoundaryConditions":  {},
              "dirichletOutputFilename":      None,             # filename for a vtp file that contains the Dirichlet boundary condition nodes and their values, set to None to disable
              "checkForNanInf":               True,             # check if the solution vector contains nan or +/-inf values, if yes, an error is printed. This is a time-consuming check.
              "nAdditionalFieldVariables":    0,
              "additionalSlotNames":          [],
                  
              # cellml model of muscle spindle
              "CellML" : {
                "modelFilename":                          variables.muscle_spindle_cellml_file,           # input C++ source file or cellml XML file
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
                                     
                "handleResultFunction":                   None, #handle_result,              # callback function that gets all current values and can do something with them
                "handleResultCallInterval":               1,                          # interval in which handle_result will be called
                "handleResultFunctionAdditionalParameter": None,                      # additional last argument for handle_result
                
                "mappings":                               variables.muscle_spindle_mappings,              # mappings between parameters and algebraics/constants and between outputConnectorSlots and states, algebraics or parameters, they are defined in helper.py
                "parametersInitialValues":                variables.muscle_spindle_parameters_initial_values,  # initial values for the parameters, either once for all instances, or for all instances in array of struct ordering with nParameters_ parameters per instance: [inst0p0, inst0p1, ... inst0pn, inst1p0, inst1p1, ...]
                
                "meshName":                               "muscleSpindleMesh",                                       
                "stimulationLogFilename":                 "out/stimulation.log",

                # output writer for states, algebraics and parameters                
                "OutputWriter" : [
                  {"format": "PythonFile", "outputInterval": (int)(2./variables.dt_muscle_spindles*variables.output_timestep_neurons), "filename": "out/muscle_spindles", "binary": True, "fixedFormat": False, "combineFiles": True, "onlyNodalValues": True, "fileNumbering": "incremental"}
                ]
              }
            }
          }
        },
        "Term2": {
          # processed Golgi tendon sensor solver
          "Coupling": {
            "description":            "processed Golgi tendon sensors",      # description that will be shown in solver structure visualization
            "endTime":                variables.end_time,
            "timeStepWidth":          variables.dt_neuron_transfer,
            "logTimeStepWidthAsKey":  "dt_golgi_tendon_organs",
            "durationLogKey":         "duration_golgi_tendon",
            "timeStepOutputInterval": 100,
            "connectedSlotsTerm1To2": None,
            "connectedSlotsTerm2To1": None,
            
            "Term1": {
              # mapping Golgi tendon organs -> interneurons
              "MapDofs": {
                "description":                "golgi_tendon_o._to_interneurons",       # description that will be shown in solver structure visualization
                "nAdditionalFieldVariables":  1,                             # number of additional field variables that are defined by this object. They have 1 component, use the templated function space and mesh given by meshName.
                "additionalSlotNames":        ["in_g"],
                "meshName":                   "interneuronMesh",             # the mesh on which the additional field variables will be defined
                "beforeComputation": None,
                "afterComputation": [                                        # transfer/mapping of dofs that will be performed after the computation of the nested solver
                  {                                                 
                    "fromConnectorSlot":                0,
                    "toConnectorSlots":                 2,
                    "fromSlotConnectorArrayIndex":      0,                   # which fiber/compartment
                    "toSlotConnectorArrayIndex":        0,
                    "mode":                             "callback",          # "copyLocal", "copyLocalIfPositive", "localSetIfAboveThreshold" or "communicate"
                    "fromDofNosNumbering":              "local",
                    "toDofNosNumbering":                "local",
                    "dofsMapping":                      None,
                    "inputDofs":                        list(range(variables.n_golgi_tendon_organs)),   # [0,1,...,n_golgi_tendon_organs], this is for mesh "golgiTendonOrganMesh"
                    "outputDofs":                       [list(range(variables.n_interneurons))],          # [0,1,...,n_interneurons]         this is for mesh "interneuronMesh"
                    "callback":                         variables.callback_golgi_tendon_organs_to_interneurons,
                    #"thresholdValue":                   20,                 # if mode is "localSetIfAboveThreshold", this is the threshold, if the value is above it, set the value `valueToSet`
                    #"valueToSet":                       20,                 # if mode is "localSetIfAboveThreshold", this is the value to set the target dof to, if the source dof is above thresholdValue.
                  }
                ],
                # Golgi tendon organ model solver
                "Heun" : {
                  "description":                  "Golgi tendon organs",
                  "timeStepWidth":                variables.dt_golgi_tendon_organs,
                  "logTimeStepWidthAsKey":        "dt_golgi_tendon_organs",
                  "durationLogKey":               "duration_golgi_tendon_organ",
                  "initialValues":                [],
                  "timeStepOutputInterval":       1e4,
                  "inputMeshIsGlobal":            True,
                  "dirichletBoundaryConditions":  {},
                  "dirichletOutputFilename":      None,             # filename for a vtp file that contains the Dirichlet boundary condition nodes and their values, set to None to disable
                  "checkForNanInf":               True,             # check if the solution vector contains nan or +/-inf values, if yes, an error is printed. This is a time-consuming check.
                  "nAdditionalFieldVariables":    0,
                  "additionalSlotNames":          [],
                      
                  # cellml model of golgi tendon organs
                  "CellML" : {
                    "modelFilename":                          variables.golgi_tendon_organ_cellml_file,       # input C++ source file or cellml XML file
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
                    
                    "mappings":                               variables.golgi_tendon_organ_mappings,          # mappings between parameters and algebraics/constants and between connectorSlots and states, algebraics or parameters, they are defined in helper.py
                    "parametersInitialValues":                variables.golgi_tendon_organ_parameters_initial_values,    # # initial values for the parameters, either once for all instances, or for all instances in array of struct ordering with nParameters_ parameters per instance: [inst0p0, inst0p1, ... inst0pn, inst1p0, inst1p1, ...]
                    
                    "meshName":                               "golgiTendonOrganMesh",                         
                    "stimulationLogFilename":                 "out/stimulation.log",

                    # output writer for states, algebraics and parameters                
                    "OutputWriter" : [
                      {"format": "PythonFile", "outputInterval": (int)(2./variables.dt_golgi_tendon_organs*variables.output_timestep_neurons), "filename": "out/golgi_tendon_organs", "binary": True, "fixedFormat": False, "combineFiles": True, "onlyNodalValues": True, "fileNumbering": "incremental"}
                    ]
                  }
                }
              }
            },
            "Term2": {
              # mapping interneurons -> input for motor neurons
              "MapDofs": {
                "description":                "interneurons_to_motoneurons",  # description that will be shown in solver structure visualization
                "nAdditionalFieldVariables":  1,                              # number of additional field variables that are defined by this object. They have 1 component, use the templated function space and mesh given by meshName.
                "additionalSlotNames":        ["msin_i"],
                "meshName":                   "muscleSpindleAndInterneuronMesh", # the mesh on which the additional field variables will be defined
                "beforeComputation": None,                                    # transfer/mapping of dofs that will be performed before the computation of the nested solver
                "afterComputation": [                                         # transfer/mapping of dofs that will be performed after the computation of the nested solver
                  {                                                 
                    "fromConnectorSlot":          0,
                    "toConnectorSlots":           2,
                    "fromSlotConnectorArrayIndex":    0,                    # which fiber/compartment
                    "toSlotConnectorArrayIndex":      0,
                    "mode":                             "callback",           # "copyLocal", "copyLocalIfPositive", "localSetIfAboveThreshold" or "communicate"
                    "fromDofNosNumbering":              "local",
                    "toDofNosNumbering":                "local",
                    "dofsMapping":                      None,
                    "inputDofs":                        list(range(variables.n_interneurons)),   # [0,1,...,n_interneurons]
                    "outputDofs":                       [list(range(variables.n_muscle_spindles,variables.n_muscle_spindles+variables.n_interneurons))],          # [n_muscle_spindles,n_muscle_spindles+1,...,n_muscle_spindles+n_interneurons]
                    "callback":                         variables.callback_interneurons_to_motoneurons,
                    #"thresholdValue":                   20,                    # if mode is "localSetIfAboveThreshold", this is the threshold, if the value is above it, set the value `valueToSet`
                    #"valueToSet":                       20,                   # if mode is "localSetIfAboveThreshold", this is the value to set the target dof to, if the source dof is above thresholdValue.
                  }
                ],
                       
                # interneuron solver
                "Heun" : {
                  "description":                  "interneurons",
                  "timeStepWidth":                variables.dt_interneuron,
                  "logTimeStepWidthAsKey":        "dt_interneuron",
                  "durationLogKey":               "duration_interneuron",
                  "initialValues":                [],
                  "timeStepOutputInterval":       1e4,
                  "inputMeshIsGlobal":            True,
                  "dirichletBoundaryConditions":  {},
                  "dirichletOutputFilename":      None,             # filename for a vtp file that contains the Dirichlet boundary condition nodes and their values, set to None to disable
                  "checkForNanInf":               True,             # check if the solution vector contains nan or +/-inf values, if yes, an error is printed. This is a time-consuming check.
                  "nAdditionalFieldVariables":    0,
                  "additionalSlotNames":          [],
                      
                  # cellml model of interneurons
                  "CellML" : {
                    "modelFilename":                          variables.interneuron_cellml_file,              # input C++ source file or cellml XML file
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
                    
                    "mappings":                               variables.interneuron_mappings,                 # mappings between parameters and algebraics/constants and between outputConnectorSlots and states, algebraics or parameters, they are defined in helper.py
                    "parametersInitialValues":                variables.interneuron_parameters_initial_values, # # initial values for the parameters, either once for all instances, or for all instances in array of struct ordering with nParameters_ parameters per instance: [inst0p0, inst0p1, ... inst0pn, inst1p0, inst1p1, ...]
                    
                    "meshName":                               "interneuronMesh",                           
                    "stimulationLogFilename":                 "out/stimulation.log",

                    # output writer for states, algebraics and parameters                
                    "OutputWriter" : [
                      {"format": "PythonFile", "outputInterval": (int)(2./variables.dt_interneuron*variables.output_timestep_neurons), "filename": "out/interneurons", "binary": True, "fixedFormat": False, "combineFiles": True, "onlyNodalValues": True, "fileNumbering": "incremental"}
                    ]
                  }
                }
              }
            }
          }
        }
      }
    },
    
    # motoneuron with electro-mechanics
    "Term2": {
      
      # motoneuron with electro-mechanics
      "Coupling": {
        "description":            "motoneuron with electro-mechanics",    # description that will be shown in solver structure visualization
        "endTime":                variables.end_time,
        "timeStepWidth":          variables.dt_stimulation_check,
        "logTimeStepWidthAsKey":  "dt_stimulation_check",
        "durationLogKey":         "duration_multidomain",
        "timeStepOutputInterval": 100,
        "connectedSlotsTerm1To2": None,
        "connectedSlotsTerm2To1": None,
        
        # motoneuron
        "Term1": { 
          # mapping motor neuron signals + cortical input to actual inputs
          "MapDofs": {
            "description":                "motoneurons_input",   # description that will be shown in solver structure visualization
            "nAdditionalFieldVariables":  1,                              # number of additional field variables that are defined by this object. They have 1 component, use the templated function space and mesh given by meshName.
            "additionalSlotNames":        ["msin_m"],
            "meshName":                   "muscleSpindleAndInterneuronMesh",               # the mesh on which the additional field variables will be defined
            "beforeComputation": [                                        # transfer/mapping of dofs that will be performed before the computation of the nested solver
              {                                                 
                "fromConnectorSlot":                2,
                "toConnectorSlots":                 1,
                "fromSlotConnectorArrayIndex":      0,
                "toSlotConnectorArrayIndex":        0,
                "mode":                             "callback",          # "copyLocal", "copyLocalIfPositive", "localSetIfAboveThreshold" or "communicate"
                "fromDofNosNumbering":              "local",
                "toDofNosNumbering":                "local",
                "dofsMapping":                      None,
                "inputDofs":                        list(range(variables.n_muscle_spindles+variables.n_interneurons)),   # [0,1,...,n_muscle_spindles,...,n_muscle_spindles+n_interneurons]
                "outputDofs":                       [list(range(variables.n_motoneurons))],   # [0,1,...,n_motoneurons]
                "callback":                         variables.callback_motoneurons_input,
                #"thresholdValue":                   20,                   # if mode is "localSetIfAboveThreshold", this is the threshold, if the value is above it, set the value `valueToSet`
                #"valueToSet":                       20,                   # if mode is "localSetIfAboveThreshold", this is the value to set the target dof to, if the source dof is above thresholdValue.
              }
            ],
            "afterComputation": None,
                            
            # motoneuron solver
            "Heun" : {
              "description":                  "motoneurons",
              "timeStepWidth":                variables.dt_motoneuron,
              "logTimeStepWidthAsKey":        "dt_motoneuron",
              "durationLogKey":               "duration_motoneuron",
              "initialValues":                [],
              "timeStepOutputInterval":       1e4,
              "inputMeshIsGlobal":            True,
              "dirichletBoundaryConditions":  {},
              "dirichletOutputFilename":      None,             # filename for a vtp file that contains the Dirichlet boundary condition nodes and their values, set to None to disable
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
                  {"format": "PythonFile", "outputInterval": (int)(2./variables.dt_motoneuron*variables.output_timestep_motoneuron), "filename": "out/motoneurons", "binary": True, "fixedFormat": False, "combineFiles": True, "onlyNodalValues": True, "fileNumbering": "incremental"}
                ]
              }
            }
          }
        },
        
        # electrophysiology and mechanics solver
        "Term2": {
          # map from λ in the 3D mesh to muscle spindles input
          "MapDofs": {
            "description":                "muscle_spindles_input", # description that will be shown in solver structure visualization
            "nAdditionalFieldVariables":  1,                              # number of additional field variables that are defined by this object. They have 1 component, use the templated function space and mesh given by meshName.
            "additionalSlotNames":        ["ms"],
            "meshName":                   "muscleSpindleMesh",            # the mesh on which the additional field variables will be defined
            "beforeComputation": [                                        # transfer/mapping of dofs that will be performed before the computation of the nested solver
              {                                                 
                "fromConnectorSlot":                1,
                "toConnectorSlots":                 4,
                "fromSlotConnectorArrayIndex":      0,                    # which fiber/compartment, this does not matter here because all compartment meshes have the same displacements
                "toSlotConnectorArrayIndex":        0,
                "mode":                             "callback",           # "copyLocal", "copyLocalIfPositive", "localSetIfAboveThreshold" or "communicate"
                "fromDofNosNumbering":              "global",
                "toDofNosNumbering":                "local",
                "dofsMapping":                      None,
                "inputDofs":                        None,          # nodes are set at bottom in helper.py
                "outputDofs":                       [list(range(variables.n_muscle_spindles))],   # [0,1,...,n_muscle_spindles]
                "callback":                         variables.callback_muscle_spindles_input,
                #"thresholdValue":                   20,                  # if mode is "localSetIfAboveThreshold", this is the threshold, if the value is above it, set the value `valueToSet`
                #"valueToSet":                       20,                  # if mode is "localSetIfAboveThreshold", this is the value to set the target dof to, if the source dof is above thresholdValue.
              }
            ],
            "afterComputation":  None,
              
            # map from λ in the 3D mesh to golgi tendon organs
            "MapDofs": {
              "description":                "golgi_tendon_organs_input",      # description that will be shown in solver structure visualization
              "nAdditionalFieldVariables":  1,                              # number of additional field variables that are defined by this object. They have 1 component, use the templated function space and mesh given by meshName.
              "additionalSlotNames":        ["gt"],
              "meshName":                   "golgiTendonOrganMesh",               # the mesh on which the additional field variables will be defined
              "beforeComputation":          None, 
              "afterComputation": [                                        # transfer/mapping of dofs that will be performed before the computation of the nested solver
                {                                                 
                  "fromConnectorSlot":                1,
                  "toConnectorSlots":                 3,
                  "fromSlotConnectorArrayIndex":      0,                   # which fiber/compartment, this does not matter here because all compartment meshes have the same displacements
                  "toSlotConnectorArrayIndex":        0,
                  "mode":                             "callback",          # "copyLocal", "copyLocalIfPositive", "localSetIfAboveThreshold" or "communicate"
                  "fromDofNosNumbering":              "global",
                  "toDofNosNumbering":                "local",
                  "dofsMapping":                      None,
                  "inputDofs":                        None,          # nodes are set at bottom in helper.py
                  "outputDofs":                       [list(range(variables.n_golgi_tendon_organs))],   # [0,1,...,n_golgi_tendon_organs]
                  "callback":                         variables.callback_golgi_tendon_organs_input,
                  #"thresholdValue":                   20,                 # if mode is "localSetIfAboveThreshold", this is the threshold, if the value is above it, set the value `valueToSet`
                  #"valueToSet":                       20,                 # if mode is "localSetIfAboveThreshold", this is the value to set the target dof to, if the source dof is above thresholdValue.
                }
              ],
                  
              # map from motoneuronMesh to stimulated nodes
              "MapDofs": {
                "description":                "motoneurons->stimulated nodes",  # description that will be shown in solver structure visualization
                "nAdditionalFieldVariables":  1,                              # number of additional field variables that are defined by this object. They have 1 component, use the templated function space and mesh given by meshName.
                "additionalSlotNames":        ["mn"],
                "meshName":                   "motoneuronMesh",               # the mesh on which the additional field variables will be defined
                "beforeComputation":          None, 
                "afterComputation": [                                        # transfer/mapping of dofs that will be performed before the computation of the nested solver
                  {
                    "fromConnectorSlot":                2,
                    "toConnectorSlots":                 0,
                    "fromSlotConnectorArrayIndex":      0,
                    "toSlotConnectorArrayIndex":        compartment_no,      # which motor unit
                    "mode":                             "localSetIfAboveThreshold",          # "copyLocal", "copyLocalIfPositive", "localSetIfAboveThreshold" or "communicate"
                    "fromDofNosNumbering":              "local",
                    "toDofNosNumbering":                "global",
                    "dofsMapping":                      None,
                    "inputDofs":                        compartment_no,
                    "outputDofs":                       None,
                    #"callback":                         variables.callback_motoneuron_output,
                    "thresholdValue":                   20,                    # if mode is "localSetIfAboveThreshold", this is the threshold, if the value is above it, set the value `valueToSet`
                    "valueToSet":                       20,                   # if mode is "localSetIfAboveThreshold", this is the value to set the target dof to, if the source dof is above thresholdValue.
                  }
                for compartment_no in range(variables.n_compartments)],
                  
                # electro-mechanics solver
                "PrescribedValues": {
                  "meshName":               "3Dmesh",      # reference to the multidomain mesh
                  "numberTimeSteps":        1,             # number of timesteps to call the callback functions subsequently, this is usually 1 for prescribed values, because it is enough to set the reaction term only once per time step
                  "slotNames":              [],
                  "timeStepOutputInterval": 20,            # if the time step should be written to console, a value > 10 produces no output
                  
                  # a list of field variables that will get values assigned in every timestep, by the provided callback function
                  "fieldVariables1": [
                    {"name": "Vm_input", "callback": None},
                    {"name": "stress", "callback": None},
                  ],
                  "fieldVariables2":     [],
                  "additionalArgument":  None,         # a custom argument to the fieldVariables callback functions, this will be passed on as the last argument
                  
                  "OutputWriter" : [
                    #{"format": "Paraview", "outputInterval": 1, "filename": "out/" + variables.scenario_name + "/prescribed_stress{}".format(compartment_no), "binary": True, "fixedFormat": False, "combineFiles": True, "fileNumbering": "incremental"}
                  ]
                }
              }
            }
          }
        }
      }
    }
  }
}

with open("config.py", "w") as f:
  f.write(str(config))

# stop timer and calculate how long parsing lasted
if rank_no == 0:
  t_stop_script = timeit.default_timer()
  print("Python config parsed in {:.1f}s.".format(t_stop_script - t_start_script))

