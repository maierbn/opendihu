# Multicompartment 3D, biceps
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
    print("Warning: There is no variables file, e.g:\n ./multidomain_with_fat ../settings_multidomain_with_fat.py ramp_emg.py\n")
  exit(0)
  
# -------------------------------------------------------- begin parameters ---------------------------------------------------------
# -------------------------------------------------------- end parameters ------------------------------------------------------------
parser = argparse.ArgumentParser(description='multidomain_with_fat')
parser.add_argument('--n_subdomains', nargs=3,               help='Number of subdomains in x,y,z direction.',    type=int)
parser.add_argument('--n_subdomains_x', '-x',                help='Number of subdomains in x direction.',        type=int, default=variables.n_subdomains_x)
parser.add_argument('--n_subdomains_y', '-y',                help='Number of subdomains in y direction.',        type=int, default=variables.n_subdomains_y)
parser.add_argument('--n_subdomains_z', '-z',                help='Number of subdomains in z direction.',        type=int, default=variables.n_subdomains_z)
parser.add_argument('--paraview_output',                     help='Enable the paraview output writer.',          default=variables.paraview_output, action='store_true')
parser.add_argument('--adios_output',                        help='Enable the MegaMol/ADIOS output writer.',     default=variables.adios_output, action='store_true')
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
parser.add_argument('--sampling_stride_x',                   help='Stride to select the mesh points in x direction.',         type=int, default=variables.sampling_stride_x)
parser.add_argument('--sampling_stride_y',                   help='Stride to select the mesh points in y direction.',         type=int, default=variables.sampling_stride_y)
parser.add_argument('--sampling_stride_z',                   help='Stride to select the mesh points in z direction.',         type=int, default=variables.sampling_stride_z)
parser.add_argument('--dt_0D',                               help='Timestep width of subcellular problem.',                   type=float, default=variables.dt_0D)
parser.add_argument('--dt_multidomain',                      help='Timestep width of multidomain (diffusion) problem.',       type=float, default=variables.dt_multidomain)
parser.add_argument('--dt_splitting',                        help='Timestep width of Strang splitting.',                      type=float, default=variables.dt_splitting)

if variables.scenario_name:
  parser.add_argument('--scenario_name',                       help='Name of the scenario in the log file.',                    type=str, default=variables.scenario_name)

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

if variables.scenario_name is None:
  variables.scenario_name = "{}_{},hypre_dt{}_atol{}_rtol{}_theta{}_sym{}_lump{}_{}mus".format(variables.multidomain_solver_type, variables.multidomain_preconditioner_type, variables.dt_splitting, variables.multidomain_absolute_tolerance, variables.multidomain_relative_tolerance, variables.theta, variables.use_symmetric_preconditioner_matrix, variables.use_lumped_mass_matrix, len(variables.motor_units))

if variables.initial_guess_nonzero is None:
  variables.initial_guess_nonzero = variables.multidomain_preconditioner_type != "lu"

if variables.multidomain_preconditioner_type == "boomeramg":
  variables.multidomain_max_iterations = 1000        # if boomeramg has bad convergence, abort after 100 iterations and use alternative preconditioner

# output information of run
if rank_no == 0:
  print("scenario_name: {},  n_subdomains: {} {} {},  n_ranks: {},  end_time: {}".format(variables.scenario_name, variables.n_subdomains_x, variables.n_subdomains_y, variables.n_subdomains_z, n_ranks, variables.end_time))
  print("dt_0D:           {:0.0e}    multidomain solver:         {} it. of {} ({} it. of {}), lumped mass matrix: {}, initial guess: {}".format(variables.dt_0D, int(variables.multidomain_max_iterations), variables.multidomain_solver_type, int(variables.multidomain_alternative_solver_max_iterations), variables.multidomain_alternative_solver_type, variables.use_lumped_mass_matrix, "0" if not variables.initial_guess_nonzero else "previous solution"))
  print("dt_multidomain:  {:0.0e}    multidomain preconditioner: {} ({}), symmetric precond.: {}".format(variables.dt_multidomain, variables.multidomain_preconditioner_type, variables.multidomain_alternative_preconditioner_type, variables.use_symmetric_preconditioner_matrix))
  print("dt_splitting:    {:0.0e}    theta: {}, solver tolerances, abs: {}, rel: {}".format(variables.dt_splitting, variables.theta, variables.multidomain_absolute_tolerance, variables.multidomain_relative_tolerance))
  print("fiber_file:              {}".format(variables.fiber_file))
  print("fat_mesh_file:           {}".format(variables.fat_mesh_file))
  print("cellml_file:             {}".format(variables.cellml_file))
  print("firing_times_file:       {}".format(variables.firing_times_file))
  print("********************************************************************************")
  
  # start timer to measure duration of parsing of this script  
  t_start_script = timeit.default_timer()
    
# initialize all helper variables
from helper import *

# settings for the multidomain solver
multidomain_solver = {
  "timeStepWidth":                    variables.dt_multidomain,             # time step width of the subcellular problem
  "endTime":                          variables.end_time,                   # end time, this is not relevant because it will be overridden by the splitting scheme
  "timeStepOutputInterval":           100,                                  # how often the output timestep should be printed
  "durationLogKey":                   "duration_multidomain",               # key for duration in log.csv file
  "slotNames":                        ["vm_old", "vm_new", "g_mu", "g_tot"],  # names of the data connector slots, maximum length per name is 10 characters. g_mu is gamma (active stress) of the compartment, g_tot is the total gamma
  
  # material parameters for the compartments
  "nCompartments":                    variables.n_compartments,             # number of compartments
  "compartmentRelativeFactors":       variables.relative_factors.tolist(),  # list of lists of (the factors for all dofs), because "inputIsGlobal": True, this contains the global dofs
  "inputIsGlobal":                    True,                                 # if values and dofs correspond to the global numbering
  "am":                               [variables.get_am(mu_no) for mu_no in range(variables.n_compartments)],   # Am parameter for every motor unit (ration of surface to volume of fibers)
  "cm":                               [variables.get_cm(mu_no) for mu_no in range(variables.n_compartments)],   # Cm parameter for every motor unit (capacitance of the cellular membrane)
  
  # solver options
  "solverName":                       "multidomainLinearSolver",            # reference to the solver used for the global linear system of the multidomain eq.
  "alternativeSolverName":            "multidomainAlternativeLinearSolver", # reference to the alternative solver, which is used when the normal solver diverges
  "subSolverType":                    "gamg",                               # sub solver when block jacobi preconditioner is used
  "subPreconditionerType":            "none",                               # sub preconditioner when block jacobi preconditioner is used
  #"subPreconditionerType":            "boomeramg",                          # sub preconditioner when block jacobi preconditioner is used, boomeramg is the AMG preconditioner of HYPRE

  "hypreOptions":                     "-pc_hypre_boomeramg_strong_threshold 0.7",       # additional options if a hypre preconditioner is selected
  "theta":                            variables.theta,                      # weighting factor of implicit term in Crank-Nicolson scheme, 0.5 gives the classic, 2nd-order Crank-Nicolson scheme, 1.0 gives implicit euler
  "useLumpedMassMatrix":              variables.use_lumped_mass_matrix,     # which formulation to use, the formulation with lumped mass matrix (True) is more stable but approximative, the other formulation (False) is exact but needs more iterations
  "useSymmetricPreconditionerMatrix": variables.use_symmetric_preconditioner_matrix,    # if the diagonal blocks of the system matrix should be used as preconditioner matrix
  "initialGuessNonzero":              variables.initial_guess_nonzero,      # if the initial guess for the 3D system should be set as the solution of the previous timestep, this only makes sense for iterative solvers
  "enableFatComputation":             True,                                 # disabling the computation of the fat layer is only for debugging and speeds up computation. If set to False, the respective matrix is set to the identity
  "showLinearSolverOutput":           variables.show_linear_solver_output,  # if convergence information of the linear solver in every timestep should be printed, this is a lot of output for fast computations
  "updateSystemMatrixEveryTimestep":  False,                                # if this multidomain solver will update the system matrix in every first timestep, us this only if the geometry changed, e.g. by contraction
  "recreateLinearSolverInterval":     0,                                    # how often the Petsc KSP object (linear solver) should be deleted and recreated. This is to remedy memory leaks in Petsc's implementation of some solvers. 0 means disabled.
  "rescaleRelativeFactors":           True,                                 # if all relative factors should be rescaled such that max Σf_r = 1  
  "setDirichletBoundaryConditionPhiE":False,                                # (set to False) if the last dof of the extracellular space (variable phi_e) should have a 0 Dirichlet boundary condition. However, this makes the solver converge slower.
  "setDirichletBoundaryConditionPhiB":False,                                # (set to False) if the last dof of the fat layer (variable phi_b) should have a 0 Dirichlet boundary condition. However, this makes the solver converge slower.
  "resetToAverageZeroPhiE":           True,                                 # if a constant should be added to the phi_e part of the solution vector after every solve, such that the average is zero
  "resetToAverageZeroPhiB":           True,                                 # if a constant should be added to the phi_b part of the solution vector after every solve, such that the average is zero
  
  "PotentialFlow": {
    "FiniteElementMethod" : {  
      "meshName":                     "3Dmesh",
      "solverName":                   "potentialFlowSolver",
      "prefactor":                    1.0,
      "dirichletBoundaryConditions":  variables.potential_flow_dirichlet_bc,
      "dirichletOutputFilename":      "out/" + variables.scenario_name + "/dirichlet_potential_flow",               # output filename for the dirichlet boundary conditions, set to "" to have no output
      "neumannBoundaryConditions":    [],
      "inputMeshIsGlobal":            True,
      "slotName":                     "",
    },
  },
  "Activation": {
    "FiniteElementMethod" : {  
      "meshName":                     "3Dmesh",
      "solverName":                   "multidomainLinearSolver",
      "prefactor":                    1.0,
      "inputMeshIsGlobal":            True,
      "dirichletBoundaryConditions":  {},
      "dirichletOutputFilename":      None,               # output filename for the dirichlet boundary conditions, set to "" to have no output
      "neumannBoundaryConditions":    [],
      "slotName":                     "",
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
      "dirichletOutputFilename":      None,               # output filename for the dirichlet boundary conditions, set to "" to have no output
      "neumannBoundaryConditions":    [],
      "slotName":                     "",
    },
  },
  
  "OutputWriter" : [
    {"format": "Paraview", "outputInterval": (int)(1./variables.dt_multidomain*variables.output_timestep_multidomain), "filename": "out/" + variables.scenario_name + "/multidomain", "binary": True, "fixedFormat": False, "combineFiles": True, "fileNumbering": "incremental"},
    #{"format": "ExFile", "filename": "out/fiber_"+str(i), "outputInterval": 1./dt_1D*output_timestep, "sphereSize": "0.02*0.02*0.02", "fileNumbering": "incremental"},
    #{"format": "PythonFile", "filename": "out/fiber_"+str(i), "outputInterval": int(1./dt_1D*output_timestep), "binary":True, "onlyNodalValues":True, "fileNumbering": "incremental"},
  ]
}
  
config = {
  "scenarioName":          variables.scenario_name,
  "logFormat":                      "csv",                                # "csv" or "json", format of the lines in the log file, csv gives smaller files
  "solverStructureDiagramFile":     "out/"+variables.scenario_name+"/solver_structure.txt",               # output file of a diagram that shows data connection between solvers
  "mappingsBetweenMeshesLogFile":   "out/"+variables.scenario_name+"/mappings_between_meshes_log.txt",    # log file for mappings 
  "meta": {                                                               # additional fields that will appear in the log
    "partitioning":         [variables.n_subdomains_x, variables.n_subdomains_y, variables.n_subdomains_z]
  },
  "Meshes":                variables.meshes,
  "MappingsBetweenMeshes": variables.mappings_between_meshes,
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
    "multidomainLinearSolver": {
      "relativeTolerance":  variables.multidomain_relative_tolerance,
      "absoluteTolerance":  variables.multidomain_absolute_tolerance,         # 1e-10 absolute tolerance of the residual          
      "maxIterations":      variables.multidomain_max_iterations,             # maximum number of iterations
      "solverType":         variables.multidomain_solver_type,
      "preconditionerType": variables.multidomain_preconditioner_type,
      "hypreOptions":       "-pc_hypre_boomeramg_strong_threshold 0.7",       # additional options if a hypre preconditioner is selected
      "dumpFormat":         "matlab",
      "dumpFilename":       "",
          
      # gamg specific options:
      "gamgType":                         "agg",                          # one of agg, geo, or classical 
      "cycleType":                        "cycleW",                             # either cycleV or cycleW
      "nLevels":                          25,      
    },
    "multidomainAlternativeLinearSolver": {                                   # the alternative solver is used when the normal solver diverges
      "relativeTolerance":  variables.multidomain_relative_tolerance,
      "absoluteTolerance":  variables.multidomain_absolute_tolerance,         # 1e-10 absolute tolerance of the residual          
      "maxIterations":      variables.multidomain_alternative_solver_max_iterations,   # maximum number of iterations
      "solverType":         variables.multidomain_alternative_solver_type,
      "preconditionerType": variables.multidomain_alternative_preconditioner_type,
      "hypreOptions":       "",                                   # additional options if a hypre preconditioner is selected
      "dumpFormat":         "matlab",
      "dumpFilename":       "",
    }
  },
  "StrangSplitting": {
    "timeStepWidth":          variables.dt_splitting,
    "logTimeStepWidthAsKey":  "dt_splitting",
    "durationLogKey":         "duration_total",
    "timeStepOutputInterval": 100,
    "endTime":                variables.end_time,
    "connectedSlotsTerm1To2": [0],          # CellML V_mk (0) <=> Multidomain V_mk^(i) (0)
    "connectedSlotsTerm2To1": [None, 0],    # Multidomain V_mk^(i+1) (1) -> CellML V_mk (0)

    "Term1": {      # CellML
      "MultipleInstances": {
        "nInstances": variables.n_compartments,  
        "instances": [        # settings for each motor unit, `i` is the index of the motor unit
        {
          "ranks": list(range(n_ranks)),
          "Heun" : {
            "timeStepWidth":                variables.dt_0D,  # 5e-5
            "logTimeStepWidthAsKey":        "dt_0D",
            "durationLogKey":               "duration_0D",
            "initialValues":                [],
            "timeStepOutputInterval":       100,
            "inputMeshIsGlobal":            True,
            "dirichletBoundaryConditions":  {},
            "dirichletOutputFilename":      None,               # output filename for the dirichlet boundary conditions, set to "" to have no output
            "checkForNanInf":               True,             # check if the solution vector contains nan or +/-inf values, if yes, an error is printed. This is a time-consuming check.
            "nAdditionalFieldVariables":    0,
            "additionalSlotNames":          [],               # slot names of the additional slots, maximum 10 characters per name
                
            "CellML" : {
              "modelFilename":                          variables.cellml_file,                          # input C++ source file or cellml XML file
              "initializeStatesToEquilibrium":          False,                                          # if the equilibrium values of the states should be computed before the simulation starts
              "initializeStatesToEquilibriumTimestepWidth": 1e-4,                                       # if initializeStatesToEquilibrium is enable, the timestep width to use to solve the equilibrium equation
              
              # optimization parameters
              "optimizationType":                       "vc",                                           # "vc", "simd", "openmp" type of generated optimizated source file
              "approximateExponentialFunction":         True,                                           # if optimizationType is "vc", whether the exponential function exp(x) should be approximate by (1+x/n)^n with n=1024
              "compilerFlags":                          "-fPIC -O3 -march=native -shared ",             # compiler flags used to compile the optimized model code
              "maximumNumberOfThreads":                 0,                                              # if optimizationType is "openmp", the maximum number of threads to use. Default value 0 means no restriction.
              
              # stimulation callbacks
              "setSpecificStatesFunction":              set_specific_states,                                             # callback function that sets states like Vm, activation can be implemented by using this method and directly setting Vm values, or by using setParameters/setSpecificParameters
              #"setSpecificStatesCallInterval":         2*int(1./variables.stimulation_frequency/variables.dt_0D),       # set_specific_states should be called variables.stimulation_frequency times per ms, the factor 2 is needed because every Heun step includes two calls to rhs
              "setSpecificStatesCallInterval":          0,                                                               # 0 means disabled
              "setSpecificStatesCallFrequency":         variables.get_specific_states_call_frequency(compartment_no),   # set_specific_states should be called variables.stimulation_frequency times per ms
              "setSpecificStatesFrequencyJitter":       variables.get_specific_states_frequency_jitter(compartment_no), # random value to add or substract to setSpecificStatesCallFrequency every stimulation, this is to add random jitter to the frequency
              "setSpecificStatesRepeatAfterFirstCall":  0.01,                                                            # [ms] simulation time span for which the setSpecificStates callback will be called after a call was triggered
              "setSpecificStatesCallEnableBegin":       variables.get_specific_states_call_enable_begin(compartment_no),# [ms] first time when to call setSpecificStates
              "additionalArgument":                     compartment_no,
              
              "mappings":                               variables.mappings,                             # mappings between parameters and algebraics/constants and between outputConnectorSlots and states, algebraics or parameters, they are defined in helper.py
              "parametersInitialValues":                variables.parameters_initial_values,            #[0.0, 1.0],      # initial values for the parameters: I_Stim, l_hs
              
              "meshName":                               "3Dmesh",
              "stimulationLogFilename":                 "out/" + variables.scenario_name + "/stimulation.log",
            },
						"OutputWriter" : [
              {"format": "Paraview", "outputInterval": (int)(1./variables.dt_0D*variables.output_timestep_0D_states), "filename": "out/" + variables.scenario_name + "/0D_states_compartment_{}".format(compartment_no), "binary": False, "fixedFormat": False, "combineFiles": True, "fileNumbering": "incremental"}
            ] if variables.states_output else []
          }
        } for compartment_no in range(variables.n_compartments)]
      },
    },
    "Term2": {     # Diffusion, i.e. Multidomain
      "MultidomainSolver" : multidomain_solver,
      "OutputSurface": {        # version for fibers_emg_2d_output
        "OutputWriter": [
          {"format": "Paraview",   "outputInterval": int(1./variables.dt_multidomain*variables.output_timestep_surface), "filename": "out/" + variables.scenario_name + "/surface_emg", "binary": True, "fixedFormat": False, "onlyNodalValues": True, "combineFiles": True, "fileNumbering": "incremental"},
          {"format": "PythonFile", "outputInterval": int(1./variables.dt_multidomain*variables.output_timestep_surface), "filename": "out/" + variables.scenario_name + "/surface_emg", "binary": True, "fixedFormat": False, "onlyNodalValues": True, "combineFiles": True, "fileNumbering": "incremental"},
        ],
        #"face":                    ["1+","0+"],         # which faces of the 3D mesh should be written into the 2D mesh
        "face":                     ["1+"],              # which faces of the 3D mesh should be written into the 2D mesh
        "samplingPoints":           variables.hdemg_electrode_positions,    # the electrode positions, they are created in the helper.py script
        "updatePointPositions":     False,               # the electrode points should be initialize in every timestep (set to False for the static case). This makes a difference if the muscle contracts, then True=fixed electrodes, False=electrodes moving with muscle.
        "filename":                 "out/{}/electrodes.csv".format(variables.scenario_name),
        "enableCsvFile":            True,                # if the values at the sampling points should be written to csv files
        "enableVtpFile":            False,               # if the values at the sampling points should be written to vtp files
        "enableGeometryInCsvFile":  False,               # if the csv output file should contain geometry of the electrodes in every time step. This increases the file size and only makes sense if the geometry changed throughout time, i.e. when computing with contraction
        "enableGeometryFiles":      False,               # if there should be extra files of the locations of the electrodes on every rank
        "xiTolerance":              0.3,                 # tolerance for element-local coordinates xi, for finding electrode positions inside the elements. Increase or decrease this numbers if not all electrode points are found.
        "MultidomainSolver":        multidomain_solver,
      }
    }
  }
}

# stop timer and calculate how long parsing lasted
if rank_no == 0:
  t_stop_script = timeit.default_timer()
  print("Python config parsed in {:.1f}s.".format(t_stop_script - t_start_script))

