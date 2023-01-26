import sys, os
import timeit
import argparse
import importlib
import distutils.util

# parse rank arguments
rank_no = int(sys.argv[-2])
n_ranks = int(sys.argv[-1])

# add variables subfolder to python path where the variables script is located
script_path = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, script_path)
sys.path.insert(0, os.path.join(script_path,'variables'))

import variables              # file variables.py, defined default values for all parameters, you can set the parameters there


variables.n_subdomains = variables.n_subdomains_x*variables.n_subdomains_y*variables.n_subdomains_z



# automatically initialize partitioning if it has not been set
if n_ranks != variables.n_subdomains:

  # create all possible partitionings to the given number of ranks
  optimal_value = n_ranks**(1/3)
  possible_partitionings = []
  for i in range(1,n_ranks+1):
    for j in range(1,n_ranks+1):
      if i*j <= n_ranks and n_ranks % (i*j) == 0:
        k = int(n_ranks / (i*j))
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
  print("cellml_file:             {}".format(variables.cellml_file))
  print("firing_times_file:       {}".format(variables.firing_times_file))
  print("********************************************************************************")

  # start timer to measure duration of parsing of this script
  t_start_script = timeit.default_timer()

# initialize all helper variables
from helper import *

variables.n_subdomains_xy = variables.n_subdomains_x * variables.n_subdomains_y
variables.n_fibers_total = variables.n_fibers_x * variables.n_fibers_y



#### set Dirichlet BC and Neumann BC for the free side of the muscle

[nx, ny, nz] = [elem + 1 for elem in variables.n_elements_muscle1]
[mx, my, mz] = [elem // 2 for elem in variables.n_elements_muscle1] # quadratic elements consist of 2 linear elements along each axis

k = 0 #free side of the muscle

for j in range(ny):
    for i in range(nx):
      variables.elasticity_dirichlet_bc[k*nx*ny + j*nx + i] = [0.0, 0.0, 0.0, None,None,None] # displacement ux uy uz, velocity vx vy vz


# meshes
variables.meshes.update(fiber_meshes)


def dbg(x, name=None):
  if name:
    print(name, end=': ')
  print(x)
  return x

# define the config dict
config = {
"scenarioName":                   variables.scenario_name,    # scenario nam
"logFormat":                      "csv",                      # "csv" or "json", format of the lines in the log file, csv gives smaller files
"solverStructureDiagramFile":     "solver_structure.txt",     # output file of a diagram that shows data connection between solvers
"mappingsBetweenMeshesLogFile":   "out/" + "/mappings_between_meshes.txt",
"meta": {                 # additional fields that will appear in the log
  "partitioning": [variables.n_subdomains_x, variables.n_subdomains_y, variables.n_subdomains_z]
},
"Meshes":                variables.meshes,
"Solvers": {
  "diffusionTermSolver": {# solver for the implicit timestepping scheme of the diffusion time step
    "relativeTolerance":  variables.diffusion_solver_reltol,
    "absoluteTolerance":  1e-10,         # 1e-10 absolute tolerance of the residual
    "maxIterations":      variables.diffusion_solver_maxit,
    "solverType":         variables.diffusion_solver_type,
    "preconditionerType": variables.diffusion_preconditioner_type,
    "dumpFilename":       "",   # "out/dump_"
    "dumpFormat":         "matlab",
  },
},
"Coupling": {
    "timeStepWidth":          variables.dt_elasticity,  # 1e-1
    "logTimeStepWidthAsKey":  "dt_splitting",
    "durationLogKey":         "duration_total",
    "timeStepOutputInterval": 1,
    "endTime":                variables.end_time,
    "connectedSlotsTerm1To2": None,        
    "connectedSlotsTerm2To1":  None,       
    "Term1": {        # monodomain, fibers
      "MultipleInstances": {
      "logKey":                     "duration_subdomains_xy",
      "ranksAllComputedInstances":  list(range(n_ranks)),
      "nInstances":                 variables.n_subdomains_xy,
      "instances":
      [{
      "ranks": list(range(subdomain_coordinate_y*variables.n_subdomains_x + subdomain_coordinate_x, n_ranks, variables.n_subdomains_x*variables.n_subdomains_y)),
      "StrangSplitting": {
        # "numberTimeSteps":        1,             # number of timesteps to call the callback functions subsequently, this is usually 1 for prescribed values, because it is enough to set the reaction term only once per time step
        "timeStepWidth":          variables.dt_splitting,  # 1e-1
        "endTime":                variables.end_time,
        "logTimeStepWidthAsKey":  "dt_splitting",
        "durationLogKey":         "duration_monodomain",
        "timeStepOutputInterval": 100,
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
                "timeStepWidth":                variables.dt_0D,                         # timestep width of 0D problem
                "logTimeStepWidthAsKey":        "dt_0D",                                 # key under which the time step width will be written to the log file
                "durationLogKey":               "duration_0D",                           # log key of duration for this solver
                "timeStepOutputInterval":       1e4,                                     # how often to print the current timestep
                "initialValues":                [],                                      # no initial values are specified
                "dirichletBoundaryConditions":  {},                                      # no Dirichlet boundary conditions are specified
                "dirichletOutputFilename":      None,                                    # filename for a vtp file that contains the Dirichlet boundary condition nodes and their values, set to None to disable

                "inputMeshIsGlobal":            True,                                    # the boundary conditions and initial values would be given as global numbers
                "checkForNanInf":               True,                                    # abort execution if the solution contains nan or inf values
                "nAdditionalFieldVariables":    0,                                       # number of additional field variables
                "additionalSlotNames":          [],                                      # names for the additional slots

                "CellML" : {
                  "modelFilename":                          variables.cellml_file,                          # input C++ source file or cellml XML file
                  #"statesInitialValues":                   [],                                             # if given, the initial values for the the states of one instance
                  "initializeStatesToEquilibrium":          False,                                          # if the equilibrium values of the states should be computed before the simulation starts
                  "initializeStatesToEquilibriumTimestepWidth": 1e-4,                                       # if initializeStatesToEquilibrium is enable, the timestep width to use to solve the equilibrium equation

                  # optimization parameters
                  "optimizationType":                       variables.optimization_type,           # "vc", "simd", "openmp" or "gpu", type of generated optimizated source file
                  "approximateExponentialFunction":         True,                                           # if optimizationType is "vc" or "gpu", whether the exponential function exp(x) should be approximate by (1+x/n)^n with n=1024
                  "compilerFlags":                          "-fPIC -O3 -march=native -shared ",             # compiler flags used to compile the optimized model code
                  "maximumNumberOfThreads":                 0,                                              # if optimizationType is "openmp", the maximum number of threads to use. Default value 0 means no restriction.

                  # stimulation callbacks
                  #"libraryFilename":                       "cellml_simd_lib.so",                           # compiled library
                  #"setSpecificParametersFunction":         set_specific_parameters,                        # callback function that sets parameters like stimulation current
                  #"setSpecificParametersCallInterval":     int(1./variables.stimulation_frequency/variables.dt_0D),         # set_specific_parameters should be called every 0.1, 5e-5 * 1e3 = 5e-2 = 0.05
                  "setSpecificStatesFunction":              None,                                             # callback function that sets states like Vm, activation can be implemented by using this method and directly setting Vm values, or by using setParameters/setSpecificParameters
                  #"setSpecificStatesCallInterval":         2*int(1./variables.stimulation_frequency/variables.dt_0D),       # set_specific_states should be called variables.stimulation_frequency times per ms, the factor 2 is needed because every Heun step includes two calls to rhs
                  "setSpecificStatesCallInterval":          0,                                                               # 0 means disabled
                  "setSpecificStatesCallFrequency":         variables.get_specific_states_call_frequency(fiber_no, motor_unit_no),   # set_specific_states should be called variables.stimulation_frequency times per ms
                  "setSpecificStatesFrequencyJitter":       variables.get_specific_states_frequency_jitter(fiber_no, motor_unit_no), # random value to add or substract to setSpecificStatesCallFrequency every stimulation, this is to add random jitter to the frequency
                  "setSpecificStatesRepeatAfterFirstCall":  variables.dt_elasticity,                                                            # [ms] simulation time span for which the setSpecificStates callback will be called after a call was triggered
                  "setSpecificStatesCallEnableBegin":       variables.get_specific_states_call_enable_begin(fiber_no, motor_unit_no),# [ms] first time when to call setSpecificStates
                  "additionalArgument":                     fiber_no,                                       # last argument that will be passed to the callback functions set_specific_states, set_specific_parameters, etc.

                  # parameters to the cellml model
                  "mappings":                               variables.muscle1_mappings,                             # mappings between parameters and algebraics/constants and between outputConnectorSlots and states, algebraics or parameters, they are defined in helper.py
                  "parametersInitialValues":                variables.parameters_initial_values,            #[0.0, 1.0],      # initial values for the parameters: I_Stim, l_hs

                  "meshName":                               "MeshFiber_{}".format(fiber_no),                # reference to the fiber mesh
                  "stimulationLogFilename":                 "out/stimulation.log",                          # a file that will contain the times of stimulations
                },
                "OutputWriter" : [
                  {"format": "Paraview", "outputInterval": 1, "filename": "out/" + "/0D_states({},{})".format(fiber_in_subdomain_coordinate_x,fiber_in_subdomain_coordinate_y), "binary": True, "fixedFormat": False, "combineFiles": True}
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
              "ranks":                         list(range(variables.n_subdomains_z)),   # these rank nos are local nos to the outer instance of MultipleInstances, i.e. from 0 to number of ranks in z direction
              "ImplicitEuler" : {
                "initialValues":               [],                                      # no initial values are given
                #"numberTimeSteps":            1,
                "timeStepWidth":               variables.dt_1D,                         # timestep width for the diffusion problem
                "timeStepWidthRelativeTolerance": 1e-10,
                "logTimeStepWidthAsKey":       "dt_1D",                                 # key under which the time step width will be written to the log file
                "durationLogKey":              "duration_1D",                           # log key of duration for this solver
                "timeStepOutputInterval":      1e8,                                     # how often to print the current timestep
                "dirichletBoundaryConditions": {},                                      # old Dirichlet BC that are not used in FastMonodomainSolver: {0: -75.0036, -1: -75.0036},
                "dirichletOutputFilename":     None,                                    # filename for a vtp file that contains the Dirichlet boundary condition nodes and their values, set to None to disable
                "inputMeshIsGlobal":           True,                                    # initial values would be given as global numbers
                "solverName":                  "diffusionTermSolver",                   # reference to the linear solver
                "nAdditionalFieldVariables":   2,
                "additionalSlotNames":         [],
                "checkForNanInf":              False,

                "FiniteElementMethod" : {
                  "inputMeshIsGlobal":         True,
                  "meshName":                  "MeshFiber_{}".format(fiber_no),
                  "solverName":                "diffusionTermSolver",
                  "prefactor":                 get_diffusion_prefactor(fiber_no, motor_unit_no),  # resolves to Conductivity / (Am * Cm)
                  "slotName":                  "vm",
                },
                "OutputWriter" : [
                ]
              },
            } for fiber_in_subdomain_coordinate_y in range(n_fibers_in_subdomain_y(subdomain_coordinate_y)) \
                for fiber_in_subdomain_coordinate_x in range(n_fibers_in_subdomain_x(subdomain_coordinate_x)) \
                  for fiber_no in [get_fiber_no(subdomain_coordinate_x, subdomain_coordinate_y, fiber_in_subdomain_coordinate_x, fiber_in_subdomain_coordinate_y)] \
                    for motor_unit_no in [get_motor_unit_no(fiber_no)]],
            "OutputWriter" : [
              {"format": "Paraview", "outputInterval": 1, "filename": "out/" + "/CoupledMeshFibers", "binary": True, "fixedFormat": False, "onlyNodalValues":True,  "combineFiles": True, "fileNumbering": "incremental"}
            ]
          },
        },
      }
      } if (subdomain_coordinate_x,subdomain_coordinate_y) == (variables.own_subdomain_coordinate_x,variables.own_subdomain_coordinate_y) else None
      for subdomain_coordinate_y in range(variables.n_subdomains_y)
        for subdomain_coordinate_x in range(variables.n_subdomains_x)]
      },
      "fiberDistributionFile":    variables.fiber_distribution_file,   # for FastMonodomainSolver, e.g. MU_fibre_distribution_3780.txt
      "firingTimesFile":          variables.firing_times_file,         # for FastMonodomainSolver, e.g. MU_firing_times_real.txt
      "onlyComputeIfHasBeenStimulated": variables.fast_monodomain_solver_optimizations,                          # only compute fibers after they have been stimulated for the first time
      "disableComputationWhenStatesAreCloseToEquilibrium": variables.fast_monodomain_solver_optimizations,       # optimization where states that are close to their equilibrium will not be computed again
      "valueForStimulatedPoint":  variables.vm_value_stimulated,       # to which value of Vm the stimulated node should be set
      "neuromuscularJunctionRelativeSize": 0.1,                          # range where the neuromuscular junction is located around the center, relative to fiber length. The actual position is draws randomly from the interval [0.5-s/2, 0.5+s/2) with s being this option. 0 means sharply at the center, 0.1 means located approximately at the center, but it can vary 10% in total between all fibers.
      "generateGPUSource":        True,                                # (set to True) only effective if optimizationType=="gpu", whether the source code for the GPU should be generated. If False, an existing source code file (which has to have the correct name) is used and compiled, i.e. the code generator is bypassed. This is useful for debugging, such that you can adjust the source code yourself. (You can also add "-g -save-temps " to compilerFlags under CellMLAdapter)
      "useSinglePrecision":       False,                               # only effective if optimizationType=="gpu", whether single precision computation should be used on the GPU. Some GPUs have poor double precision performance. Note, this drastically increases the error and, in consequence, the timestep widths should be reduced.
      #"preCompileCommand":        "bash -c 'module load argon-tesla/gcc/11-20210110-openmp; module list; gcc --version",     # only effective if optimizationType=="gpu", system command to be executed right before the compilation
      #"postCompileCommand":       "'",   # only effective if optimizationType=="gpu", system command to be executed right after the compilation
  }  
}
}

# stop   timer and calculate how long parsing lasted
if rank_no == 0:
  t_stop_script = timeit.default_timer()
  print("Python config parsed in {:.1f}s.".format(t_stop_script - t_start_script))
