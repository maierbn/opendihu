# Multiple 1D fibers (monodomain) with 3D dynamic mooney rivlin with active contraction term, on biceps geometry

import sys, os
import timeit
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
from helper import *

# modify variables according to specific scenario

variables.scenario_name = "muscle_left"

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
  print("dt_3D:           {:0.0e}, paraview_output: {}".format(variables.dt_3D, variables.paraview_output))
  print("output_timestep: {:0.0e}  stimulation_frequency: {} 1/ms = {} Hz".format(variables.output_timestep, variables.stimulation_frequency, variables.stimulation_frequency*1e3))
  print("cellml_file:             {}".format(variables.cellml_file))
  print("fiber_distribution_file: {}".format(variables.fiber_distribution_file))
  print("firing_times_file:       {}".format(variables.firing_times_file))
  print("********************************************************************************")
  
  print("prefactor: sigma_eff/(Am*Cm) = {} = {} / ({}*{})".format(variables.Conductivity/(variables.Am*variables.Cm), variables.Conductivity, variables.Am, variables.Cm))
  

# initialize all helper variables
variables.n_subdomains_xy = variables.n_subdomains_x * variables.n_subdomains_y
variables.n_fibers_total = variables.n_fibers_x * variables.n_fibers_y

# add meshes
meshes_muscle_left = {
  # no `nodePositions` fields as the nodes are created internally
  "muscle_left_Mesh": {
    "nElements" :         variables.n_elements_muscle,
    "physicalExtent":     variables.muscle_left_extent,
    "physicalOffset":     variables.muscle_left_offset,
    "logKey":             "muscle_left",
    "inputMeshIsGlobal":  True,
    "nRanks":             n_ranks
  },
  # needed for mechanics solver
  "muscle_left_Mesh_quadratic": {
    "nElements" :         [elems // 2 for elems in variables.n_elements_muscle],
    "physicalExtent":     variables.muscle_left_extent,
    "physicalOffset":     variables.muscle_left_offset,
    "logKey":             "muscle_left_quadratic",
    "inputMeshIsGlobal":  True,
    "nRanks":             n_ranks,
  }
}
variables.meshes.update(meshes_muscle_left)
variables.meshes.update(fiber_meshes)

#############################

# parameters for the main simulation
# ---------------------------------------------
# https://uni-stuttgart.conceptboard.com/board/i45e-9bz9-qzb0-ppek-n9s2
#          ----> z
#          +----------+-.  tendon  .-+----------+
#  u_z=0 & | muscle 1 |~:~~~~~~~~~~:~| muscle 2 | u_z=0 &
# one      +----------+-'    ^     '-+----------+  one
# edge                       |                     edge
# u=0                    u_x=u_y=0                 u=0
#

#### set Dirichlet BC and Neumann BC for the free side of the muscle

[nx, ny, nz] = [elem + 1 for elem in variables.n_elements_muscle]
[mx, my, mz] = [elem // 2 for elem in variables.n_elements_muscle] # quadratic elements consist of 2 linear elements along each axis

variables.elasticity_dirichlet_bc = {}

k = 0 #free side of the muscle

for j in range(ny):
    for i in range(nx):
      variables.elasticity_dirichlet_bc[k*nx*ny + j*nx + i] = [None,None,0.0, None,None,None] # displacement ux uy uz, velocity vx vy vz

# fix edge, note: the multidomain simulation does not work without this (linear solver finds no solution)
for i in range(nx):
    variables.elasticity_dirichlet_bc[k*nx*ny + 0*nx + i] = [0.0,0.0,0.0, None,None,None]
    
# fix corner completely
variables.elasticity_dirichlet_bc[k*nx*ny + 0] = [0.0,0.0,0.0, None,None,None]

# # guide right end of muscle along z axis
# # muscle mesh
# for j in range(ny):
#     for i in range(nx):
#       variables.muscle1_elasticity_dirichlet_bc[(nz-1)*nx*ny + j*nx + i] = [0.0,0.0,None, None,None,None]

# initial Neumann BC at bottom nodes, traction along z axis
# will be set by tendon
# k = 0 #0 or mz-1
# variables.force = 1.0
# traction_vector = [0, 0, -variables.force]     # the traction force in specified in the reference configuration
# face = "2-"
# variables.elasticity_neumann_bc = [{"element": k*mx*my + j*mx + i, "constantVector": traction_vector, "face": face} for j in range(my) for i in range(mx)]

# define the config dict
config = {
  "scenarioName":          variables.scenario_name,
  "logFormat":             "csv",
  "solverStructureDiagramFile":     "out/solver_structure.txt",     # output file of a diagram that shows data connection between solvers
  "mappingsBetweenMeshesLogFile":   "out/mappings_between_meshes.txt",  # log file of when mappings between meshes occur
  "Meshes":                variables.meshes,
  "MappingsBetweenMeshes": {"muscle_left_fiber_{}".format(f) : ["muscle_left_Mesh", "muscle_left_Mesh_quadratic"] for f in range(variables.n_fibers_total)},
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
    "mechanicsSolver": {   # solver for the dynamic mechanics problem
      "relativeTolerance":  1e-10,           # 1e-10 relative tolerance of the linear solver
      "absoluteTolerance":  1e-10,          # 1e-10 absolute tolerance of the residual of the linear solver
      "solverType":         "preonly",      # type of the linear solver: cg groppcg pipecg pipecgrr cgne nash stcg gltr richardson chebyshev gmres tcqmr fcg pipefcg bcgs ibcgs fbcgs fbcgsr bcgsl cgs tfqmr cr pipecr lsqr preonly qcg bicg fgmres pipefgmres minres symmlq lgmres lcd gcr pipegcr pgmres dgmres tsirm cgls
      "preconditionerType": "lu",           # type of the preconditioner
      "maxIterations":       1e4,           # maximum number of iterations in the linear solver
      "snesMaxFunctionEvaluations": 1e8,    # maximum number of function iterations
      "snesMaxIterations":   14,            # maximum number of iterations in the nonlinear solver
      "snesRelativeTolerance": 1e-5,        # relative tolerance of the nonlinear solver
      "snesAbsoluteTolerance": 1e-5,        # absolute tolerance of the nonlinear solver
      "snesLineSearchType": "l2",           # type of linesearch, possible values: "bt" "nleqerr" "basic" "l2" "cp" "ncglinear"
      "snesRebuildJacobianFrequency": 5,    # how often the jacobian should be recomputed, -1 indicates NEVER rebuild, 1 means rebuild every time the Jacobian is computed within a single nonlinear solve, 2 means every second time the Jacobian is built etc. -2 means rebuild at next chance but then never again 
      "dumpFilename":        "",            # dump system matrix and right hand side after every solve
      "dumpFormat":          "matlab",      # default, ascii, matlab
    }
  },
  "PreciceAdapter": {        # precice adapter for muscle
    "couplingEnabled":          True,
    "outputOnlyConvergedTimeSteps": False, #default is true
    "scalingFactor":            1,

    "timeStepOutputInterval":   100,                        # interval in which to display current timestep and time in console
    "timestepWidth":            1,                          # coupling time step width, must match the value in the precice config
    "preciceConfigFilename":    "precice_config_dirichlet_neumann_cubic_geometry.xml",    # the preCICE configuration file
    "preciceParticipantName":   "MuscleSolverLeft",             # name of the own precice participant, has to match the name given in the precice xml config file
    "preciceMeshes": [                                      # the precice meshes get created as the top or bottom surface of the main geometry mesh of the nested solver
      {
        "preciceMeshName":      "MuscleMeshLeft",         # precice name of the 2D coupling mesh
        "face":                 "2+",                       # face of the 3D mesh where the 2D mesh is located, "2-" = left, "2+" = right (z-coordinate)
      }
    ],
    "preciceData": [
      {
        "mode":                 "read-displacements-velocities",    # mode is one of "read-displacements-velocities", "read-traction", "write-displacements-velocities", "write-traction"
        "preciceMeshName":      "MuscleMeshLeft",                 # name of the precice coupling surface mesh, as given in the precice xml settings file
        "displacementsName":    "Displacement",                     # name of the displacements "data", i.e. field variable, as given in the precice xml settings file
        "velocitiesName":       "Velocity",                         # name of the velocity "data", i.e. field variable, as given in the precice xml settings file
      },
      {
        "mode":                 "write-traction",                   # mode is one of "read-displacements-velocities", "read-traction", "write-displacements-velocities", "write-traction"
        "preciceMeshName":      "MuscleMeshLeft",                 # name of the precice coupling surface mesh, as given in the precice xml settings 
        "tractionName":         "Traction",                         # name of the traction "data", i.e. field variable, as given in the precice xml settings file
      }
    ],
    
    "Coupling": {
      "description":            "fibers and contraction",
      "timeStepWidth":          variables.dt_3D,  # 1e-1
      "logTimeStepWidthAsKey":  "dt_3D",
      "durationLogKey":         "duration_total",
      "timeStepOutputInterval": 1,
      "endTime":                variables.end_time,
      "connectedSlotsTerm1To2": {1:2},          # transfer gamma to MuscleContractionSolver, the receiving slots are λ, λdot, γ
      "connectedSlotsTerm2To1":  None,       # transfer nothing back
      
      "Term1": {        # monodomain, fibers
        "MultipleInstances": {
          "logKey":                     "duration_subdomains_xy",
          "ranksAllComputedInstances":  list(range(n_ranks)),
          "nInstances":                 variables.n_subdomains_xy,
          "instances": 
          [{
            "ranks": list(range(subdomain_coordinate_y*variables.n_subdomains_x + subdomain_coordinate_x, n_ranks, variables.n_subdomains_x*variables.n_subdomains_y)),
            
            # this is for the actual model with fibers
            "StrangSplitting": {
              #"numberTimeSteps": 1,
              "timeStepWidth":          variables.dt_splitting,  # 1e-1
              "logTimeStepWidthAsKey":  "dt_splitting",
              "durationLogKey":         "duration_monodomain",
              "timeStepOutputInterval": 100,
              "endTime":                variables.dt_splitting,
              "connectedSlotsTerm1To2": [0,1,2],   # transfer slot 0 = state Vm from Term1 (CellML) to Term2 (Diffusion)
              "connectedSlotsTerm2To1": [0,None,2],   # transfer the same back, this avoids data copy

              "Term1": {      # CellML, i.e. reaction term of Monodomain equation
                "MultipleInstances": {
                  "logKey":             "duration_subdomains_z",
                  "nInstances":         n_fibers_in_subdomain_x(subdomain_coordinate_x)*n_fibers_in_subdomain_y(subdomain_coordinate_y),
                  "instances": 
                  [{
                    "ranks":                          list(range(variables.n_subdomains_z)),   # these rank nos are local nos to the outer instance of MultipleInstances, i.e. from 0 to number of ranks in z direction
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
                      "additionalSlotNames":          [],

                      "CellML" : {
                        "modelFilename":                          variables.cellml_file,                          # input C++ source file or cellml XML file
                        #"statesInitialValues":                   [],                                             # if given, the initial values for the the states of one instance
                        "statesInitialValues":                    variables.states_initial_values,                # initial values for new_slow_TK
                        "initializeStatesToEquilibrium":          False,                                          # if the equilibrium values of the states should be computed before the simulation starts
                        "initializeStatesToEquilibriumTimestepWidth": 1e-4,                                       # if initializeStatesToEquilibrium is enable, the timestep width to use to solve the equilibrium equation
                        
                        # optimization parameters
                        "optimizationType":                       "vc",                                           # "vc", "simd", "openmp" type of generated optimizated source file
                        "approximateExponentialFunction":         True,                                           # if optimizationType is "vc", whether the exponential function exp(x) should be approximate by (1+x/n)^n with n=1024
                        "compilerFlags":                          "-fPIC -O3 -march=native -shared ",             # compiler flags used to compile the optimized model code
                        "maximumNumberOfThreads":                 0,                                              # if optimizationType is "openmp", the maximum number of threads to use. Default value 0 means no restriction.
                        
                        # stimulation callbacks
                        #"libraryFilename":                       "cellml_simd_lib.so",                           # compiled library
                        #"setSpecificParametersFunction":         set_specific_parameters,                        # callback function that sets parameters like stimulation current
                        #"setSpecificParametersCallInterval":     int(1./variables.stimulation_frequency/variables.dt_0D),         # set_specific_parameters should be called every 0.1, 5e-5 * 1e3 = 5e-2 = 0.05
                        "setSpecificStatesFunction":              set_specific_states,                                             # callback function that sets states like Vm, activation can be implemented by using this method and directly setting Vm values, or by using setParameters/setSpecificParameters
                        #"setSpecificStatesCallInterval":         2*int(1./variables.stimulation_frequency/variables.dt_0D),       # set_specific_states should be called variables.stimulation_frequency times per ms, the factor 2 is needed because every Heun step includes two calls to rhs
                        "setSpecificStatesCallInterval":          0,                                                               # 0 means disabled
                        "setSpecificStatesCallFrequency":         variables.get_specific_states_call_frequency(fiber_no, motor_unit_no),   # set_specific_states should be called variables.stimulation_frequency times per ms
                        "setSpecificStatesFrequencyJitter":       variables.get_specific_states_frequency_jitter(fiber_no, motor_unit_no), # random value to add or substract to setSpecificStatesCallFrequency every stimulation, this is to add random jitter to the frequency
                        "setSpecificStatesRepeatAfterFirstCall":  0.01,                                                            # [ms] simulation time span for which the setSpecificStates callback will be called after a call was triggered
                        "setSpecificStatesCallEnableBegin":       variables.get_specific_states_call_enable_begin(fiber_no, motor_unit_no),# [ms] first time when to call setSpecificStates
                        "additionalArgument":                     fiber_no,                                       # last argument that will be passed to the callback functions set_specific_states, set_specific_parameters, etc.
                        
                        # parameters to the cellml model
                        "mappings":                               variables.mappings,                             # mappings between parameters and algebraics/constants and between outputConnectorSlots and states, algebraics or parameters, they are defined in helper.py
                        "parametersInitialValues":                variables.parameters_initial_values,            #[0.0, 1.0],      # initial values for the parameters: I_Stim, l_hs
                        
                        "meshName":                               "muscle_left_fiber_{}".format(fiber_no),                # reference to the fiber mesh
                        "stimulationLogFilename":                 "out/muscle_left/stimulation.log",                          # a file that will contain the times of stimulations
                      },      
                      "OutputWriter" : [
                        {"format": "Paraview", "outputInterval": 1, "filename": "out/" + variables.scenario_name + "/0D_states({},{})".format(fiber_in_subdomain_coordinate_x,fiber_in_subdomain_coordinate_y), "binary": True, "fixedFormat": False, "combineFiles": True, "fileNumbering": "incremental"}
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
                    "CrankNicolson" : {
                      "initialValues":               [],                                      # no initial values are given
                      #"numberTimeSteps":            1,
                      "timeStepWidth":               variables.dt_1D,                         # timestep width for the diffusion problem
                      "timeStepWidthRelativeTolerance": 1e-10,
                      "logTimeStepWidthAsKey":       "dt_1D",                                 # key under which the time step width will be written to the log file
                      "durationLogKey":              "duration_1D",                           # log key of duration for this solver
                      "timeStepOutputInterval":      1e4,                                     # how often to print the current timestep
                      "dirichletBoundaryConditions": {},                                      # old Dirichlet BC that are not used in FastMonodomainSolver: {0: -75.0036, -1: -75.0036},
                      "dirichletOutputFilename":      None,                                    # filename for a vtp file that contains the Dirichlet boundary condition nodes and their values, set to None to disable
                      "inputMeshIsGlobal":           True,                                    # initial values would be given as global numbers
                      "solverName":                  "diffusionTermSolver",                   # reference to the linear solver
                      "nAdditionalFieldVariables":   2,                                       # number of additional field variables that will be written to the output file, here for stress
                      "checkForNanInf":              True,                                    # abort execution if the solution contains nan or inf values
                      "additionalSlotNames":          [],

                      "FiniteElementMethod" : {
                        "inputMeshIsGlobal":         True,
                        "meshName":                  "muscle_left_fiber_{}".format(fiber_no),
                        "solverName":                "diffusionTermSolver",
                        "prefactor":                 get_diffusion_prefactor(fiber_no, motor_unit_no),  # resolves to Conductivity / (Am * Cm)
                        "slotName":                  ""
                      },
                      "OutputWriter" : [
                        #{"format": "Paraview", "outputInterval": int(1./variables.dt_3D*variables.output_timestep_3D), "filename": "out/"+variables.scenario_name+"/muscle_left_fiber_"+str(fiber_no), "binary": True, "fixedFormat": False, "combineFiles":True, "fileNumbering": "incremental"},

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
                  "OutputWriter" : [
                    {"format": "Paraview", "outputInterval": 1000, "filename": "out/"+variables.scenario_name+"/fibers", "binary": True, "fixedFormat": False, "combineFiles":True, "fileNumbering": "incremental"},
                  ]
                },
              },
            },
            
           
              
          } if (subdomain_coordinate_x,subdomain_coordinate_y) == (variables.own_subdomain_coordinate_x,variables.own_subdomain_coordinate_y) else None
          for subdomain_coordinate_y in range(variables.n_subdomains_y)
              for subdomain_coordinate_x in range(variables.n_subdomains_x)]
        },
        "fiberDistributionFile":    variables.fiber_distribution_file,   # for FastMonodomainSolver, e.g. MU_fibre_distribution_3780.txt
        "firingTimesFile":          variables.firing_times_file,         # for FastMonodomainSolver, e.g. MU_firing_times_real.txt
        "onlyComputeIfHasBeenStimulated": True,                          # only compute fibers after they have been stimulated for the first time
        "disableComputationWhenStatesAreCloseToEquilibrium": True,       # optimization where states that are close to their equilibrium will not be computed again      
        "valueForStimulatedPoint":  variables.vm_value_stimulated,       # to which value of Vm the stimulated node should be set      
      },
      "Term2": {        # solid mechanics
        "MuscleContractionSolver": {
          "numberTimeSteps":              1,                         # only use 1 timestep per interval
          "timeStepOutputInterval":       100,                       # do not output time steps
          "Pmax":                         variables.pmax,            # maximum PK2 active stress
          "slotNames":                    [],                        # names of the data connector slots
          "OutputWriter" : [
            {"format": "Paraview", "outputInterval": 10, "filename": "out/" + variables.scenario_name + "/mechanics_3D", "binary": True, "fixedFormat": False, "onlyNodalValues":True, "combineFiles":True, "fileNumbering": "incremental"},
          ],
          "mapGeometryToMeshes":          ["muscle_left_Mesh"] + [key for key in fiber_meshes.keys()],    # the mesh names of the meshes that will get the geometry transferred
          "reverseMappingOrder":          True,                      # if the mapping target->own mesh should be used instead of own->target mesh. This gives better results in some cases.
          "mapGeometryToMeshes":          [],                        # the mesh names of the meshes that will get the geometry transferred
          "dynamic":                      True,                      # if the dynamic solid mechanics solver should be used, else it computes the quasi-static problem
          "enableForceLengthRelation":    True,
          "lambdaDotScalingFactor":       1,
          # the actual solid mechanics solver, this is either "DynamicHyperelasticitySolver" or "HyperelasticitySolver", depending on the value of "dynamic"
          "DynamicHyperelasticitySolver": {
            "timeStepWidth":              variables.dt_3D,           # time step width 
            "durationLogKey":             "nonlinear",               # key to find duration of this solver in the log file
            "timeStepOutputInterval":     1,                         # how often the current time step should be printed to console
            "totalForceLogFilename":      "",
            "materialParameters":         variables.material_parameters,  # material parameters of the Mooney-Rivlin material
            "density":                    variables.rho,             # density of the material
            "displacementsScalingFactor": 1.0,                       # scaling factor for displacements, only set to sth. other than 1 only to increase visual appearance for very small displacements
            "residualNormLogFilename":    "log_residual_norm.txt",   # log file where residual norm values of the nonlinear solver will be written
            "useAnalyticJacobian":        True,                      # whether to use the analytically computed jacobian matrix in the nonlinear solver (fast)
            "useNumericJacobian":         False,                     # whether to use the numerically computed jacobian matrix in the nonlinear solver (slow), only works with non-nested matrices, if both numeric and analytic are enable, it uses the analytic for the preconditioner and the numeric as normal jacobian
              
            "dumpDenseMatlabVariables":   False,                     # whether to have extra output of matlab vectors, x,r, jacobian matrix (very slow)
            # if useAnalyticJacobian,useNumericJacobian and dumpDenseMatlabVariables all all three true, the analytic and numeric jacobian matrices will get compared to see if there are programming errors for the analytic jacobian
            
            # mesh
            "inputMeshIsGlobal":          True,                     # the mesh is given locally
            "meshName":                   "muscle_left_Mesh_quadratic",        # name of the 3D mesh, it is defined under "Meshes" at the beginning of this config
            "fiberMeshNames":             [],  # fiber meshes that will be used to determine the fiber direction, for multidomain there are no fibers so this would be empty list
            "fiberDirection":             [0,0,1],                  # if fiberMeshNames is empty, directly set the constant fiber direction, in element coordinate system
      
            # solving
            "solverName":                 "mechanicsSolver",         # name of the nonlinear solver configuration, it is defined under "Solvers" at the beginning of this config
            #"loadFactors":                [0.25, 0.66, 1.0],                # load factors for every timestep
            "loadFactors":                [],                        # no load factors, solve problem directly
            "loadFactorGiveUpThreshold":  1e-5,                      # a threshold for the load factor, when to abort the solve of the current time step. The load factors are adjusted automatically if the nonlinear solver diverged. If the load factors get too small, it aborts the solve.
            "nNonlinearSolveCalls":       1,                         # how often the nonlinear solve should be repeated
            
            # boundary and initial conditions
            "dirichletBoundaryConditions": variables.elasticity_dirichlet_bc,   # the initial Dirichlet boundary conditions that define values for displacements u and velocity v
            "neumannBoundaryConditions":   variables.elasticity_neumann_bc,     # Neumann boundary conditions that define traction forces on surfaces of elements
            "divideNeumannBoundaryConditionValuesByTotalArea": True,            # if the given Neumann boundary condition values under "neumannBoundaryConditions" are total forces instead of surface loads and therefore should be scaled by the surface area of all elements where Neumann BC are applied
            "updateDirichletBoundaryConditionsFunction": None,                  # function that updates the dirichlet BCs while the simulation is running
            "updateDirichletBoundaryConditionsFunctionCallInterval": 1,         # every which step the update function should be called, 1 means every time step
            "updateNeumannBoundaryConditionsFunction":   None,                    # function that updates the Neumann BCs while the simulation is running
            "updateNeumannBoundaryConditionsFunctionCallInterval": 1,           # every which step the update function should be called, 1 means every time step

            #TODO: avoid hard coded 45. Smt like mx*my*mz
            "initialValuesDisplacements":  [[0.0,0.0,0.0] for _ in range(nx*ny*nz)],     # the initial values for the displacements, vector of values for every node [[node1-x,y,z], [node2-x,y,z], ...]
            "initialValuesVelocities":     [[0.0,0.0,0.0] for _ in range(nx*ny*nz)],     # the initial values for the velocities, vector of values for every node [[node1-x,y,z], [node2-x,y,z], ...]
            "extrapolateInitialGuess":     True,                                # if the initial values for the dynamic nonlinear problem should be computed by extrapolating the previous displacements and velocities
            "constantBodyForce":           variables.constant_body_force,       # a constant force that acts on the whole body, e.g. for gravity
            
            "dirichletOutputFilename":     "out/"+variables.scenario_name+"/dirichlet_boundary_conditions_muscle",    # filename for a vtp file that contains the Dirichlet boundary condition nodes and their values, set to None to disable
            # "totalForceLogFilename":       "out/muscle_force.csv",              # filename of a log file that will contain the total (bearing) forces and moments at the top and bottom of the volume
            # "totalForceLogOutputInterval": 10,                                  # output interval when to write the totalForceLog file
            # "totalForceBottomElementNosGlobal":  [j*nx + i for j in range(ny) for i in range(nx)],                  # global element nos of the bottom elements used to compute the total forces in the log file totalForceLogFilename
            # "totalForceTopElementNosGlobal":     [(nz-1)*ny*nx + j*nx + i for j in range(ny) for i in range(nx)],   # global element nos of the top elements used to compute the total forces in the log file totalForceTopElementsGlobal
      
            # define which file formats should be written
            # 1. main output writer that writes output files using the quadratic elements function space. Writes displacements, velocities and PK2 stresses.
            "OutputWriter" : [
              
              # Paraview files
              #{"format": "Paraview", "outputInterval": 1, "filename": "out/"+variables.scenario_name+"/u", "binary": True, "fixedFormat": False, "onlyNodalValues":True, "combineFiles":True, "fileNumbering": "incremental"},
              
              # Python callback function "postprocess"
              #{"format": "PythonCallback", "outputInterval": 1, "callback": postprocess, "onlyNodalValues":True, "filename": ""},
            ],
            # 2. additional output writer that writes also the hydrostatic pressure
            "pressure": {   # output files for pressure function space (linear elements), contains pressure values, as well as displacements and velocities
              "OutputWriter" : [
                #{"format": "Paraview", "outputInterval": 1, "filename": "out/"+variables.scenario_name+"/p", "binary": True, "fixedFormat": False, "onlyNodalValues":True, "combineFiles":True, "fileNumbering": "incremental"},
              ]
            },
            # 3. additional output writer that writes virtual work terms
            "dynamic": {    # output of the dynamic solver, has additional virtual work values 
              "OutputWriter" : [   # output files for displacements function space (quadratic elements)
                #{"format": "Paraview", "outputInterval": int(output_interval/dt), "filename": "out/dynamic", "binary": False, "fixedFormat": False, "onlyNodalValues":True, "combineFiles":True, "fileNumbering": "incremental"},
                {"format": "Paraview", "outputInterval": 10, "filename": "out/"+variables.scenario_name+"/virtual_work", "binary": True, "fixedFormat": False, "onlyNodalValues":True, "combineFiles":True, "fileNumbering": "incremental"},
              ],
            },
            # 4. output writer for debugging, outputs files after each load increment, the geometry is not changed but u and v are written
            "LoadIncrements": {   
              "OutputWriter" : [
                #{"format": "Paraview", "outputInterval": 1, "filename": "out/load_increments", "binary": False, "fixedFormat": False, "onlyNodalValues":True, "combineFiles":True, "fileNumbering": "incremental"},
              ]
            },
          }
        }
      }
    }
  }
}



