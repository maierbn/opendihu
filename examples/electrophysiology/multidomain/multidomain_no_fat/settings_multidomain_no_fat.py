# Multicompartment 3D, biceps
#

import numpy as np
import scipy.stats
import pickle
import sys,os 
import struct
sys.path.insert(0, "..")

# global parameters
PMax = 7.3              # maximum stress [N/cm^2]
Conductivity = 3.828    # sigma, conductivity [mS/cm]
Am = 500.0              # surface area to volume ratio [cm^-1]
Cm = 0.58               # membrane capacitance [uF/cm^2]

# timing parameters
end_time = 4000.0                   # [ms] end time of the simulation
stimulation_frequency = 100*1e-3    # [ms^-1] sampling frequency of stimuli in firing_times_file, in stimulations per ms, number before 1e-3 factor is in Hertz.
dt_0D = 1e-3                        # [ms] timestep width of ODEs (1e-3)
dt_multidomain = 1e-3               # [ms] timestep width of multidomain solver
dt_splitting = 1e-3                 # [ms] timestep width of splitting
output_timestep = 1e-1              # [ms] timestep for output files

solver_tolerance = 1e-10
#Am = 0.2   # mesh_small
#Am = 0.1

# input files, these are old files, better use "left_biceps_brachii_*"
#fiber_file = "../../../input/scaled_mesh_tiny"
#fiber_file = "../../../input/scaled_mesh_small"
fiber_file = "../../../input/scaled_mesh_normal"
#fiber_file = "../../../input/scaled_mesh_big"

#fiber_file = "../../../input/laplace3d_structured_linear"
#fiber_file = "../../../input/7x7fibers.bin"
fiber_file = "../../../input/left_biceps_brachii_7x7fibers.bin"

# stride which points to select for the 3D mesh, along the muscle (z-direction)
sampling_stride_z = 10


cellml_file = "../../../input/hodgkin_huxley_1952.c"
fiber_distribution_file = "../../../input/MU_fibre_distribution_3780.txt"
#firing_times_file = "../../../input/MU_firing_times_real.txt"
firing_times_file = "../../../input/MU_firing_times_immediately.txt"

# motor unit parameters
motor_units = [
  {"fiber_no": 10, "standard_deviation": 20.0, "maximum": 0.5},
  {"fiber_no": 30, "standard_deviation": 20.0, "maximum": 0.4},
  {"fiber_no": 40, "standard_deviation": 30.0, "maximum": 0.6},
]


# for debugging use the following, non-physiological values. This produces a fast simulation
Am = 1.0
sampling_stride_z = 74
motor_units = motor_units[0:2]    # only 2 motor units
solver_tolerance = 1e-10


n_compartments = len(motor_units)

# own MPI rank no and number of MPI ranks
rank_no = (int)(sys.argv[-2])
n_ranks = (int)(sys.argv[-1])

# load MU distribution and firing times
fiber_distribution = np.genfromtxt(fiber_distribution_file, delimiter=" ")
firing_times = np.genfromtxt(firing_times_file)
  
# load mesh
import helper
(mesh_node_positions,fiber_data,bottom_node_indices,top_node_indices,n_linear_elements_per_coordinate_direction) = helper.load_mesh(fiber_file, sampling_stride_z, rank_no)

# load relative factors for motor units
relative_factors_file = "{}.compartment_relative_factors".format(os.path.basename(fiber_file))

# determine relative factor fields fr(x) for compartments
if os.path.exists(relative_factors_file):
  with open(relative_factors_file, "rb") as f:
    if rank_no == 0:
      print("load relative factors from file \"{}\"".format(relative_factors_file))
    relative_factors = pickle.load(f, encoding='latin1')

else:
  sys.path.append(os.path.abspath(".."))
  n_points_x = n_linear_elements_per_coordinate_direction[0]+1
  n_points_y = n_linear_elements_per_coordinate_direction[1]+1
  n_points_z = n_linear_elements_per_coordinate_direction[2]+1
  n_points_xy = n_points_x * n_points_y
  
  relative_factors = helper.compute_compartment_relative_factors(mesh_node_positions, n_points_xy, n_points_z, fiber_data, motor_units)
  if rank_no == 0:
    print("save relative factors to file \"{}\"".format(relative_factors_file))
    with open(relative_factors_file, "wb") as f:
      pickle.dump(relative_factors, f)

# set variable mappings for cellml model
if "hodgkin_huxley" in cellml_file:
  # parameters: I_stim
  mappings = {
    ("parameter", 0):           ("constant", "membrane/i_Stim"),      # parameter 0 is constant 2 = I_stim
    ("connectorSlot", 0): ("state", "membrane/V"),              # expose state 0 = Vm to the operator splitting
  }
  parameters_initial_values = [0.0]                         # initial value for stimulation current
  nodal_stimulation_current = 40.                           # not used
  vm_value_stimulated = 20.                                 # to which value of Vm the stimulated node should be set (option "valueForStimulatedPoint" of FastMonodomainSolver)
  
  #Cm = 1.0

# for debugging
#n_compartments = 1
#relative_factors = np.ones((n_compartments, len(mesh_node_positions)))   # each row is one compartment
#print("DEBUG settings")

# debugging output
if rank_no == 0:
  for i,factors_list in enumerate(relative_factors.tolist()):
    print("MU {}, maximum fr: {}".format(i,max(factors_list)))

# ---------------------
# callback functions
def get_motor_unit_no(fiber_no):
  return int(fiber_distribution[fiber_no % len(fiber_distribution)]-1)

def compartment_gets_stimulated(compartment_no, current_time):
  # determine motor unit
  mu_no = (int)(get_motor_unit_no(compartment_no)*0.8)
  
  # determine if MU fires now
  index = int(current_time * stimulation_frequency)
  n_firing_times = np.size(firing_times,0)
  return firing_times[index % n_firing_times, mu_no] == 1
  
def set_parameters(n_nodes_global, time_step_no, current_time, parameters, dof_nos_global, compartment_no):
  
  # determine if fiber gets stimulated at the current time
  is_compartment_gets_stimulated = compartment_gets_stimulated(compartment_no, current_time)
  
  # determine nodes to stimulate (center node, left and right neighbour)
  
  n_nodes_x = n_linear_elements_per_coordinate_direction[0]+1
  n_nodes_y = n_linear_elements_per_coordinate_direction[1]+1
  n_nodes_z = n_linear_elements_per_coordinate_direction[2]+1
  z_index_center = (int)(n_nodes_z/2)
  y_index_center = (int)(n_nodes_y/2)
  x_index_center = (int)(n_nodes_x/2)
  
  #nodes_to_stimulate_global = [k*n_nodes_y*n_nodes_x + j*n_nodes_y + i for i in range(n_nodes_x) for j in range(n_nodes_y) for k in [z_index_center-1, z_index_center, z_index_center+1]]

  # stimulation value
  if is_compartment_gets_stimulated:
    stimulation_current = 400.
  else:
    stimulation_current = 0.
  
  for dof_no_local,dof_no_global in enumerate(dof_nos_global):
    k = (int)(dof_no_global / (n_nodes_x*n_nodes_y))
    j = (int)((dof_no_global % (n_nodes_x*n_nodes_y)) / n_nodes_x)
    i = (int)(dof_no_global % n_nodes_x)
  
    if z_index_center-1 <= k <= z_index_center+1:
      if y_index_center-1 <= j <= y_index_center+1:
        if x_index_center-1 <= i <= x_index_center+1:
      
          parameters[dof_no_local] = stimulation_current
  
          if is_compartment_gets_stimulated:
            print("       {}: set stimulation for local dof {}".format(rank_no, dof_no_local))
  
  #print("       {}: setParameters at timestep {}, t={}, n_nodes_global={}, range: [{},{}], fiber no {}, MU {}, stimulated: {}".\
        #format(rank_no, time_step_no, current_time, n_nodes_global, first_dof_global, last_dof_global, fiber_no, get_motor_unit_no(fiber_no), compartment_gets_stimulated))
    
  #wait = input("Press any key to continue...")

# callback function that can set states, i.e. prescribed values for stimulation
def set_specific_states(n_nodes_global, time_step_no, current_time, states, compartment_no):
  
  # determine if fiber gets stimulated at the current time
  is_compartment_gets_stimulated = compartment_gets_stimulated(compartment_no, current_time)

  if is_compartment_gets_stimulated:  
    
    n_nodes_x = n_linear_elements_per_coordinate_direction[0]+1
    n_nodes_y = n_linear_elements_per_coordinate_direction[1]+1
    n_nodes_z = n_linear_elements_per_coordinate_direction[2]+1
    z_index_center = (int)(n_nodes_z/2)
    y_index_center = (int)(n_nodes_y/2)
    x_index_center = (int)(n_nodes_x/2)
    
    for k in range(n_nodes_z):
      if z_index_center-1 <= k <= z_index_center+1:
        for j in range(n_nodes_y):
          if y_index_center-1 <= j <= y_index_center+1:
            for i in range(n_nodes_x):
              if x_index_center-1 <= i <= x_index_center+1:
                key = ((i,j,k),0,0)        # key: ((x,y,z),nodal_dof_index,state_no)
                states[key] = 20.0
                #print("set states at ({},{},{}) to 40".format(i,j,k))

    #print("states: {}".format(states))
    #print("n_nodes: ({},{},{})".format(n_nodes_x, n_nodes_y, n_nodes_z))
    #print("n_nodes_global: {}, time_step_no: {}, current_time: {}, compartment_no: {}".format(n_nodes_global, time_step_no, current_time, compartment_no))
    #wait = input("Press any key to continue...")
    
# boundary conditions for potential flow
potential_flow_bc = {}
for bottom_node_index in bottom_node_indices:
  potential_flow_bc[bottom_node_index] = 0.0
  
for top_node_index in top_node_indices:
  potential_flow_bc[top_node_index] = 1.0
  
# settings for the multidomain solver
multidomain_solver = {
  "nCompartments":                    n_compartments,                     # number of compartments
  "am":                               Am,                                 # Am parameter (ration of surface to volume of fibers)
  "cm":                               Cm,                                 # Cm parameter (capacitance of the cellular membrane)
  "timeStepWidth":                    dt_multidomain,                     # time step width of the diffusion, i.e. the global linear system in the multidomain solver
  "endTime":                          end_time,                           # end time, this is not relevant because it will be overridden by the splitting scheme
  "timeStepOutputInterval":           100,                                # how often the output timestep should be printed
  "solverName":                       "activationSolver",                 # reference to the solver used for the global linear system of the multidomain eq.
  "slotNames":                        [],
  "initialGuessNonzero":              True,                               # if the initial guess for the 3D system should be set as the solution of the previous timestep, this only makes sense for iterative solvers
  "inputIsGlobal":                    True,                               # if values and dofs correspond to the global numbering
  "showLinearSolverOutput":           False,                              # if convergence information of the linear solver in every timestep should be printed, this is a lot of output for fast computations
  "compartmentRelativeFactors":       relative_factors.tolist(),          # list of lists of the factors for every dof, because "inputIsGlobal": True, this contains the global dofs
  "updateSystemMatrixEveryTimestep":  False,                              # if this multidomain solver will update the system matrix in every first timestep, us this only if the geometry changed, e.g. by contraction
  "useSymmetricPreconditionerMatrix": True,                               # if the diagonal blocks of the system matrix should be used as preconditioner matrix
  
  "PotentialFlow": {
    "FiniteElementMethod" : {  
      "meshName":                     "mesh",
      "solverName":                   "potentialFlowSolver",
      "prefactor":                    1.0,
      "slotName":                     "",
      "dirichletBoundaryConditions":  potential_flow_bc,
      "dirichletOutputFilename":      None,               # output filename for the dirichlet boundary conditions, set to "" to have no output
      "neumannBoundaryConditions":    [],
      "inputMeshIsGlobal":            True,
    },
  },
  "Activation": {
    "FiniteElementMethod" : {  
      "meshName":                     "mesh",
      "solverName":                   "activationSolver",
      "prefactor":                    1.0,
      "slotName":                     "",
      "inputMeshIsGlobal":            True,
      "dirichletBoundaryConditions":  {},
      "dirichletOutputFilename":      None,               # output filename for the dirichlet boundary conditions, set to "" to have no output
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
  
  "OutputWriter" : [
    {"format": "Paraview", "outputInterval": (int)(1./dt_multidomain*output_timestep), "filename": "out/output", "binary": True, "fixedFormat": False, "combineFiles": True, "fileNumbering": "incremental"},
    #{"format": "ExFile", "filename": "out/fiber_"+str(i), "outputInterval": 1./dt_1D*output_timestep, "sphereSize": "0.02*0.02*0.02"},
    #{"format": "PythonFile", "filename": "out/fiber_"+str(i), "outputInterval": int(1./dt_1D*output_timestep), "binary":True, "onlyNodalValues":True},
  ]
}
  
config = {
  "scenarioName":              "multidomain_no_fat",        # name of the scenario such that the run can be identified in the log file
  "logFormat":                  "csv",                      # "csv" or "json", format of the lines in the log file, csv gives smaller files
  "solverStructureDiagramFile": "solver_structure.txt",     # output file of a diagram that shows data connection between solvers
  "mappingsBetweenMeshesLogFile": "out/mappings_between_meshes.txt",  # name of a log file for the mapping between meshes
  "Meshes": {
    "mesh": {
      "nElements":             n_linear_elements_per_coordinate_direction,
      "nodePositions":         mesh_node_positions,
      "inputMeshIsGlobal":     True,
      "setHermiteDerivatives": False
    }
  },
  "Solvers": {
    "potentialFlowSolver": {
      "relativeTolerance":  1e-10,
      "absoluteTolerance":  1e-10,         # 1e-10 absolute tolerance of the residual          
      "maxIterations":      1e5,
      "solverType":         "gmres",
      "preconditionerType": "none",
      "dumpFormat":         "default",
      "dumpFilename":       "",
    },
    "activationSolver": {
      "relativeTolerance":  1e-9,
      "absoluteTolerance":  solver_tolerance,         # 1e-10 absolute tolerance of the residual          
      "maxIterations":      1e3,
      "solverType":         "gmres",
      "preconditionerType": "ilu" if n_ranks == 1 else "euclid",
      "hypreOptions":       "",              # additional command line type options of petsc for the hypre solver
      "dumpFormat":         "matlab",
      "dumpFilename":       "", #"out/no",
    }
  },
  "StrangSplitting": {
    "timeStepWidth":          dt_splitting,  # 1e-1
    "logTimeStepWidthAsKey":  "dt_splitting",
    "durationLogKey":         "duration_total",
    "timeStepOutputInterval": 100,
    "endTime":                end_time,
    "connectedSlotsTerm1To2": [0],          # CellML V_mk (0) <=> Multidomain V_mk^(i) (0)
    "connectedSlotsTerm2To1": [None, 0],    # Multidomain V_mk^(i+1) (1) -> CellML V_mk (0)

    "Term1": {      # CellML
      "MultipleInstances": {
        "nInstances": n_compartments,  
        "instances": [        # settings for each motor unit, `i` is the index of the motor unit
        {
          "ranks": list(range(n_ranks)),
          "Heun" : {
            "timeStepWidth": dt_0D,  # 5e-5
            "logTimeStepWidthAsKey":        "dt_0D",
            "durationLogKey":               "duration_0D",
            "initialValues":                [],
            "timeStepOutputInterval":       1e4,
            "inputMeshIsGlobal":            True,
            "dirichletBoundaryConditions":  {},
            "dirichletOutputFilename":      None,               # output filename for the dirichlet boundary conditions, set to "" to have no output
            "nAdditionalFieldVariables":    0,
            "additionalSlotNames":          [],
            "checkForNanInf":               True,
                
            "CellML" : {
              "modelFilename":                          cellml_file,                            # input C++ source file or cellml XML file
              "initializeStatesToEquilibrium":          False,                                  # if the equilibrium values of the states should be computed before the simulation starts
              "initializeStatesToEquilibriumTimestepWidth": 1e-4,                               # if initializeStatesToEquilibrium is enable, the timestep width to use to solve the equilibrium equation
              
              # optimization parameters
              "optimizationType":                       "vc",                                   # "vc", "simd", "openmp" type of generated optimizated source file
              "approximateExponentialFunction":         True,                                   # if optimizationType is "vc", whether the exponential function exp(x) should be approximate by (1+x/n)^n with n=1024
              "compilerFlags":                          "-fPIC -O3 -march=native -shared ",     # compiler flags used to compile the optimized model code
              "maximumNumberOfThreads":                 0,                                      # if optimizationType is "openmp", the maximum number of threads to use. Default value 0 means no restriction.
              
              # stimulation callbacks
              #"statesInitialValues": [],
              "setSpecificStatesFunction":              set_specific_states,                    # callback function that sets states like Vm, activation can be implemented by using this method and directly setting Vm values, or by using setParameters/setSpecificParameters
              "setSpecificStatesCallInterval":          int(1./stimulation_frequency/dt_0D),    # set_specific_states should be called every 0.1, 5e-5 * 1e3 = 5e-2 = 0.05
              "setSpecificStatesCallFrequency":         0,                                      # set_specific_states should be called  stimulation_frequency times per ms
              "setSpecificStatesFrequencyJitter":       0,                                      # random value to add or substract to setSpecificStatesCallFrequency every stimulation, this is to add random jitter to the frequency
              "setSpecificStatesRepeatAfterFirstCall":  0.01,                                   # [ms] simulation time span for which the setSpecificStates callback will be called after a call was triggered
              "setSpecificStatesCallEnableBegin":       0,                                      # [ms] first time when to call setSpecificStates
              "additionalArgument":                     compartment_no,                         # the compartment no is the last argument to set_specific_states function such that it knows which compartments to stimulate when
              #"setParametersFunction":                 set_parameters,                         # callback function that sets parameters like stimulation current
              #"setParametersCallInterval":             int(1./stimulation_frequency/dt_0D),    # set_parameters should be called every 0.1, 5e-5 * 1e3 = 5e-2 = 0.05
              #"setParametersFunctionAdditionalParameter": compartment_no,

              # parameters to the cellml model
              "mappings":                               mappings,
              "parametersInitialValues":                parameters_initial_values,              #[0.0, 1.0],      # initial values for the parameters: I_Stim, l_hs
              
              "meshName": "mesh",
              "stimulationLogFilename": "out/stimulation.log",
            }
          }
        } for compartment_no in range(n_compartments)]
      },
    },
    "Term2": {     # Diffusion, i.e. Multidomain
      "MultidomainSolver" : multidomain_solver,
      "OutputSurface": {        # version for fibers_emg_2d_output
        "OutputWriter": [
          {"format": "Paraview", "outputInterval": (int)(1./dt_multidomain*output_timestep), "filename": "out/surface", "binary": True, "fixedFormat": False, "combineFiles": True, "fileNumbering": "incremental"},
        ],
        "face": "1-",
        "MultidomainSolver" : multidomain_solver,
      }
    }
  }
}

print("Linear solver type: {}".format(config["Solvers"]["activationSolver"]["solverType"]))
