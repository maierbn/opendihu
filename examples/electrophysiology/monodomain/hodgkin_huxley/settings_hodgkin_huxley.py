# Electrophysiology
# Monodomain with Hodgkin-Huxley model as rhs. This also demonstrates how use output writers and how to pass on algebraic values to the diffusion solver.
# The paraview output files "strang_*.vtp" of the output writer under diffusion will include the selected algebraic values no. 0,1 and 2 of the CellML problem. 
# Adjust "algebraicsForTransfer" and "statesForTransfer" to select different states, adjust "connectedSlotsTerm1To2" and connectedSlotsTerm2To1 accordingly.
# There are two more output writers, one directly in the CellML adapter (produces cellml_*.vtp) which writes all states and all algebraics,
# and a second in the Heun solver, which writes all states ("states_*.py)
# 
#
# parameters: [<scenario_name>]

import sys

end_time = 100   # [ms] end time of simulation
n_elements = 100

# global parameters
PMax = 7.3              # maximum stress [N/cm^2]
Conductivity = 3.828    # sigma, conductivity [mS/cm]
Am = 500.0              # surface area to volume ratio [cm^-1]
Cm = 0.58               # membrane capacitance [uF/cm^2]
innervation_zone_width = 0.  # cm
cellml_file = "../../../input/hodgkin_huxley_1952.c"
solver_type = "gmres"

print("prefactor: {}".format(Conductivity/(Am*Cm)))

# timing parameters
stimulation_frequency = 10.0      # stimulations per ms
dt_0D = 1e-3                      # timestep width of ODEs
dt_1D = 1e-3                      # timestep width of diffusion
dt_splitting = 1e-3                      # overall timestep width of splitting

dt_0D = 2e-3                     # timestep width of ODEs
dt_1D = 4e-3                     # timestep width of diffusion
dt_splitting = dt_1D                      # overall timestep width of splitting

output_timestep = 1e0            # timestep for output files

# input files
#cellml_file = "../../../input/shorten_ocallaghan_davidson_soboleva_2007.c"
#cellml_file = "../../../input/shorten.cpp"
cellml_file = "../../../input/hodgkin_huxley_1952.c"

#fibre_file = "../../../input/laplace3d_structured_quadratic"
fibre_file = "../../../input/laplace3d_structured_linear"
#fibre_file = "../../../input1000/laplace3d_structured_quadratic"

fibre_distribution_file = "../../../input/MU_fibre_distribution_3780.txt"
firing_times_file = "../../../input/MU_firing_times_real.txt"
#firing_times_file = "../../../input/MU_firing_times_immediately.txt"

# import needed packages
import sys
import numpy as np

# parse command line options (scenario name)
scenario_name = ""
if len(sys.argv) <= 2:
  scenario_name = ""
else:
  scenario_name = sys.argv[0]

rank_no = (int)(sys.argv[-2])
n_ranks = (int)(sys.argv[-1])

if rank_no == 0:
  print("scenario_name: {}".format(scenario_name))

  print("n elements: {}, end time: {}".format(n_elements,end_time))
  print("prefactor: ",Conductivity/(Am*Cm))

# set values for cellml model
if "shorten" in cellml_file:
  # parameters: stimulation current I_stim, fiber stretch λ
  mappings = {
    ("parameter", 0):           ("algebraic", "wal_environment/I_HH"), # parameter is algebraic 32
    ("parameter", 1):           ("constant", "razumova/L_x"),             # parameter is constant 65, fiber stretch λ, this indicates how much the fiber has stretched, 1 means no extension
    ("outputConnectorSlot", 0): ("state", "wal_environment/vS"),          # expose state

  }
  parameters_initial_values = [0.0, 1.0]
  nodal_stimulation_current = 400.
  
elif "hodgkin_huxley" in cellml_file:
  mappings = {
    ("parameter", 0):           ("constant", "membrane/i_Stim"),      # parameter 0 is constant 2 = I_stim
    ("outputConnectorSlot", 0): ("state", "membrane/V"),              # expose state 0 = Vm to the operator splitting
    ("outputConnectorSlot", 1): ("state", "sodium_channel_m_gate/m"),     # expose state 1 = m
    ("outputConnectorSlot", 2): ("state", "sodium_channel_h_gate/h"),     # expose state 2 = h
    ("outputConnectorSlot", 3): ("state", "potassium_channel_n_gate/n"),  # expose state 3 = n
    ("outputConnectorSlot", 4): ("algebraic", "leakage_current/i_L"),  # expose algebraic 8 = leakage current
  }
  parameters_initial_values = [0.0]
  nodal_stimulation_current = 40.

# load MU distribution and firing times
fibre_distribution = np.genfromtxt(fibre_distribution_file, delimiter=" ")
firing_times = np.genfromtxt(firing_times_file)

def get_motor_unit_no(fibre_no):
  return int(fibre_distribution[fibre_no % len(fibre_distribution)]-1)

def fibre_gets_stimulated(fibre_no, frequency, current_time):

  # determine motor unit
  mu_no = (int)(get_motor_unit_no(fibre_no)*0.8)
  
  # determine if fibre fires now
  index = int(current_time * frequency)
  n_firing_times = np.size(firing_times,0)
  return firing_times[index % n_firing_times, mu_no] == 1
  
# determine if fibre gets stimulation at given time
def fibre_gets_stimulated_constantly(current_time):
  a = current_time * stimulation_frequency
  
  if a - int(a) < 0.1 and current_time < 5:
    return True
  else:
    return False
  
# callback function that can set parameters, i.e. stimulation current
def set_parameters(n_nodes_global, time_step_no, current_time, parameters, dof_nos_global, fibre_no):
  
  # determine if fibre gets stimulated at the current time
  is_fibre_gets_stimulated = fibre_gets_stimulated(fibre_no, stimulation_frequency, current_time)
  
  # determine nodes to stimulate (center node, left and right neighbour)
  innervation_zone_width_n_nodes = innervation_zone_width*100  # 100 nodes per cm
  innervation_node_global = int(n_nodes_global / 2)  # + np.random.randint(-innervation_zone_width_n_nodes/2,innervation_zone_width_n_nodes/2+1)
  nodes_to_stimulate_global = [innervation_node_global]
  
  for k in range(10):
    if innervation_node_global-k >= 0:
      nodes_to_stimulate_global.insert(0, innervation_node_global-k)
    if innervation_node_global+k <= n_nodes_global-1:
      nodes_to_stimulate_global.append(innervation_node_global+k)
  
  # stimulation value
  if is_fibre_gets_stimulated:
    stimulation_current = nodal_stimulation_current
  else:
    stimulation_current = 0.
  
  first_dof_global = dof_nos_global[0]
  last_dof_global = dof_nos_global[-1]
    
  for node_no_global in nodes_to_stimulate_global:
    if first_dof_global <= node_no_global <= last_dof_global:
      # get local no for global no (1D)
      dof_no_local = node_no_global - first_dof_global
      parameters[dof_no_local] = stimulation_current

# callback function that can set parameters, i.e. stimulation current
def set_specific_parameters(n_nodes_global, time_step_no, current_time, parameters, fibre_no):
  
  # determine if fibre gets stimulated at the current time
  is_fibre_gets_stimulated = fibre_gets_stimulated(fibre_no, stimulation_frequency, current_time)
  
  # determine nodes to stimulate (center node, left and right neighbour)
  innervation_zone_width_n_nodes = innervation_zone_width*100  # 100 nodes per cm
  innervation_node_global = int(n_nodes_global / 2)  # + np.random.randint(-innervation_zone_width_n_nodes/2,innervation_zone_width_n_nodes/2+1)
  nodes_to_stimulate_global = [innervation_node_global]
  
  for k in range(10):
    if innervation_node_global-k >= 0:
      nodes_to_stimulate_global.insert(0, innervation_node_global-k)
    if innervation_node_global+k <= n_nodes_global-1:
      nodes_to_stimulate_global.append(innervation_node_global+k)
  
  # stimulation value
  if is_fibre_gets_stimulated:
    stimulation_current = 40.
  else:
    stimulation_current = 0.

  for node_no_global in nodes_to_stimulate_global:
    parameters[(node_no_global,0,0)] = stimulation_current   # key: ((x,y,z),nodal_dof_index,parameter_no)

# callback function that can set states, i.e. prescribed values for stimulation
def set_specific_states(n_nodes_global, time_step_no, current_time, states, fibre_no):
  
  # determine if fibre gets stimulated at the current time
  is_fibre_gets_stimulated = fibre_gets_stimulated(fibre_no, stimulation_frequency, current_time)
  
  if is_fibre_gets_stimulated:  
    # determine nodes to stimulate (center node, left and right neighbour)
    innervation_zone_width_n_nodes = innervation_zone_width*100  # 100 nodes per cm
    innervation_node_global = int(n_nodes_global / 2)  # + np.random.randint(-innervation_zone_width_n_nodes/2,innervation_zone_width_n_nodes/2+1)
    nodes_to_stimulate_global = [innervation_node_global]

    for node_no_global in nodes_to_stimulate_global:
      states[(node_no_global,0,0)] = 20.0   # key: ((x,y,z),nodal_dof_index,state_no)

# callback function from output writer
def callback(data, shape, nEntries, dim, timeStepNo, currentTime, null):
  pass
    
config = {
  "scenarioName": scenario_name,
  "solverStructureDiagramFile":     "solver_structure.txt",     # output file of a diagram that shows data connection between solvers
  "mappingsBetweenMeshesLogFile":   None,
  "Meshes": {
    "MeshFiber": {
      "nElements": n_elements,
      "physicalExtent": n_elements/10.,
      "logKey": "Fiber",
      "inputMeshIsGlobal": True,
    },
  },
  "Solvers": {
    "implicitSolver": {
      "maxIterations": 1e4,
      "relativeTolerance": 1e-10,
      "absoluteTolerance": 1e-10,         # 1e-10 absolute tolerance of the residual    
      "solverType": solver_type,
      "preconditionerType": "none",
      "dumpFormat": "default",
      "dumpFilename": "",   # dump of rhs and system matrix disabled (no filename specified)
    }
  },
  "StrangSplitting": {
    #"numberTimeSteps": 1,
    "timeStepWidth":          dt_splitting,  # 1e-1
    "endTime":                end_time,
    "connectedSlotsTerm1To2": [0,1,2,3,4],   # Transfer slot 0 = state Vm from Term1 (CellML) to Term2 (Diffusion), slots 1-3: algebraics that should only be transferred to Diffusion because of the output writer (such that they will be included in the output files), not for actual computation.
    "connectedSlotsTerm2To1": [0,1,2,3,4],   # Transfer the same values back. Use None for slots that should not be connected. In case of the algebraics it is good to have them connected both directions, 1->2 and 2->1, only then copying will be avoided (variables are reused) because it is asserted that Term 2 does not change the values.
    "logTimeStepWidthAsKey":  "dt_splitting",
    "durationLogKey":         "duration_total",
    "timeStepOutputInterval": 1000,
    "Term1": {      # CellML
      "Heun" : {
        "timeStepWidth":                dt_0D,  # 5e-5
        "initialValues":                [],
        "timeStepOutputInterval":       1e4,
        "logTimeStepWidthAsKey":        "dt_0D",
        "durationLogKey":               "duration_0D",
        "inputMeshIsGlobal":            True,
        "dirichletBoundaryConditions":  {},
        "nAdditionalFieldVariables":    0,
        "checkForNanInf":               True,                                                    
        
        "CellML" : {
          "modelFilename":                          cellml_file,                          # input C++ source file or cellml XML file
          #"statesInitialValues":                   [],                                             # if given, the initial values for the the states of one instance
          "initializeStatesToEquilibrium":          False,                                          # if the equilibrium values of the states should be computed before the simulation starts
          "initializeStatesToEquilibriumTimestepWidth": 1e-4,                                       # if initializeStatesToEquilibrium is enable, the timestep width to use to solve the equilibrium equation
          
          # optimization parameters
          "optimizationType":                       "vc",                                           # "vc", "simd", "openmp" type of generated optimizated source file
          "approximateExponentialFunction":         False,                                          # if optimizationType is "vc", whether the exponential function exp(x) should be approximate by (1+x/n)^n with n=1024
          "compilerFlags":                          "-fPIC -O3 -march=native -shared ",             # compiler flags used to compile the optimized model code
          "maximumNumberOfThreads":                 0,                                              # if optimizationType is "openmp", the maximum number of threads to use. Default value 0 means no restriction.
          
          # stimulation callbacks
          #"libraryFilename":                       "cellml_simd_lib.so",                           # compiled library
          #"setSpecificParametersFunction":         set_specific_parameters,                        # callback function that sets parameters like stimulation current
          #"setSpecificParametersCallInterval":     int(1./stimulation_frequency/dt_0D),            # set_specific_parameters should be called every 0.1, 5e-5 * 1e3 = 5e-2 = 0.05
          "setSpecificStatesFunction":              set_specific_states,                            # callback function that sets states like Vm, activation can be implemented by using this method and directly setting Vm values, or by using setParameters/setSpecificParameters
          "setSpecificStatesCallInterval":          int(1./stimulation_frequency/dt_0D),            # set_specific_states should be called variables.stimulation_frequency times per ms, the factor 2 is needed because every Heun step includes two calls to rhs
          "setSpecificStatesCallFrequency":         0,                                              # set_specific_states should be called variables.stimulation_frequency times per ms
          "setSpecificStatesFrequencyJitter":       0,                                              # random value to add or substract to setSpecificStatesCallFrequency every stimulation, this is to add random jitter to the frequency
          "setSpecificStatesRepeatAfterFirstCall":  0,                                              # [ms] simulation time span for which the setSpecificStates callback will be called after a call was triggered
          "setSpecificStatesCallEnableBegin":       0,                                              # [ms] first time when to call setSpecificStates
          "additionalArgument":                     0,                                              # last argument that will be passed to the callback functions set_specific_states, set_specific_parameters, etc.
          
          # parameters to the cellml model
          "parametersInitialValues":                parameters_initial_values,                      #[0.0, 1.0],      # initial values for the parameters: I_Stim, l_hs
          "mappings":                               mappings,                                       # mappings between parameters and algebraics/constants and between outputConnectorSlots and states, algebraics or parameters, they are defined in helper.py
          
          "meshName":                               "MeshFiber",
          "stimulationLogFilename":                 "out/stimulation.log",                          # a file that will contain the times of stimulations
          
          # output writer for states, algebraics and parameters
          "OutputWriter" : [
            {"format": "Paraview", "outputInterval": int(1./dt_1D*output_timestep), "filename": "out/cellml", "binary": True, "onlyNodalValues": True, "fixedFormat": True, "combineFiles": True, "fileNumbering": "incremental"},
          ],
        },
        
        # output writer only for states
        "OutputWriter" : [
          {"format": "PythonFile", "outputInterval": int(1./dt_1D*output_timestep), "filename": "out/states", "binary": True, "onlyNodalValues": True, "fileNumbering": "incremental"},
        ],
      },
    },
    "Term2": {     # Diffusion
      "ImplicitEuler" : {
        "initialValues": [],
        #"numberTimeSteps": 1,
        "timeStepWidth": dt_1D,
        "timeStepOutputInterval": 1e4,
        "logTimeStepWidthAsKey": "dt_1D",
        "durationLogKey": "duration_1D",
        "inputMeshIsGlobal": True,
        "dirichletBoundaryConditions": {},
        "nAdditionalFieldVariables": 5,
        "solverName": "implicitSolver",
        "checkForNanInf": False,
        
        "FiniteElementMethod" : {
          "meshName": "MeshFiber",
          "prefactor": Conductivity/(Am*Cm),
          "solverName": "implicitSolver",
          "inputMeshIsGlobal": True,
        },
        
        # output writer only for the diffusion variable (i.e. state "Vm")
        "OutputWriter" : [
          {"format": "PythonFile", "outputInterval": int(1./dt_1D*output_timestep), "filename": "out/strang", "binary": True, "onlyNodalValues": False, "fileNumbering": "incremental"},
          {"format": "Paraview",   "outputInterval": int(1./dt_1D*output_timestep), "filename": "out/strang", "binary": True, "fixedFormat": False, "combineFiles": True, "fileNumbering": "incremental"},
          #{"format": "ExFile", "filename": "out/fibre", "outputInterval": 1e5, "sphereSize": "0.02*0.02*0.02"},
        ],
      },
    },
  }
}
