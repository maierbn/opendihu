# Electrophysiology
# Monodomain with Shorten model as the rhs. 
#
# parameters: [<scenario_name>]
#
# Plot in the "out" folder:
# plot strang_000*

import sys

end_time = 10.0   # [ms] end time of simulation
n_elements = 100

# global parameters
PMax = 7.3              # maximum stress [N/cm^2]
Conductivity = 3.828    # sigma, conductivity [mS/cm
#Conductivity *= 2
Am = 500.0              # surface area to volume ratio [cm^-1]
Cm = 0.58               # membrane capacitance [uF/cm^2]
innervation_zone_width = 0.  # cm
cellml_file = "../../input/hodgkin_huxley_1952.c"
solver_type = "gmres"

# timing parameters
stimulation_frequency = 10.0      # stimulations per ms
dt_1D = 7e-4                      # timestep width of diffusion
dt_0D = 7e-4                      # timestep width of ODEs
dt_splitting = 7e-4                      # overall timestep width of splitting

#dt_1D = 0.004                      # timestep width of diffusion
#dt_0D = 0.002                     # timestep width of ODEs
#dt_splitting = dt_1D                      # overall timestep width of splitting

output_timestep = 4e-1             # timestep for output files
output_timestep = 1e-1             # timestep for output files

# input files
cellml_file = "../../input/shorten_ocallaghan_davidson_soboleva_2007.c"
cellml_file = "../../input/shorten_ocallaghan_davidson_soboleva_2007.cellml"
#cellml_file = "../../input/shorten.cpp"
#cellml_file = "../../input/hodgkin_huxley_1952.c"

#fiber_file = "../../input/laplace3d_structured_quadratic"
fiber_file = "../../input/laplace3d_structured_linear"
#fiber_file = "../../input1000/laplace3d_structured_quadratic"

fiber_distribution_file = "../../input/MU_fibre_distribution_3780.txt"
firing_times_file = "../../input/MU_firing_times_real.txt"
#firing_times_file = "../../input/MU_firing_times_immediately.txt"

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
  parameters_used_as_intermediate = [32]
  parameters_used_as_constant = [65]
  parameters_initial_values = [0.0, 1.0]
  nodal_stimulation_current = 1200.
  
elif "hodgkin_huxley" in cellml_file:
  parameters_used_as_intermediate = []
  parameters_used_as_constant = [2]
  parameters_initial_values = [0.0]
  nodal_stimulation_current = 40.

# load MU distribution and firing times
fiber_distribution = np.genfromtxt(fiber_distribution_file, delimiter=" ")
firing_times = np.genfromtxt(firing_times_file)

def get_motor_unit_no(fiber_no):
  """
  get the no. of the motor unit which fiber fiber_no is part of
  """
  return int(fiber_distribution[fiber_no % len(fiber_distribution)]-1)

def fiber_gets_stimulated(fiber_no, frequency, current_time):
  """
  determine if fiber fiber_no gets stimulated at simulation time current_time
  """

  # determine motor unit
  alpha = 1.0   # 0.8
  mu_no = (int)(get_motor_unit_no(fiber_no)*alpha)

  # determine if fiber fires now
  index = int(np.round(current_time * frequency))
  n_firing_times = np.size(firing_times,0)

  #if firing_times[index % n_firing_times, mu_no] == 1:
  #  print("{}: fiber {} is mu {}, t = {}, row: {}, stimulated: {} {}".format(rank_no, fiber_no, mu_no, current_time, (index % n_firing_times), firing_times[index % n_firing_times, mu_no], "true" if firing_times[index % n_firing_times, mu_no] == 1 else "false"))

  return firing_times[index % n_firing_times, mu_no] == 1

# callback function that can set states, i.e. prescribed values for stimulation
def set_specific_states(n_nodes_global, time_step_no, current_time, states, fiber_no):

  #print("call set_specific_states at time {}".format(current_time))

  # determine if fiber gets stimulated at the current time
  is_fiber_gets_stimulated = fiber_gets_stimulated(fiber_no, stimulation_frequency, current_time)

  if is_fiber_gets_stimulated:
    # determine nodes to stimulate (center node, left and right neighbour)
    innervation_zone_width_n_nodes = innervation_zone_width*100  # 100 nodes per cm
    innervation_node_global = int(n_nodes_global / 2)  # + np.random.randint(-innervation_zone_width_n_nodes/2,innervation_zone_width_n_nodes/2+1)
    nodes_to_stimulate_global = [innervation_node_global]
    if innervation_node_global > 0:
      nodes_to_stimulate_global.insert(0, innervation_node_global-1)
    if innervation_node_global < n_nodes_global-1:
      nodes_to_stimulate_global.append(innervation_node_global+1)
    print("rank {}, t: {}, stimulate fiber {} at nodes {}".format(rank_no, current_time, fiber_no, nodes_to_stimulate_global))

    for node_no_global in nodes_to_stimulate_global:
      states[(node_no_global,0,0)] = 20.0   # key: ((x,y,z),nodal_dof_index,state_no)

# callback function from output writer
def callback(data, shape, nEntries, dim, timeStepNo, currentTime, null):
  pass
    
config = {
  "scenarioName": scenario_name,
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
      "relativeTolerance": 1e-5,
      "solverType": solver_type,
      "preconditionerType": "none",
      "dumpFormat": "default",
      "dumpFilename": "",
    }
  },
  "GodunovSplitting": {
    #"numberTimeSteps": 1,
    "timeStepWidth": dt_splitting,  # 1e-1
    "endTime": end_time,
    "logTimeStepWidthAsKey": "dt_splitting",
    "durationLogKey": "duration_total",
    "timeStepOutputInterval": 1000,
    "connectedSlotsTerm1To2": [0],   # transfer slot 0 = state Vm from Term1 (CellML) to Term2 (Diffusion)
    "connectedSlotsTerm2To1": [0],   # transfer the same back
    
    "Term1": {      # CellML
      "ExplicitEuler" : {
        "timeStepWidth": dt_0D,  # 5e-5
        "initialValues": [],
        "timeStepOutputInterval": 1e4,
        "logTimeStepWidthAsKey": "dt_0D",
        "durationLogKey": "duration_0D",
        "inputMeshIsGlobal": True,
        "dirichletBoundaryConditions": {},
        "nAdditionalFieldVariables": 0,
        
        "CellML" : {
          "modelFilename": cellml_file,             # input C++ source file, can be either generated by OpenCMISS or OpenCOR from cellml model
          #"statesInitialValues": [],
          #"setParametersFunction": set_parameters,    # callback function that sets parameters like stimulation current
          #"setParametersCallInterval": int(1./stimulation_frequency/dt_0D),     # set_parameters should be called every 0.1, 5e-5 * 1e3 = 5e-2 = 0.05
          #"setSpecificParametersFunction": set_specific_parameters,    # callback function that sets parameters like stimulation current
          #"setSpecificParametersCallInterval": int(1./stimulation_frequency/dt_0D),     # set_parameters should be called every 0.1, 5e-5 * 1e3 = 5e-2 = 0.05
          "setSpecificStatesFunction": set_specific_states,    # callback function that sets states like Vm, activation can be implemented by using this method and directly setting Vm values, or by using setParameters/setSpecificParameters
          "setSpecificStatesCallInterval": int(1./stimulation_frequency/dt_0D),     # set_specific_states should be called every 0.1, 5e-5 * 1e3 = 5e-2 = 0.05
          "additionalArgument": 0,   # use fiber 0
          #"handleResultFunction": handleResult,
          #"handleResultCallInterval": 2e3,
          
          "statesForTransfer": 0,     # state 0 = Vm
          "intermediatesForTransfer": [], # no intermediates are reused in different solvers
          "parametersUsedAsIntermediate": parameters_used_as_intermediate,  #[32],       # list of intermediate value indices, that will be set by parameters. Explicitely defined parameters that will be copied to intermediates, this vector contains the indices of the algebraic array. This is ignored if the input is generated from OpenCMISS generated c code.
          "parametersUsedAsConstant": parameters_used_as_constant,          #[65],           # list of constant value indices, that will be set by parameters. This is ignored if the input is generated from OpenCMISS generated c code.
          "parametersInitialValues": parameters_initial_values,            #[0.0, 1.0],      # initial values for the parameters: I_Stim, l_hs
          "meshName": "MeshFiber",
          "stimulationLogFilename": "out/stimulation.log",
        },
        
        "OutputWriter" : [
          #{"format": "PythonFile", "outputInterval": 1e4, "filename": "out/states", "binary": True, "onlyNodalValues": True},
        ],
      },
    },
    "Term2": {     # Diffusion
      "ExplicitEuler" : {
        "initialValues": [],
        "timeStepWidth": dt_1D,
        "timeStepOutputInterval": 1e4,
        "logTimeStepWidthAsKey": "dt_1D",
        "durationLogKey": "duration_1D",
        "inputMeshIsGlobal": True,
        "dirichletBoundaryConditions": {},
        "solverName": "implicitSolver",
        "nAdditionalFieldVariables": 0,
        
        "FiniteElementMethod" : {
          "meshName": "MeshFiber",
          "prefactor": Conductivity/(Am*Cm),
          "solverName": "implicitSolver",
        },
        "OutputWriter" : [
          {"format": "PythonFile", "outputInterval": int(1./dt_1D*output_timestep), "filename": "out/godunov", "binary": True, "onlyNodalValues": False},
          {"format": "Paraview", "outputInterval": int(1./dt_1D*output_timestep), "filename": "out/godunov", "binary": True, "fixedFormat": False, "combineFiles": True},
          #{"format": "ExFile", "filename": "out/fibe", "outputInterval": 1e5, "sphereSize": "0.02*0.02*0.02"},
        ],
      },
    },
  },
    
  "StrangSplitting": {
    #"numberTimeSteps": 1,
    "timeStepWidth": dt_splitting,  # 1e-1
    "endTime": end_time,
    "logTimeStepWidthAsKey": "dt_splitting",
    "durationLogKey": "duration_total",
    "timeStepOutputInterval": 1000,
    "connectedSlotsTerm1To2": [0],   # transfer slot 0 = state Vm from Term1 (CellML) to Term2 (Diffusion)
    "connectedSlotsTerm2To1": [0],   # transfer the same back
    
    "Term1": {      # CellML
      "Heun" : {
        "timeStepWidth": dt_0D,  # 5e-5
        "initialValues": [],
        "timeStepOutputInterval": 1e4,
        "logTimeStepWidthAsKey": "dt_0D",
        "durationLogKey": "duration_0D",
        "inputMeshIsGlobal": True,
        "dirichletBoundaryConditions": {},
        "nAdditionalFieldVariables": 0,
        
        "CellML" : {
          "modelFilename": "../../input/new_slow_TK_2014_12_08.c", #cellml_file,             # input C++ source file, can be either generated by OpenCMISS or OpenCOR from cellml model
          "optimizationType": "vc",                # "vc", "simd", "openmp" type of generated optimizated source file
          "approximateExponentialFunction": False,
          #"libraryFilename": "cellml_simd_lib.so",   # compiled library
          #"statesInitialValues": [],
          #"setParametersFunction": set_parameters,    # callback function that sets parameters like stimulation current
          #"setParametersCallInterval": int(1./stimulation_frequency/dt_0D),     # set_parameters should be called every 0.1, 5e-5 * 1e3 = 5e-2 = 0.05
          #"setSpecificParametersFunction": set_specific_parameters,    # callback function that sets parameters like stimulation current
          #"setSpecificParametersCallInterval": int(1./stimulation_frequency/dt_0D),     # set_parameters should be called every 0.1, 5e-5 * 1e3 = 5e-2 = 0.05
          "setSpecificStatesFunction": set_specific_states,    # callback function that sets states like Vm, activation can be implemented by using this method and directly setting Vm values, or by using setParameters/setSpecificParameters
          "setSpecificStatesCallInterval": 0,   # int(1./stimulation_frequency/dt_0D),     # set_specific_states should be called every 0.1, 5e-5 * 1e3 = 5e-2 = 0.05
          "setSpecificStatesCallFrequency": stimulation_frequency,     # set_specific_states should be called stimulation_frequency times per ms, the factor 2 is needed because every Heun step includes two calls to rhs
          "setSpecificStatesFrequencyJitter":  0,     # random value to add or substract to setSpecificStatesCallFrequency every stimulation, this is to add random jitter to the frequency
          "setSpecificStatesRepeatAfterFirstCall": 0.01,      # simulation time span for which the setSpecificStates callback will be called after a call was triggered
          "setSpecificStatesCallEnableBegin": 0,         # [ms] first time when to call setSpecificStates
          "additionalArgument": 0,
           #"handleResultFunction": handleResult,
           #"handleResultCallInterval": 2e3,
          "compilerFlags": "-O3 -march=native -fPIC -ftree-vectorize -fopt-info-vec-optimized=vectorizer_optimized.log -shared -lm",
          
          "statesForTransfer": 0,     # state 0 = Vm
          "intermediatesForTransfer": [], # no intermediates are reused in different solvers
          "parametersUsedAsIntermediate": parameters_used_as_intermediate,  #[32],       # list of intermediate value indices, that will be set by parameters. Explicitely defined parameters that will be copied to intermediates, this vector contains the indices of the algebraic array. This is ignored if the input is generated from OpenCMISS generated c code.
          "parametersUsedAsConstant": parameters_used_as_constant,          #[65],           # list of constant value indices, that will be set by parameters. This is ignored if the input is generated from OpenCMISS generated c code.
          "parametersInitialValues": parameters_initial_values,            #[0.0, 1.0],      # initial values for the parameters: I_Stim, l_hs
          "meshName": "MeshFiber",
          "stimulationLogFilename": "out/stimulation.log",
        },
        
        "OutputWriter" : [
          {"format": "PythonFile", "outputInterval": int(1./dt_1D*output_timestep), "filename": "out/states", "binary": True, "onlyNodalValues": True},
        ],
      },
    },
    "Term2": {     # Diffusion
      "CrankNicolson" : {
        "initialValues": [],
        #"numberTimeSteps": 1,
        "timeStepWidth": dt_1D,
        "timeStepOutputInterval": 1e4,
        "logTimeStepWidthAsKey": "dt_1D",
        "durationLogKey": "duration_1D",
        "inputMeshIsGlobal": True,
        "dirichletBoundaryConditions": {},
        "solverName": "implicitSolver",
        "nAdditionalFieldVariables": 0,
        
        "FiniteElementMethod" : {
          "meshName": "MeshFiber",
          "prefactor": Conductivity/(Am*Cm),
          "solverName": "implicitSolver",
        },
        "OutputWriter" : [
          {"format": "PythonFile", "outputInterval": int(1./dt_1D*output_timestep), "filename": "out/strang", "binary": True, "onlyNodalValues": False},
          #{"format": "Paraview", "outputInterval": int(1./dt_1D*output_timestep), "filename": "out/strang", "binary": True, "fixedFormat": False, "combineFiles": True},
          #{"format": "ExFile", "filename": "out/fiber", "outputInterval": 1e5, "sphereSize": "0.02*0.02*0.02"},
        ],
      },
    },
  }
}
