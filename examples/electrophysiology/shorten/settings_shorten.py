# Electrophysiology
# Monodomain with Shorten model as the rhs. 
#
# parameters: [<scenario_name>]
#
# Plot in the "out" folder:
# plot strang_000*

import sys
import numpy as np

end_time = 1000   # [ms] end time of simulation
n_elements = 200

# global parameters
Conductivity = 3.828    # sigma, conductivity [mS/cm]
Am = 500.0              # surface area to volume ratio [cm^-1]
Cm = 0.58               # membrane capacitance [uF/cm^2]
innervation_zone_width = 0.  # [cm]
solver_type = "gmres"

# timing parameters
stimulation_frequency = 10.0      # stimulations per ms
dt_1D = 7e-4                      # timestep width of diffusion
dt_0D = 7e-4                      # timestep width of ODEs
dt_splitting = 7e-4                      # overall timestep width of splitting

dt_0D = 2e-3                        # [ms] timestep width of ODEs
dt_1D = 4e-3                        # [ms] timestep width of diffusion
dt_splitting = 4e-3                 # [ms] overall timestep width of strang splitting

#dt_1D = 0.004                      # timestep width of diffusion
#dt_0D = 0.002                     # timestep width of ODEs
#dt_splitting = dt_1D                      # overall timestep width of splitting

output_timestep = 4e-1             # timestep for output files
#output_timestep = 1e-1             # timestep for output files

# input files
#cellml_file = "../../input/shorten_ocallaghan_davidson_soboleva_2007.c"
#cellml_file = "../../input/shorten_ocallaghan_davidson_soboleva_2007.cellml"
#cellml_file = "../../input/shorten.cpp"
#cellml_file = "../../input/hodgkin_huxley_1952.c"
cellml_file = "../../input/new_slow_TK_2014_12_08.c"

fiber_distribution_file = "../../input/MU_fibre_distribution_3780.txt"
firing_times_file = "../../input/MU_firing_times_real.txt"
#firing_times_file = "../../input/MU_firing_times_immediately.txt"

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
if "shorten" in cellml_file or "TK" in cellml_file:
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
  "solverStructureDiagramFile": "solver_structure.txt",     # filename of file that will contain a visualization of the solver structure and data mapping
  "Meshes": {
    "MeshFiber": {
      "nElements": n_elements,
      "physicalExtent": n_elements/100.,      # 100 elements per cm
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
          "modelFilename":                          cellml_file,                                    # input C++ source file or cellml XML file
          #"statesInitialValues":                   [],                                             # if given, the initial values for the the states of one instance
          "statesInitialValues": [-81.5938, -81.5407, 7.166, 150.916, 6.03857, 12.6618, 131.577, 132.928, 0.00751472, 0.996183, 0.0292367, 0.569413, 0.731601, 0.0075721, 0.996084, 0.0294341, 0.567101, 0.730931, 1.75811e-06, 5.75735e-06, 7.07019e-06, 3.85884e-06, 7.89791e-07, 0.879053, 0.115147, 0.00565615, 0.000123483, 1.01094e-06, -916.582, 0.0284792, 56.5564, 0.0284779, 1687.31, 2.98725, 615, 615, 811, 811, 1342.65, 17807.7, 0.107772, 0.10777, 7243.03, 7243.03, 756.867, 756.867, 956.975, 956.975, 0.0343398, 0.0102587, 0.0136058, 0.0314258, 0.0031226, 0.00249808, 0.223377, 0.264145, 1.74046e-06],
          "initializeStatesToEquilibrium":          False,                                           # if the equilibrium values of the states should be computed before the simulation starts
          "initializeStatesToEquilibriumTimestepWidth": 1e-1,                                       # if initializeStatesToEquilibrium is enable, the timestep width to use to solve the equilibrium equation
          
          # optimization parameters
          "optimizationType":                       "vc",                                           # "vc", "simd", "openmp" type of generated optimizated source file
          "approximateExponentialFunction":         False,                                          # if optimizationType is "vc", whether the exponential function exp(x) should be approximate by (1+x/n)^n with n=1024
          "compilerFlags":                          "-fPIC -O3 -march=native -shared ",             # compiler flags used to compile the optimized model code
          "maximumNumberOfThreads":                 0,                                              # if optimizationType is "openmp", the maximum number of threads to use. Default value 0 means no restriction.
          
          # stimulation callbacks
          #"libraryFilename":                       "cellml_simd_lib.so",                           # compiled library
          #"setParametersFunction":                 set_parameters,                                 # callback function that sets parameters like stimulation current
          #"setParametersCallInterval":             int(1./stimulation_frequency/dt_0D),            # set_parameters should be called every 0.1, 5e-5 * 1e3 = 5e-2 = 0.05
          #"setSpecificParametersFunction":         set_specific_parameters,                        # callback function that sets parameters like stimulation current
          #"setSpecificParametersCallInterval":     int(1./stimulation_frequency/dt_0D),            # set_parameters should be called every 0.1, 5e-5 * 1e3 = 5e-2 = 0.05
          "setSpecificStatesFunction":              set_specific_states,                            # callback function that sets states like Vm, activation can be implemented by using this method and directly setting Vm values, or by using setParameters/setSpecificParameters
          "setSpecificStatesCallInterval":          0,                                              # int(1./stimulation_frequency/dt_0D),     # set_specific_states should be called every 0.1, 5e-5 * 1e3 = 5e-2 = 0.05
          "setSpecificStatesCallFrequency":         stimulation_frequency,                          # set_specific_states should be called stimulation_frequency times per ms, the factor 2 is needed because every Heun step includes two calls to rhs
          "setSpecificStatesFrequencyJitter":       0,                                              # random value to add or substract to setSpecificStatesCallFrequency every stimulation, this is to add random jitter to the frequency
          "setSpecificStatesRepeatAfterFirstCall":  0.01,                                           # simulation time span for which the setSpecificStates callback will be called after a call was triggered
          "setSpecificStatesCallEnableBegin":       0,                                              # [ms] first time when to call setSpecificStates
          "additionalArgument":                     0,                                              # last argument that will be passed to the callback functions set_specific_states, set_specific_parameters, etc.
          #"handleResultFunction": handleResult,
          #"handleResultCallInterval": 2e3,
          
          # output connector slots
          "statesForTransfer": 0,                                                                   # which state values to use in further computation, Shorten / Hodgkin Huxley: state 0 = Vm
          "intermediatesForTransfer": [],                                                           # no intermediates are reused in different solvers
          
          # parameters to the cellml model
          "parametersUsedAsIntermediate":           parameters_used_as_intermediate,  #[32],        # list of intermediate value indices, that will be set by parameters. Explicitely defined parameters that will be copied to intermediates, this vector contains the indices of the algebraic array.
          "parametersUsedAsConstant":               parameters_used_as_constant,      #[65],        # list of constant value indices, that will be set by parameters
          "parametersInitialValues":                parameters_initial_values,        #[0.0, 1.0],  # initial values for the parameters: I_Stim, l_hs
          "meshName": "MeshFiber",                                                                  # name of the mesh to use, this has to be defined under config["Meshes"]
          "stimulationLogFilename": "out/stimulation.log",                                          # a file that will contain the times of stimulations
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
          "inputMeshIsGlobal": True,
        },
        "OutputWriter" : [
          {"format": "PythonFile", "outputInterval": int(1./dt_1D*output_timestep), "filename": "out/vm", "binary": True, "onlyNodalValues": False},
          #{"format": "Paraview", "outputInterval": int(1./dt_1D*output_timestep), "filename": "out/vm", "binary": True, "fixedFormat": False, "combineFiles": True},
          #{"format": "ExFile", "filename": "out/fiber", "outputInterval": 1e5, "sphereSize": "0.02*0.02*0.02"},
        ],
      },
    },
  }
}
