# 2 fibers, biceps, example for load_balancing
#

end_time = 20.0     # end time for the simulation

import numpy as np
import pickle
import sys

# global parameters
Conductivity = 3.828    # sigma, conductivity [mS/cm]
Am = 500.0              # surface area to volume ratio [cm^-1]
Cm = 0.58               # membrane capacitance [uF/cm^2]
innervation_zone_width = 1.  # cm
innervation_zone_width = 0.  # cm
solver_type = "gmres"   # solver for the linear system

# timing parameters
stimulation_frequency = 10.0      # stimulations per ms
dt_1D = 1e-4                      # timestep width of diffusion
dt_0D = 1e-5                      # timestep width of ODEs
dt_3D = 1e-3                      # overall timestep width of splitting
output_timestep = 1               # timestep for output files

# input files
fiber_file = "../../../input/laplace3d_structured_linear"
fiber_distribution_file = "../../../input/MU_fibre_distribution_3780.txt"
firing_times_file = "../input/MU_firing_times_load_balancing.txt"
cellml_file = "../../../input/hodgkin_huxley_1952.c"

# get own rank no and number of ranks
rank_no = (int)(sys.argv[-2])
n_ranks = (int)(sys.argv[-1])

timestep_adapt_option = "regular"
if len(sys.argv) == 2:
  if rank_no == 0:
    print("usage: mpirun -n 4 ./load_balancing ../settings_load_balancing.py <regular or modified>")
else:
  timestep_adapt_option = sys.argv[0]

if timestep_adapt_option != "regular" and timestep_adapt_option != "modified":
  if rank_no == 0:
    print("timestep_adapt_option has to be either regular or modified!")
  timestep_adapt_option = "regular"

if rank_no == 0:
  print("using timestep_adapt_option={}".format(timestep_adapt_option))
  print("end_time: {}".format(end_time))

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
    # print("rank {}, t: {}, stimulate fiber {} at nodes {}".format(rank_no, current_time, fiber_no, nodes_to_stimulate_global))

    for node_no_global in nodes_to_stimulate_global:
      states[(node_no_global,0,0)] = vm_value_stimulated  # key: ((x,y,z),nodal_dof_index,state_no)

# load fiber meshes, called streamlines
with open(fiber_file, "rb") as f:
  streamlines = pickle.load(f)

# assign fiber meshes to names streamline0, streamline1
if len(streamlines) < 2:
  print("Error: input file {} only contains {} fibers, 2 needed".format(fiber_file, len(streamlines)))
streamline0 = streamlines[0]
#streamline0 = streamlines[0][(int)(len(streamlines[0])/2)-5:(int)(len(streamlines[0])/2)+6]
streamline1 = streamlines[1]
    
# load MU distribution and firing times
fiber_distribution = np.genfromtxt(fiber_distribution_file, delimiter=" ")
firing_times = np.genfromtxt(firing_times_file)

# determine when the fibers will fire, this is only for debugging output and is not used by the config
if rank_no == 0:
  print("Debugging output about fiber firing: Taking input from file \"{}\"".format(firing_times_file))
  
  n_firing_times = np.size(firing_times,0)
  for fiber_no_index in range(2):
    first_stimulation = None
    for current_time in np.linspace(0,1./stimulation_frequency*n_firing_times,n_firing_times):
      if fiber_gets_stimulated(fiber_no_index, stimulation_frequency, current_time):
        first_stimulation = current_time
        break
  
    mu_no = get_motor_unit_no(fiber_no_index)
    print("   fiber {} is of MU {} and will be stimulated for the first time at {}".format(fiber_no_index, mu_no, first_stimulation))

# configure on which ranks fibers 0 and 1 will run
ranks = [
  [0,1],       # rank nos that will compute fiber 0
  [2,3]        # rank nos that will compute fiber 1
]

config = {
  "scenarioName": "load_balancing",
  "solverStructureDiagramFile":     "solver_structure.txt",     # output file of a diagram that shows data connection between solvers
  "mappingsBetweenMeshesLogFile": None,
  "logFormat": "csv",
  # the 1D meshes for each fiber
  "Meshes": {
    "MeshFiber0":
    {
      "nElements": len(streamline0)-1,
      "nodePositions": streamline0,
      "inputMeshIsGlobal": True,
      "setHermiteDerivatives": False
    },
    "MeshFiber1":
    {
      "nElements": len(streamline1)-1,
      "nodePositions": streamline1,
      "inputMeshIsGlobal": True,
      "setHermiteDerivatives": False
    }
  },
  # the linear solver for the diffusion problem
  "Solvers": {
    "implicitSolver": {
      "maxIterations": 1e4,
      "relativeTolerance": 1e-10,
      "absoluteTolerance": 1e-10,         # 1e-10 absolute tolerance of the residual    
      "solverType": solver_type,
      "preconditionerType": "none",
      "dumpFormat": "default",
      "dumpFilename": "",
    }
  },
  # control class that configures multiple instances of the fiber model
  "MultipleInstances": {
    "nInstances": 2,      # number of fibers
    "instances": [        # settings for each fiber, `i` is the index of the fiber (0 or 1)
    {
      "ranks": ranks[i],
      
      # config for strang splitting
      "StrangSplitting": {
        "timeStepWidth": dt_3D,  # 1e-1
        "logTimeStepWidthAsKey": "dt_3D",
        "durationLogKey": "duration_total",
        "timeStepOutputInterval" : 1000,
        "endTime": end_time,
        "connectedSlotsTerm1To2": [0],
        "connectedSlotsTerm2To1": [0],

        "Term1": {      # CellML
          "HeunAdaptive": {
            "timeStepWidth": dt_0D,  # 5e-5
            "tolerance": 1e-4,
            "minTimeStepWidth": 1e-5,
            "timeStepAdaptOption": timestep_adapt_option,
            "lowestMultiplier": 1000,
            "logTimeStepWidthAsKey": "dt_0D",
            "durationLogKey": "duration_0D",
            "initialValues": [],
            "timeStepOutputInterval": 1e4,
            "inputMeshIsGlobal": True,
            "checkForNanInf": True,
            "dirichletBoundaryConditions": {},
            "dirichletOutputFilename":     None,                # filename for a vtp file that contains the Dirichlet boundary condition nodes and their values, set to None to disable
            "nAdditionalFieldVariables": 0,
            "additionalSlotNames": [],
            "timeStepWidthsLogFilename": "out/log_dt.{}.csv".format(rank_no),      # filename for timestep widths, set to None to not capture timestep widths
              
            "CellML" : {
              "modelFilename":                          cellml_file,                                    # input C++ source file or cellml XML file
              "initializeStatesToEquilibrium":          False,                                          # if the equilibrium values of the states should be computed before the simulation starts
              "initializeStatesToEquilibriumTimestepWidth": 1e-4,                                       # if initializeStatesToEquilibrium is enable, the timestep width to use to solve the equilibrium equation
              
              # optimization parameters
              "optimizationType":                       "vc",                                           # "vc", "simd", "openmp" type of generated optimizated source file
              "approximateExponentialFunction":         True,                                           # if optimizationType is "vc", whether the exponential function exp(x) should be approximate by (1+x/n)^n with n=1024
              "compilerFlags":                          "-fPIC -O3 -march=native -shared ",             # compiler flags used to compile the optimized model code
              "maximumNumberOfThreads":                 0,                                              # if optimizationType is "openmp", the maximum number of threads to use. Default value 0 means no restriction.
              
              # stimulation callbacks
              "setSpecificStatesFunction":              set_specific_states,                                             # callback function that sets states like Vm, activation can be implemented by using this method and directly setting Vm values, or by using setParameters/setSpecificParameters
              #"setSpecificStatesCallInterval":         2*int(1./stimulation_frequency/dt_0D),          # set_specific_states should be called variables.stimulation_frequency times per ms, the factor 2 is needed because every Heun step includes two calls to rhs
              "setSpecificStatesCallInterval":          0,                                                               # 0 means disabled
              "setSpecificStatesCallFrequency":         stimulation_frequency,                                           # set_specific_states should be called variables.stimulation_frequency times per ms
              "setSpecificStatesFrequencyJitter":       None,                                                            # random value to add or substract to setSpecificStatesCallFrequency every stimulation, this is to add random jitter to the frequency
              "setSpecificStatesRepeatAfterFirstCall":  0.01,                                                            # [ms] simulation time span for which the setSpecificStates callback will be called after a call was triggered
              "setSpecificStatesCallEnableBegin":       0,                                                               # [ms] first time when to call setSpecificStates
              "additionalArgument":                     i,
              
              "mappings":                               mappings,
              "parametersInitialValues":                parameters_initial_values,
              
              "meshName":                               "MeshFiber"+str(i),
              "stimulationLogFilename":                 "out/stimulation.log",
            },
          },
        },
        
        "Term2": {    # Diffusion
          "ImplicitEuler" : {
            "initialValues": [],
            #"numberTimeSteps": 1,
            "timeStepWidth": dt_1D,  # 1e-5
            "timeStepWidthRelativeTolerance": 1e-10,
            "logTimeStepWidthAsKey": "dt_1D",
            "durationLogKey": "duration_1D",
            "timeStepOutputInterval": 1e4,
            "dirichletBoundaryConditions": {0: -75, -1: -75},      # set first and last value of fiber to -75
            "dirichletOutputFilename":     None,                # filename for a vtp file that contains the Dirichlet boundary condition nodes and their values, set to None to disable
            "timeStepWidthRelativeTolerance": 1e-10,
            "inputMeshIsGlobal": True,
            "solverName": "implicitSolver",
            "checkForNanInf": True,
            "nAdditionalFieldVariables": 0,
            "additionalSlotNames": [],      

            "FiniteElementMethod" : {
              "inputMeshIsGlobal": True,
              "timeStepWidthRelativeTolerance": 1e-10,
              "meshName": "MeshFiber"+str(i),
              "prefactor": Conductivity/(Am*Cm),
              "solverName": "implicitSolver",
              "slotName": "vm",
            },
            "OutputWriter" : [
              {"format": "Paraview", "outputInterval": int(1./dt_1D*output_timestep), "filename": "out/fiber_"+str(i), "binary": True, "fixedFormat": False, "combineFiles": False, "fileNumbering": "incremental"},
              #{"format": "Paraview", "outputInterval": 1./dt_1D*output_timestep, "filename": "out/fiber_"+str(i)+"_txt", "binary": False, "fixedFormat": False, "fileNumbering": "incremental"},
              #{"format": "ExFile", "filename": "out/fiber_"+str(i), "outputInterval": 1./dt_1D*output_timestep, "sphereSize": "0.02*0.02*0.02", "fileNumbering": "incremental"},
              {"format": "PythonFile", "filename": "out/fiber_"+str(i), "outputInterval": int(1./dt_1D*output_timestep), "binary":True, "onlyNodalValues":True, "fileNumbering": "incremental"},
            ]
          },
        },
      }
    }
    for i in range(2)
    ]
  }
}
