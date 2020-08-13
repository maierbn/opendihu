# Electrophysiology
# Monodomain with Hodgkin-Huxley model as rhs
#
# This file was created by Nehzat.
# parameters: [<scenario_name>]

import sys

end_time = 5   # [ms] end time of simulation
n_elements = 100 # elements of pysical mesh
n_total = 403 # rows of the snapshot matrix
n_reduced = 99 # number of reduced bases, columns of the left singular vector, is equal to n_reduced+1
snapshots_file = "./out/snapshots.csv"

# global parameters
PMax = 7.3              # maximum stress [N/cm^2]
Conductivity = 3.828    # sigma, conductivity [mS/cm]
Am = 500.0              # surface area to volume ratio [cm^-1]
Cm = 0.58               # membrane capacitance [uF/cm^2]
innervation_zone_width = 0.  # cm
cellml_file = "../../../input/hodgkin_huxley_1952.c"
solver_type = "gmres"

# timing parameters
stimulation_frequency = 10.0      # stimulations per ms
dt_1D = 1e-3                      # timestep width of diffusion
dt_0D = 1e-3                      # timestep width of ODEs
dt_3D = 1e-3                      # overall timestep width of splitting
output_timestep = 1e0            # timestep for output files

# input files
#cellml_file = "../../../input/shorten_ocallaghan_davidson_soboleva_2007.c"
#cellml_file = "../../../input/shorten.cpp"
cellml_file = "../../../input/hodgkin_huxley_1952.c"

#fibre_file = "../../../input/laplace3d_structured_quadratic"
fibre_file = "../../../input/laplace3d_structured_linear"
#fibre_file = "../../../input1000/laplace3d_structured_quadratic"

fibre_distribution_file = "../../../input/MU_fibre_distribution_3780.txt"
#firing_times_file = "../../../input/MU_firing_times_real.txt"
firing_times_file = "../../../input/MU_firing_times_immediately.txt"

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
  parameters_used_as_algebraic = [32]
  parameters_used_as_constant = [65]
  parameters_initial_values = [0.0, 1.0]
  nodal_stimulation_current = 400.
  
elif "hodgkin_huxley" in cellml_file:
  parameters_used_as_algebraic = []
  parameters_used_as_constant = [2]
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
  "logFormat": "csv",
  "Meshes": {
    "MeshFibre": {
      "nElements": n_elements,
      "physicalExtent": n_elements/10.,
      "logKey": "Fiber",
      "inputMeshIsGlobal": True,
    },
    "MeshFibreReduced": {
      "nElements": n_reduced,
      "physicalExtent": n_reduced/10.,
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
    }
  },
  "GodunovSplitting": {
    #"numberTimeSteps": 1,
    "timeStepWidth": dt_3D,  # 1e-1
    "endTime": end_time,
    "logTimeStepWidthAsKey": "dt_3D",
    "durationLogKey": "duration_total",
    "timeStepOutputInterval": 1,
    "connectedSlotsTerm1To2": [0],   # transfer slot 0 = state Vm from Term1 (CellML) to Term2 (Diffusion)
    "connectedSlotsTerm2To1": [0],   # transfer the same back
    "Term1": {      # CellML
      "ModelOrderReduction": {
        "nRowsSnapshots" : n_total,
        "nReducedBases" : n_reduced,
        "snapshots" : snapshots_file,
        "nRowsComponents" : 1,
        "ExplicitEuler" : {
          "timeStepWidth": dt_0D,  # 5e-5
          "initialValues": [],
          "timeStepOutputInterval": 1,
          "logTimeStepWidthAsKey": "dt_0D",
          "durationLogKey": "duration_0D",
          "inputMeshIsGlobal": True,
          "dirichletBoundaryConditions": {},
          "CellML" : {
            "modelFilename": cellml_file,             # input C++ source file, can be either generated by OpenCMISS or OpenCOR from cellml model
            #"simdSourceFilename" : "simdcode.cpp",     # transformed C++ source file that gets generated from modelFilename and is ready for multiple instances
            #"libraryFilename": "cellml_simd_lib.so",   # compiled library
            "useGivenLibrary": False,
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
            "parametersUsedAsAlgebraic": parameters_used_as_algebraic,  #[32],       # list of algebraic value indices, that will be set by parameters. Explicitely defined parameters that will be copied to algebraics, this vector contains the indices of the algebraic array. This is ignored if the input is generated from OpenCMISS generated c code.
            "parametersUsedAsConstant": parameters_used_as_constant,          #[65],           # list of constant value indices, that will be set by parameters. This is ignored if the input is generated from OpenCMISS generated c code.
            "parametersInitialValues": parameters_initial_values,            #[0.0, 1.0],      # initial values for the parameters: I_Stim, l_hs
            "meshName": "MeshFibre",
          },
          
          "OutputWriter" : [
            #{"format": "PythonFile", "outputInterval": int(1./dt_0D*output_timestep), "filename": "out/states", "binary": False, "onlyNodalValues": True},
          ],
        },# ExplicitEuler
        "ExplicitEulerReduced" : {
          "timeStepWidth": dt_0D,  # 5e-5
          "initialValues": [],
          "timeStepOutputInterval": 1,
          "logTimeStepWidthAsKey": "dt_0D",
          "durationLogKey": "duration_0D",
          "inputMeshIsGlobal": True,
          "dirichletBoundaryConditions": {},
          "CellML" : {
            "modelFilename": cellml_file,             # input C++ source file, can be either generated by OpenCMISS or OpenCOR from cellml model
            #"simdSourceFilename" : "simdcode.cpp",     # transformed C++ source file that gets generated from modelFilename and is ready for multiple instances
            #"libraryFilename": "cellml_simd_lib.so",   # compiled library
            "useGivenLibrary": False,
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
          "algebraicsForTransfer": [],
            "parametersUsedAsAlgebraic": parameters_used_as_algebraic,  #[32],       # list of algebraic value indices, that will be set by parameters. Explicitely defined parameters that will be copied to algebraics, this vector contains the indices of the algebraic array. This is ignored if the input is generated from OpenCMISS generated c code.
            "parametersUsedAsConstant": parameters_used_as_constant,          #[65],           # list of constant value indices, that will be set by parameters. This is ignored if the input is generated from OpenCMISS generated c code.
            "parametersInitialValues": parameters_initial_values,            #[0.0, 1.0],      # initial values for the parameters: I_Stim, l_hs
            "meshName": "MeshFibreReduced",
            "prefactor": 1.0,
          },
          
          "OutputWriter" : [
            {"format": "PythonFile", "outputInterval": int(1./dt_0D*output_timestep), "filename": "out/statesReduced", "binary": False, "onlyNodalValues": True},
          ],
        },#ExplicitEulerReduced
      },# ModelOrderReduction
    },
    "Term2": {     # Diffusion
     "ModelOrderReduction": {
      "nRowsSnapshots" : n_total,
      "nReducedBases" : n_reduced,
      "snapshots" : snapshots_file,
      "nRowsComponents" : 1,
      "ImplicitEuler" : {
        "initialValues": [],
        "timeStepWidth": dt_1D,
        "timeStepWidthRelativeTolerance": 1e-10,
        "timeStepOutputInterval": 1,
        "logTimeStepWidthAsKey": "dt_1D",
        "durationLogKey": "duration_1D",
        "inputMeshIsGlobal": True,
        "dirichletBoundaryConditions": {},
        "solverName": "implicitSolver",
        "FiniteElementMethod" : {
          "meshName": "MeshFibre",
          "prefactor": Conductivity/(Am*Cm),
          "solverName": "implicitSolver",
        },
        "OutputWriter" : [
          #{"format": "PythonFile", "outputInterval": int(1./dt_1D*output_timestep), "filename": "out/godunov", "binary": False, "onlyNodalValues": False},
          #{"format": "Paraview", "outputInterval": int(1./dt_1D*output_timestep), "filename": "out/godunov", "binary": False, "fixedFormat": False, "combineFiles": True},
          #{"format": "ExFile", "filename": "out/fibre", "outputInterval": 1e5, "sphereSize": "0.02*0.02*0.02"},
        ],
      },# ExplicitEulerReduced
      "ImplicitEulerReduced" : {
        "initialValues": [],
        "timeStepWidth": dt_1D,
        "timeStepOutputInterval": 1,
        "logTimeStepWidthAsKey": "dt_1D",
        "durationLogKey": "duration_1D",
        "inputMeshIsGlobal": True,
        "dirichletBoundaryConditions": {},
        "solverName": "implicitSolver",
        "FiniteElementMethod" : {
          "meshName": "MeshFibreReduced",
          "prefactor": Conductivity/(Am*Cm),
          "solverName": "implicitSolver",
        },
        "OutputWriter" : [
         {"format": "PythonFile", "outputInterval": int(1./dt_1D*output_timestep), "filename": "out/godunovReduced", "binary": False, "onlyNodalValues": False},
          #{"format": "Paraview", "outputInterval": int(1./dt_1D*output_timestep), "filename": "out/godunov", "binary": False, "fixedFormat": False, "combineFiles": True},
          #{"format": "ExFile", "filename": "out/fibre", "outputInterval": 1e5, "sphereSize": "0.02*0.02*0.02"},
        ],
      },# ExplicitEulerReduced  
    },# ModelOrderReduction
   },
  },
    
  "StrangSplitting": {
    #"numberTimeSteps": 1,
    "timeStepWidth": dt_3D,  # 1e-1
    "endTime": end_time,
    "logTimeStepWidthAsKey": "dt_3D",
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
        "CellML" : {
          "modelFilename": cellml_file,             # input C++ source file, can be either generated by OpenCMISS or OpenCOR from cellml model
          #"simdSourceFilename" : "simdcode.cpp",     # transformed C++ source file that gets generated from modelFilename and is ready for multiple instances
          #"libraryFilename": "cellml_simd_lib.so",   # compiled library
          "useGivenLibrary": False,
          #"statesInitialValues": [],
          #"setParametersFunction": set_parameters,    # callback function that sets parameters like stimulation current
          #"setParametersCallInterval": int(1./stimulation_frequency/dt_0D),     # set_parameters should be called every 0.1, 5e-5 * 1e3 = 5e-2 = 0.05
          #"setSpecificParametersFunction": set_specific_parameters,    # callback function that sets parameters like stimulation current
          #"setSpecificParametersCallInterval": int(1./stimulation_frequency/dt_0D),     # set_specific_parameters should be called every 0.1, 5e-5 * 1e3 = 5e-2 = 0.05
          "setSpecificStatesFunction": set_specific_states,    # callback function that sets states like Vm, activation can be implemented by using this method and directly setting Vm values, or by using setParameters/setSpecificParameters
          "setSpecificStatesCallInterval": int(1./stimulation_frequency/dt_0D),     # set_specific_states should be called every 0.1, 5e-5 * 1e3 = 5e-2 = 0.05
          "additionalArgument": 0,
           #"handleResultFunction": handleResult,
           #"handleResultCallInterval": 2e3,
          
          "statesForTransfer": 0,     # state 0 = Vm
          "algebraicsForTransfer": [],
          "parametersUsedAsAlgebraic": parameters_used_as_algebraic,  #[32],       # list of algebraic value indices, that will be set by parameters. Explicitely defined parameters that will be copied to algebraics, this vector contains the indices of the algebraic array. This is ignored if the input is generated from OpenCMISS generated c code.
          "parametersUsedAsConstant": parameters_used_as_constant,          #[65],           # list of constant value indices, that will be set by parameters. This is ignored if the input is generated from OpenCMISS generated c code.
          "parametersInitialValues": parameters_initial_values,            #[0.0, 1.0],      # initial values for the parameters: I_Stim, l_hs
          "meshName": "MeshFibre",
          "prefactor": 1.0,
        },
        
        "OutputWriter" : [
          #{"format": "PythonFile", "outputInterval": 1e4, "filename": "out/states", "binary": True},
        ],
      },
    },
    "Term2": {     # Diffusion
      "CrankNicolson" : {
        "initialValues": [],
        #"numberTimeSteps": 1,
        "timeStepWidth": dt_1D,
        "timeStepWidthRelativeTolerance": 1e-10,
        "timeStepOutputInterval": 1e4,
        "logTimeStepWidthAsKey": "dt_1D",
        "durationLogKey": "duration_1D",
        "inputMeshIsGlobal": True,
        "dirichletBoundaryConditions": {},
        "solverName": "implicitSolver",
        "FiniteElementMethod" : {
          "meshName": "MeshFibre",
          "prefactor": Conductivity/(Am*Cm),
          "solverName": "implicitSolver",
        },
        "OutputWriter" : [
          {"format": "PythonFile", "outputInterval": int(1./dt_1D*output_timestep), "filename": "out/strang", "binary": True, "onlyNodalValues": False},
          {"format": "Paraview", "outputInterval": int(1./dt_1D*output_timestep), "filename": "out/strang", "binary": True, "fixedFormat": False, "combineFiles": True},
          #{"format": "ExFile", "filename": "out/fibre", "outputInterval": 1e5, "sphereSize": "0.02*0.02*0.02"},
        ],
      },
    },
  }
}
