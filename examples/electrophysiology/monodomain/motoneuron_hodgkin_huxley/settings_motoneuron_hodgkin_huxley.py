# Electrophysiology
# Monodomain with Hodgkin-Huxley model as rhs, stimulated by a Hodgkin-Huxley cell 
# as motoneuron with constant stimulation current that fires when Vm>0.

import sys
import numpy as np

end_time = 100   # [ms] end time of simulation
n_elements = 200
element_size = 1./100   # [cm]
#element_size = 1./10

# global parameters
Conductivity = 3.828    # sigma, conductivity [mS/cm]
Am = 500.0              # surface area to volume ratio [cm^-1]
Cm = 0.58               # membrane capacitance [uF/cm^2]
innervation_zone_width = 0.  # [cm]
solver_type = "gmres"

# timing parameters
dt_motoneuron = 1e-3             # [ms] timestep width for motoneurons 
dt_stimulation_check = 1e-2      # [ms] timestep width for when to check if the motoneurons stimulated 
dt_0D = 2e-3                     # [ms] timestep width of ODEs
dt_1D = 4e-3                     # [ms] timestep width of diffusion
dt_splitting = max(dt_0D,dt_1D)  # [ms] overall timestep width of splitting

output_timestep = 1e0            # [ms] timestep for output files
output_timestep_motoneuron = 2e-1  # [ms] timestep for output files of the motoneuron

# input files
motoneuron_cellml_file = "../../../input/motoneuron_hodgkin_huxley.cellml"
monodomain_cellml_file = "../../../input/hodgkin_huxley_1952.cellml"

# parse command line options (scenario name)
scenario_name = ""

motoneuron_stimulation_current = 5    # [mV] constant stimulation current of the motoneuron

# define the variable mappings for the motoneuron model
motoneuron_mappings = {
  ("parameter", 0):           "membrane/i_Stim",
  ("parameter", 1):           "firing_threshold/V_threshold",
  ("parameter", 2):           "firing_threshold/V_firing",
  ("parameter", 3):           "firing_threshold/V_extern_in",
  
  ("connectorSlot", 0, "v_in"): "firing_threshold/V_extern_in",
  ("connectorSlot", 1, "v_out"): "firing_threshold/V_extern_out",
}

# set values for parameters: [i_Stim, V_threshold, V_firing, V_extern_in], 
# i_Stim = stimulation current of motoneuron, V_threshold = threshold if trans-membrane voltage is above, motoneuron fires, V_firing = value of Vm to set if motoneuron fires
motoneuron_parameters_initial_values = [motoneuron_stimulation_current, 0, 20, -75]

# parse number of ranks
rank_no = (int)(sys.argv[-2])
n_ranks = (int)(sys.argv[-1])

# output some parameters
if rank_no == 0:
  print("scenario_name: {}".format(scenario_name))

  print("n elements: {}, end time: {}".format(n_elements,end_time))
  print("prefactor: ",Conductivity/(Am*Cm))
  print("dt_0D: {}, dt_1D: {}, dt_splitting: {}".format(dt_0D, dt_1D, dt_splitting))
  
# define the variable mappings for the Monodomain subcellular model
if "hodgkin_huxley" in monodomain_cellml_file:
  mappings = {
    ("parameter", 0):        ("constant", "membrane/i_Stim"),      # parameter 0 is constant 2 = I_stim
    ("connectorSlot", "vm"): ("state", "membrane/V"),              # expose state 0 = Vm to the operator splitting
  }
  parameters_initial_values = [0.0]
    
config = {
  "scenarioName":                 scenario_name,
  "logFormat":                    "csv",     # "csv" or "json", format of the lines in the log file, csv gives smaller files
  "solverStructureDiagramFile":   "solver_structure.txt",     # filename of file that will contain a visualization of the solver structure and data mapping
  "mappingsBetweenMeshesLogFile": "mappings_between_meshes_log.txt",    # log file for mappings 
  
  "Meshes": {
    "MeshFiber": {
      "nElements":          n_elements,
      "physicalExtent":     n_elements*element_size,      # 100 elements per cm
      "physicalOffset":     [0,0,0],
      "logKey":             "Fiber",
      "inputMeshIsGlobal":  True,
    },
    "motoneuronMesh": {
      "nElements" :         0 if n_ranks == 1 else n_ranks,
      "physicalExtent":     1,
      "physicalOffset":     0,
      "logKey":             "motoneuron",
      "inputMeshIsGlobal":  True,
      "nRanks":             n_ranks
    }
  },
  
  "Solvers": {
    "implicitSolver": {
      "maxIterations":      1e4,
      "relativeTolerance":  1e-5,
      "absoluteTolerance":  1e-10,         # 1e-10 absolute tolerance of the residual          
      "solverType":         solver_type,
      "preconditionerType": "none",
      "dumpFormat":         "default",
      "dumpFilename":       "",   # dump of rhs and system matrix disabled (no filename specified)
    }
  },
  "Coupling": {
    "endTime":                end_time,
    "timeStepWidth":          dt_stimulation_check,
    "logTimeStepWidthAsKey":  "dt_stimulation_check",
    "durationLogKey":         "duration_multidomain",
    "timeStepOutputInterval": 100,
    "connectedSlotsTerm1To2": None,      # slot connections between motoneuron and MapDofs object which maps values to nodes of the fiber mesh 
    "connectedSlotsTerm2To1": None,
        
    # motoneuron
    "Term1": {    
      "Heun" : {
        "timeStepWidth":                dt_motoneuron,
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
            
        "CellML" : {
          "modelFilename":                          motoneuron_cellml_file,                          # input C++ source file or cellml XML file
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
          
          "mappings":                               motoneuron_mappings,                             # mappings between parameters and algebraics/constants and between outputConnectorSlots and states, algebraics or parameters, they are defined in helper.py
          "parametersInitialValues":                motoneuron_parameters_initial_values,            #[0.0, 1.0],      # initial values for the parameters: I_Stim, l_hs
          
          "meshName":                               "motoneuronMesh",                                       # use the linear mesh, it was partitioned by the helper.py script which called opendihu/scripts/create_partitioned_meshes_for_settings.py
          "stimulationLogFilename":                 "out/stimulation.log",

          # output writer for states, algebraics and parameters                
          "OutputWriter" : [
            {"format": "PythonFile", "outputInterval": (int)(2./dt_motoneuron*output_timestep_motoneuron), "filename": "out/motoneuron", "binary": True, "fixedFormat": False, "combineFiles": True, "onlyNodalValues": True, "fileNumbering": "incremental"}
          ]
        },
      }
    },
    
    # Monodomain
    "Term2": {
      "MapDofs": {
        "nAdditionalFieldVariables":  2,                              # number of additional field variables that are defined by this object. They have 1 component, use the templated function space and mesh given by meshName.
        "additionalSlotNames":        ["v_out", "v_in"],
        "meshName":                   "motoneuronMesh",               # the mesh on which the additional field variables will be defined
        
        # mapping from motoneuronMesh which contains on every rank as many nodes as there are motoneurons to the 3D domain
        # map from motoneuronMesh (algebraics) to 3Dmesh (solution)
        "beforeComputation": [                                        # transfer/mapping of dofs that will be performed before the computation of the nested solver
          {                                                 
            "fromConnectorSlot":                2,                    # which fiber/compartment
            "toConnectorSlots":                 0,
            "fromSlotConnectorArrayIndex":      0,
            "toSlotConnectorArrayIndex":        0,
            "fromDofNosNumbering":              "local",
            "toDofNosNumbering":                "global",
            "dofsMapping":                      {0: [n_elements/2-1, n_elements/2, n_elements/2+1]},    # map from motoneuron 0 to 3 center elements of fiber
            "mode":                             "copyLocal",          # "copyLocal", "copyLocalIfPositive" or "communicate"
          }
        ],
        
        # map from 3Dmesh to motoneuronMesh
        "afterComputation": [                                         # transfer/mapping of dofs that will be performed before the computation of the nested solver
          {                                                 
            "fromConnectorSlot":                0,
            "toConnectorSlots":                 3,
            "fromSlotConnectorArrayIndex":      0,                    # which fiber/compartment
            "toSlotConnectorArrayIndex":        0,
            "fromDofNosNumbering":              "global",
            "toDofNosNumbering":                "local",
            "dofsMapping":                      {n_elements/2: 0},    # map from center elements of fiber to motoneuron
            "mode":                             "copyLocal",          # "copyLocal", "copyLocalIfPositive" or "communicate"
          }
        ],
          
        "StrangSplitting": {
          "timeStepWidth":              dt_splitting,  # 1e-1
          "logTimeStepWidthAsKey":      "dt_splitting",
          "durationLogKey":             "duration_total",
          "timeStepOutputInterval":     1000,
          "connectedSlotsTerm1To2":     None,   # connection is handled by slot names
          "connectedSlotsTerm2To1":     None,  
          
          "Term1": {      # CellML
            "Heun" : {
              "timeStepWidth":                dt_0D,  # 5e-5
              "initialValues":                [],
              "timeStepOutputInterval":       1e4,
              "logTimeStepWidthAsKey":        "dt_0D",
              "durationLogKey":               "duration_0D",
              "inputMeshIsGlobal":            True,
              "dirichletBoundaryConditions":  {},
              "dirichletOutputFilename":      None,                                             # filename for a vtp file that contains the Dirichlet boundary condition nodes and their values, set to None to disable
              "nAdditionalFieldVariables":    0,
              "additionalSlotNames":          [],
              "checkForNanInf":               True,                                             # check if the solution vector contains nan or +/-inf values, if yes, an error is printed. This is a time-consuming check.                                
              
              "CellML" : {
                "modelFilename":                          monodomain_cellml_file,                         # input C++ source file or cellml XML file
                "statesInitialValues":                    [],                                             # if given, the initial values for the the states of one instance
                "initializeStatesToEquilibrium":          False,                                          # if the equilibrium values of the states should be computed before the simulation starts
                "initializeStatesToEquilibriumTimestepWidth": 1e-4,                                       # if initializeStatesToEquilibrium is enable, the timestep width to use to solve the equilibrium equation
                
                # optimization parameters
                "optimizationType":                       "vc",                                           # "vc", "simd", "openmp" type of generated optimizated source file
                "approximateExponentialFunction":         True,                                           # if optimizationType is "vc", whether the exponential function exp(x) should be approximate by (1+x/n)^n with n=1024
                "compilerFlags":                          "-fPIC -O3 -march=native -shared ",             # compiler flags used to compile the optimized model code
                "maximumNumberOfThreads":                 0,                                              # if optimizationType is "openmp", the maximum number of threads to use. Default value 0 means no restriction.
                
                # stimulation callbacks
                #"libraryFilename":                       "cellml_simd_lib.so",                           # compiled library
                
                # stimulation is done by motoneuron, no need for callbacks
                # (stimulation by setting I_Stim)
                "setSpecificParametersFunction":          None,                                           # callback function that sets parameters like stimulation current
                "setSpecificParametersCallInterval":      0,                                              # set_parameters should be called every 0.1, 5e-5 * 1e3 = 5e-2 = 0.05
                
                # (stimulation by setting Vm to a prescribed value)
                "setSpecificStatesFunction":              None,                                           # callback function that sets states like Vm, activation can be implemented by using this method and directly setting Vm values, or by using setParameters/setSpecificParameters
                "setSpecificStatesCallInterval":          0,                                              # int(1./stimulation_frequency/dt_0D),     # set_specific_states should be called every 0.1, 5e-5 * 1e3 = 5e-2 = 0.05
                "setSpecificStatesCallFrequency":         0,                                              # set_specific_states should be called stimulation_frequency times per ms, the factor 2 is needed because every Heun step includes two calls to rhs
                "setSpecificStatesFrequencyJitter":       0,                                              # random value to add or substract to setSpecificStatesCallFrequency every stimulation, this is to add random jitter to the frequency
                "setSpecificStatesRepeatAfterFirstCall":  0.01,                                           # simulation time span for which the setSpecificStates callback will be called after a call was triggered
                "setSpecificStatesCallEnableBegin":       0,                                              # [ms] first time when to call setSpecificStates
                
                "additionalArgument":                     0,                                              # last argument that will be passed to the callback functions set_specific_states, set_specific_parameters, etc.
                #"handleResultFunction": handleResult,
                #"handleResultCallInterval": 2e3,
                
                # parameters to the cellml model
                "parametersInitialValues":                parameters_initial_values,                      #[0.0, 1.0],      # initial values for the parameters: I_Stim, l_hs
                "mappings":                               mappings,                                       # mappings between parameters and algebraics/constants and between outputConnectorSlots and states, algebraics or parameters, they are defined in helper.py
                
                "meshName":                               "MeshFiber",
                "stimulationLogFilename":                 "out/stimulation.log",                          # a file that will contain the times of stimulations
                
                # output writer for states, algebraics and parameters
                "OutputWriter" : [
      #            {"format": "Paraview", "outputInterval": int(1./dt_0D*output_timestep), "filename": "out/cellml", "binary": True, "onlyNodalValues": True, "fixedFormat": True, "combineFiles": True, "fileNumbering": "incremental"},
                  {"format": "PythonFile", "outputInterval": int(2./dt_0D*output_timestep), "filename": "out/cellml", "binary": True, "onlyNodalValues": True, "fixedFormat": True, "combineFiles": True, "fileNumbering": "incremental"},
                ],
              },
              
              # output writer only for states
              "OutputWriter" : [
                #{"format": "PythonFile", "outputInterval": int(1./dt_0D*output_timestep), "filename": "out/states", "binary": True, "onlyNodalValues": True, "fileNumbering": "incremental"},
              ],
            },
          },
          "Term2": {     # Diffusion
            "CrankNicolson" : {
              "initialValues": [],
              #"numberTimeSteps": 1,
              "timeStepWidth":                dt_1D,
              "timeStepWidthRelativeTolerance": 1e-10,
              "timeStepOutputInterval":       1e4,
              "logTimeStepWidthAsKey":        "dt_1D",
              "durationLogKey":               "duration_1D",
              "inputMeshIsGlobal":            True,
              "dirichletBoundaryConditions":  {},
              "dirichletOutputFilename":      None,             # filename for a vtp file that contains the Dirichlet boundary condition nodes and their values, set to None to disable
              "solverName":                   "implicitSolver",
              "checkForNanInf":               True,             # check if the solution vector contains nan or +/-inf values, if yes, an error is printed. This is a time-consuming check.
              "nAdditionalFieldVariables":    0,
              "additionalSlotNames":          [],
              
              "FiniteElementMethod" : {
                "meshName":               "MeshFiber",
                "prefactor":              Conductivity/(Am*Cm),
                "solverName":             "implicitSolver",
                "inputMeshIsGlobal":      True,
                "slotName":               "vm",
              },
              
              # output writer only for the diffusion variable (i.e. state "Vm")
              "OutputWriter" : [
                {"format": "PythonFile", "outputInterval": int(1./dt_1D*output_timestep), "filename": "out/vm", "binary": True, "onlyNodalValues": False, "fileNumbering": "incremental"},
                {"format": "Paraview",   "outputInterval": int(1./dt_1D*output_timestep), "filename": "out/vm", "binary": True, "fixedFormat": False, "combineFiles": True, "fileNumbering": "incremental"},
                #{"format": "ExFile", "filename": "out/fiber", "outputInterval": 1e5, "sphereSize": "0.02*0.02*0.02"},
              ],
            },
          },
        }
      }
    }
  }
}
