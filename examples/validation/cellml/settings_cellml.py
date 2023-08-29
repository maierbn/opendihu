# Single CellML problem, Shorten or Hodgkin-Huxley
# No time dependent stimulation, a constant stimulation current is applied.
    

import opendihu
if "shorten" in opendihu.program_name:
  model_filename = "../../../electrophysiology/input/shorten_ocallaghan_davidson_soboleva_2007_no_stim.cellml"
  stimulation_current = 50
  mappings = {
    ("parameter", 0):           "wal_environment/I_HH",                 # parameter 0 is I_stim
  }
  output_directory = "out_shorten"

elif "hodgkin_huxley" in opendihu.program_name:
  model_filename = "../../../electrophysiology/input/hodgkin_huxley_1952.cellml"
  stimulation_current = 10
  mappings = {
    ("parameter", 0):           "membrane/i_Stim",                 # parameter 0 is I_stim
  }
  output_directory = "out_hodgkin_huxley"

config = {
  "scenarioName":                   "cellml",   
  "logFormat":                      "csv", # "csv" or "json", format of the lines in the log file, csv gives smaller files
  "solverStructureDiagramFile":     "solver_structure.txt",     # output file of a diagram that shows data connection between solvers
  "mappingsBetweenMeshesLogFile":   "mappings_between_meshes.txt",   # log file for mappings between meshes
  "Heun" : {
    "timeStepWidth":          1e-5,   # dt of solver
    "endTime" :               35.0,   # end simulation time of solver
    "initialValues":          [],     # initial values (not used)
    "timeStepOutputInterval": 1e4,    # the interval when the current time will be printed in the console (e.g. 'Explicit Euler, timestep 100000/10000000, t=1')
    "inputMeshIsGlobal":      True,   # for the mesh, not relevant here as we have no elements, only one node
    "checkForNanInf":         False,   # check if the solution vector contains nan or +/-inf values, if yes, an error is printed. This is a time-consuming check.
    "nAdditionalFieldVariables": 0,   # only revelant if there are multiple nested solvers and they transfer additional data
    "additionalSlotNames": [],        # the slot names of the additional field variables
    "dirichletBoundaryConditions": {},    # we do not set dirichlet BC
    "dirichletOutputFilename":     None,  # filename for a vtp file that contains the Dirichlet boundary condition nodes and their values, set to None to disable
    
    "OutputWriter" : [
       #{"format": "Callback", "outputInterval": 1e4, "callback": callback},
      #{"format": "Paraview", "filename": "out", "binaryOutput": "false", "fixedFormat": False, "outputInterval": 1},
      {"format": "PythonFile", "filename": f"{output_directory}/result", "outputInterval": 1e4, "binary": True, "onlyNodalValues": True, "fileNumbering": "incremental"},
    ],

    "CellML" : {
      "nElements": 0,                 # information on the mesh to use, here we have no element, this means only 1 dof, i.e. 1 instance of the CellML problem
      "inputMeshIsGlobal": True,      # information on the mesh to use
      
      "modelFilename":                          model_filename,             # input C++ source file or cellml XML file        
      #"statesInitialValues": [...],                                        # use this to set initial values for all states, if not set, the default values are parsed from the modelFilename (which is usually a good idea)
      "initializeStatesToEquilibrium":          False,                      # if the equilibrium values of the states should be computed before the simulation starts
      "initializeStatesToEquilibriumTimestepWidth": 1e-4,                   # if initializeStatesToEquilibrium is enable, the timestep width to use to solve the equilibrium equation
      
      # optimization parameters
      "optimizationType":                       "vc",                       # "vc", "simd", "openmp" type of generated optimizated source file
      "approximateExponentialFunction":         False,                       # if optimizationType is "vc", whether the exponential function exp(x) should be approximate by (1+x/n)^n with n=1024
      "compilerFlags":                          "-fPIC -O3 -march=native -shared ",  # compiler flags used to compile the optimized model code
      "maximumNumberOfThreads":                 0,                          # if optimizationType is "openmp", the maximum number of threads to use. Default value 0 means no restriction.
       
      "setSpecificStatesCallFrequency":         0,                          # set_specific_states should be called stimulation_frequency times per ms
      "setSpecificStatesCallEnableBegin":       0,                          # [ms] first time when to call setSpecificStates
      "setSpecificStatesRepeatAfterFirstCall":  0,                          # [ms] simulation time span for which the setSpecificStates callback will be called after a call was triggered
                           
      "setSpecificParametersFunction":          None,                       # callback function that sets parameters like stimulation current
      "setSpecificParametersCallInterval":      1e3,                        # set_specific_parameters should be called every x ms
      "setSpecificStatesFrequencyJitter":       0,                          # random value to add or substract to setSpecificStatesCallFrequency every stimulation, this is to add random jitter to the frequency
      "additionalArgument":                     None,                       # additional last argument for set_specific_parameters
                      
      "handleResultFunction":                   None,                       # callback function that gets all current values and can do something with them, in this case, handle_result will output a plot to out.png
      "handleResultCallInterval": 1e4,                                      # interval in which handle_result will be called
      "handleResultFunctionAdditionalParameter": None,                      # additional last argument for handle_result
      
      "stimulationLogFilename":                 "out/stimulation.log",      # filename of a log file that contains all stimulation events
      
      # parameters: I_Stim (the stimulation current) is used as "algebraic" no. 32
      "parametersInitialValues":                [stimulation_current],                     # initial values for all parameters: I_Stim
      
      # If you only have the C source file of the CellML model, it is easier specify these values as indices of the ALGEBRAIC and CONSTANTS variables, as follows:
      #"parametersUsedAsAlgebraic":              [32],                       # use ALGEBRAIC[32] as parameter, list of algebraic value indices that will be set by parameters. Explicitely defined parameters that will be copied to algebraics, this vector contains the indices of the algebraic array.
      #"parametersUsedAsConstant":               [65],                       # use CONSTANTS[65] as parameter, list of constant value indices that will be set by parameters.
      
      # If you inspect the CellML in OpenCOR, it is easier to specify the names of the variables as "environment/name", as follows:
      "mappings":                               mappings,
      
      "statesForTransfer":                      [],                         # in case of coupled solvers the states to transfer to the other solver, here we have no coupled solvers therefore this is not relevant, state 0 = Vm, rate 28 = gamma
      "algebraicsForTransfer":                  [],                         # in case of coupled solvers the algebraics to transfer to the other solver, here we have no coupled solvers therefore this is not relevant
      "parametersForTransfer":                  [],                         # in case of coupled solvers the parameters to transfer to the other solver, here we have no coupled solvers therefore this is not relevant
    },
  },
}
