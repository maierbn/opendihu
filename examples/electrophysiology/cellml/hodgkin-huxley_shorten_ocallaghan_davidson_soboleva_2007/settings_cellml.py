# Single CellML problem, hodgkin-huxley_shorten_ocallaghan_davidson_soboleva_2007
# 
    
import numpy as np
import matplotlib 
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# parameters
cellml_file = "../../../input/2020_05_28_hodgkin-huxley_shorten_ocallaghan_davidson_soboleva_2007.cellml"
end_time = 30.0                   # [ms] end time of the simulation
dt = 1e-5                           # [ms] timestep width of ODEs

stimulation_frequency = 100*1e-3    # [ms^-1] sampling frequency of stimuli in firing_times_file, in stimulations per ms, number before 1e-3 factor is in Hertz.
output_timestep = 1                 # [ms] timestep for output files


# callback function that can set states, i.e. prescribed values for stimulation
def set_specific_states(n_nodes_global, time_step_no, current_time, states, additional_argument):
  # n_dofs_global:  (int) global number of dofs in the mesh, i.e. number of CellML instances to be computed
  # timestep_no:    (int)   current time step number, advances by the value of "setSpecificParametersCallInterval"
  # current_time:   (float) the current simulation time
  # global_states:  (dict)  initially an empty dict, the states to be changed should be indicated in this dict (see below)
  # additional_argument: The value of the option "additionalArgument", can be any Python object.
  
  print("set_specific_states at t: {}, {}".format(time_step_no, current_time))
    
  if 9 < current_time < 11:
    
      key = ((0),0,0)        # key: ((x,y,z),nodal_dof_index,state_no)
      states[key] = 20

# callback function that can set parameters, i.e. stimulation current. This is used to stimulate the CellML problem
def set_specific_parameters(n_nodes_global, time_step_no, current_time, parameters, additional_argument):
  
  print("set_specific_parameters at t: {}, {}".format(time_step_no, current_time))
  
  # do not stimulate for t > 1
  if current_time > 0.1:
    return
  
  # determine nodes to stimulate (center node)  (this is not relevant here, as we only have 1 node which has no. 0)
  # but this code is included here if sb. wants to copy it to a different scenario where we have multiple nodes
  center_node = int(n_nodes_global / 2)
  nodes_to_stimulate_global = [center_node]
  
  # add 10 neighbours to the left and right (not really here, as we have only 1 node)
  for k in range(10):
    if center_node-k >= 0:
      nodes_to_stimulate_global.insert(0, center_node-k)
    if center_node+k <= n_nodes_global-1:
      nodes_to_stimulate_global.append(center_node+k)
  
  # set stimulation value of 40, set parameter no. 0 of node 0. (there is only one node here)
  stimulation_current = 40.
  
  for node_no_global in nodes_to_stimulate_global:
    parameters[(node_no_global, 0, 0)] = stimulation_current   # key: ([x,y,z], nodalDofIndex, parameterNo)


xdata = []
active_stress_data = []
vm_data = []

# callback function that receives the whole result values and produces plots while the simulation is running
def handle_result(n_instances, time_step_no, current_time, states, algebraics, name_information, additional_argument):
    
  # asign some states to variables
  Vm = states[name_information["stateNames"].index("membrane/V")]
  Ca_1 = states[name_information["stateNames"].index("razumova/Ca_1")]
  A_1 = states[name_information["stateNames"].index("razumova/A_1")]
  A_2 = states[name_information["stateNames"].index("razumova/A_2")]
  x_1 = states[name_information["stateNames"].index("razumova/x_1")]
  x_2 = states[name_information["stateNames"].index("razumova/x_2")]
  
  # assign some algebraics to variables
  active_stress = algebraics[name_information["algebraicNames"].index("razumova/activestress")]
  activation = algebraics[name_information["algebraicNames"].index("razumova/activation")]
    
  print("handle_result, t: {t}, Vm: {Vm}, Ca_1: {Ca_1}, x_1: {x_1}, x_2: {x_2}, A_1:{A_1}, A_2:{A_2}, active stress: {active_stress}, activation: {activation}".format(t=current_time, Vm=Vm, Ca_1=Ca_1, A_1=A_1, A_2=A_2, x_1=x_1, x_2=x_2, active_stress=active_stress, activation=activation))
  
  # save values for plot
  xdata.append(current_time)
  active_stress_data.append(algebraics[9])
  vm_data.append(states[0])
  
  # at the end in the last call to this function, create a plot
  if current_time >= end_time - 1 - dt:
    
    # plot values of Vm and active stress
    plt.figure(1)
    plt.clf()
    ax1 = plt.gca()
    ax1.plot(xdata, vm_data, 'go-', label='$V_m$')
    plt.xlabel('t')
    plt.ylabel('$V_m$')
    ax2 = ax1.twinx()
    ax2.plot(xdata, active_stress_data, 'ro-', label='active stress, $\gamma$')
    plt.ylabel('$\gamma$')    
    
    # ask matplotlib for the plotted objects and their labels
    lines, labels = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax2.legend(lines + lines2, labels + labels2, loc=0)
    plt.tight_layout()
    
    filename = "active_stress_plot.png"
    print("created file \"{}\"".format(filename))
    plt.savefig(filename)
    #plt.draw()
  
xdata = []
gamma_data = []
vm_data = []

config = {
  "scenarioName":                   "cellml",   
  "logFormat":                      "csv",
  "solverStructureDiagramFile":     "solver_structure.txt",     # output file of a diagram that shows data connection between solvers
  "mappingsBetweenMeshesLogFile":   "mappings_between_meshes.txt",   # log file for mappings between meshes
  "ExplicitEuler" : {
    "timeStepWidth":          dt,     # dt of solver
    "endTime" :               end_time, # end simulation time of solver
    "initialValues":          [],     # initial values (not used)
    "timeStepOutputInterval": 1e5,    # the interval when the current time will be printed in the console (e.g. 'Explicit Euler, timestep 100000/10000000, t=1')
    "inputMeshIsGlobal":      True,   # for the mesh, not relevant here as we have no elements, only one node
    "checkForNanInf":         True,   # check if the solution vector contains nan or +/-inf values, if yes, an error is printed. This is a time-consuming check.
    "nAdditionalFieldVariables": 0,   # only revelant if there are multiple nested solvers and they transfer additional data
    "additionalSlotNames": [],        # the slot names of the additional field variables
    "dirichletBoundaryConditions": {},    # we do not set dirichlet BC
    "dirichletOutputFilename":     None,  # filename for a vtp file that contains the Dirichlet boundary condition nodes and their values, set to None to disable
    
    "OutputWriter" : [
       #{"format": "Callback", "outputInterval": 1e4, "callback": callback},
      #{"format": "Paraview", "filename": "out", "binaryOutput": "false", "fixedFormat": False, "outputInterval": 1},
      {"format": "PythonFile", "filename": "out/result", "outputInterval":  int(1./dt*output_timestep), "binary": True, "onlyNodalValues": True, "fileNumbering": "incremental"},
    ],

    "CellML" : {
      "nElements": 0,                 # information on the mesh to use, here we have no element, this means only 1 dof, i.e. 1 instance of the CellML problem
      "inputMeshIsGlobal": True,      # information on the mesh to use
      
      "modelFilename":                          cellml_file,                # input C++ source file or cellml XML file        
      #"statesInitialValues": [...],                                        # use this to set initial values for all states, if not set, the default values are parsed from the modelFilename (which is usually a good idea)
      "initializeStatesToEquilibrium":          False,                      # if the equilibrium values of the states should be computed before the simulation starts
      "initializeStatesToEquilibriumTimestepWidth": 1e-4,                   # if initializeStatesToEquilibrium is enable, the timestep width to use to solve the equilibrium equation
      
      # optimization parameters
      "optimizationType":                       "vc",                       # "vc", "simd", "openmp" type of generated optimizated source file
      "approximateExponentialFunction":         True,                       # if optimizationType is "vc", whether the exponential function exp(x) should be approximate by (1+x/n)^n with n=1024
      "compilerFlags":                          "-fPIC -O3 -march=native -shared ",  # compiler flags used to compile the optimized model code
      "maximumNumberOfThreads":                 0,                          # if optimizationType is "openmp", the maximum number of threads to use. Default value 0 means no restriction.
      
      # stimulation, directly setting Vm value
      "setSpecificStatesFunction":              set_specific_states,        # callback function that sets states like Vm, activation can be implemented by using this method and directly setting Vm values, or by using setParameters/setSpecificParameters
      "setSpecificStatesCallInterval":          int(1/stimulation_frequency/dt), # 0 means disabled, use CallFrequency instead
      "setSpecificStatesCallFrequency":         None, #stimulation_frequency,      # set_specific_states should be called stimulation_frequency times per ms
      "setSpecificStatesCallEnableBegin":       0,                          # [ms] first time when to call setSpecificStates
      "setSpecificStatesRepeatAfterFirstCall":  None,                       # [ms] simulation time span for which the setSpecificStates callback will be called after a call was triggered
      "setSpecificStatesFrequencyJitter":       0,                          # random value to add or substract to setSpecificStatesCallFrequency every stimulation, this is to add random jitter to the frequency
                           
      # stimulation, setting stimulation current
      "setSpecificParametersFunction":          None, #set_specific_parameters,    # callback function that sets parameters like stimulation current
      "setSpecificParametersCallInterval":      int(1/stimulation_frequency/dt), # set_specific_parameters should be called every x ms
      "additionalArgument":                     None,                       # additional last argument for set_specific_parameters
                      
      "handleResultFunction":                   handle_result,              # callback function that gets all current values and can do something with them
      "handleResultCallInterval":               1e5,                        # interval in which handle_result will be called
      "handleResultFunctionAdditionalParameter": None,                      # additional last argument for handle_result
      
      "stimulationLogFilename":                 "out/stimulation.log",      # filename of a log file that contains all stimulation events
      
      # parameters: I_Stim, l_hs. I_Stim (the stimulation current) is used as "algebraic" no. 32, l_hs (the relative half-sarcomere-length) is used as "constant" 65.
      "parametersInitialValues":                [0.0, 1.0],              # initial values for all parameters: I_Stim, l_hs
      
      # If you only have the C source file of the CellML model, it is easier specify these values as indices of the ALGEBRAIC and CONSTANTS variables, as follows:
      #"parametersUsedAsAlgebraic":              [32],                       # use ALGEBRAIC[32] as parameter, list of algebraic value indices that will be set by parameters. Explicitely defined parameters that will be copied to algebraics, this vector contains the indices of the algebraic array.
      #"parametersUsedAsConstant":               [65],                       # use CONSTANTS[65] as parameter, list of constant value indices that will be set by parameters.
      
      # If you inspect the CellML in OpenCOR, it is easier to specify the names of the variables as "environment/name", as follows:
      "mappings": {
        ("parameter", 0):           "membrane/i_Stim",                      # parameter 0 is I_stim
        ("parameter", 1):           "razumova/L_x",                         # parameter 1 is fiber stretch λ
      },
      
      "statesForTransfer":                      [],                         # in case of coupled solvers the states to transfer to the other solver, here we have no coupled solvers therefore this is not relevant, state 0 = Vm, rate 28 = gamma
      "algebraicsForTransfer":                  [],                         # in case of coupled solvers the algebraics to transfer to the other solver, here we have no coupled solvers therefore this is not relevant
      "parametersForTransfer":                  [],                         # in case of coupled solvers the parameters to transfer to the other solver, here we have no coupled solvers therefore this is not relevant
    },
  },
}
