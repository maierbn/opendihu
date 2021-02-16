# Monodomain with Hodgkin-Huxley, one fiber with 100 elements
#
# There are two versions in this example: fast_fiber and not_fast_fiber.
# The fast version uses the FastMonodomainSolver which is a more efficient specialization of the Strang<Heun,ImplicitEuler> for Hodgkin-Huxley (also for Shorten).
# The not_fast version computes the same, but without FastMonodomainSolver.
#
# You can run: 
# ./fast_fiber ../settings_fast_fiber.py          # takes 1s
# ./not_fast_fiber ../settings_fast_fiber.py      # takes 8s

import numpy as np

# timing parameters
stimulation_frequency = 10.0      # [1/ms] frequency if which stimulation current can be switched on and off
dt_0D = 2e-4                      # timestep width of ODEs, cellml integration
dt_1D = 2e-3
dt_splitting = 2e-3
dt_3D = 0.5
output_timestep = 0.5
end_time = 100.0
n_elements = 100

stimulation_frequency = 100*1e-3   # [Hz]*1e-3 = [ms^-1]
frequency_jitter = 0.1
call_enable_begin = 15.0  # [s]*1e3 = [ms]

fiber_distribution_file = "../../../input/MU_fibre_distribution_10MUs.txt"
firing_times_file = "../../../input/MU_firing_times_always.txt"

# load MU distribution and firing times
fiber_distribution = np.genfromtxt(fiber_distribution_file, delimiter=" ")
firing_times = np.genfromtxt(firing_times_file)

def fiber_gets_stimulated(fiber_no, frequency, current_time):
  """
  determine if fiber fiber_no gets stimulated at simulation time current_time
  """

  # determine motor unit
  alpha = 1.0   # 0.8
  mu_no = 0
  
  # determine if fiber fires now
  index = int(np.round(current_time * frequency))
  n_firing_times = np.size(firing_times,0)
  
  #if firing_times[index % n_firing_times, mu_no] == 1:
  #print("fiber {} is mu {}, t = {}, row: {}, stimulated: {} {}".format(fiber_no, mu_no, current_time, (index % n_firing_times), firing_times[index % n_firing_times, mu_no], "true" if firing_times[index % n_firing_times, mu_no] == 1 else "false"))
  print("fiber {} is mu {}, t = {}, row: {}, stimulated: {} {}".format(fiber_no, mu_no, current_time, (index % n_firing_times), firing_times[index % n_firing_times, mu_no], "true" if firing_times[index % n_firing_times, mu_no] == 1 else "false"))
  
  return firing_times[index % n_firing_times, mu_no] == 1
  
# callback function that can set states, i.e. prescribed values for stimulation
def set_specific_states(n_nodes_global, time_step_no, current_time, states, fiber_no):

  #print("call set_specific_states at time {}".format(current_time)) 
  
  # determine if fiber gets stimulated at the current time
  is_fiber_gets_stimulated = fiber_gets_stimulated(fiber_no, stimulation_frequency, current_time)

  if is_fiber_gets_stimulated:  
    # determine nodes to stimulate (center node, left and right neighbour)
    #innervation_zone_width_n_nodes = innervation_zone_width*100  # 100 nodes per cm
    innervation_node_global = int(n_nodes_global / 2)  # + np.random.randint(-innervation_zone_width_n_nodes/2,innervation_zone_width_n_nodes/2+1)
    nodes_to_stimulate_global = [innervation_node_global]
    if innervation_node_global > 0:
      nodes_to_stimulate_global.insert(0, innervation_node_global-1)
    if innervation_node_global < n_nodes_global-1:
      nodes_to_stimulate_global.append(innervation_node_global+1)
    print("t: {}, stimulate fiber {} at nodes {}".format(current_time, fiber_no, nodes_to_stimulate_global))

    for node_no_global in nodes_to_stimulate_global:
      states[(node_no_global,0,0)] = 20.0   # key: ((x,y,z),nodal_dof_index,state_no)


# define the config dict
config = {
  "solverStructureDiagramFile":     "solver_structure.txt",     # output file of a diagram that shows data connection between solvers
  "logFormat":                      "csv",                      # "csv" or "json", format of the lines in the log file, csv gives smaller files
  "scenarioName":                   "fast_fiber",               # scenario name to find the run in the log file
  "mappingsBetweenMeshesLogFile":   "",                         # a log file about mappings between meshes, here we do not want that because there are no mappings
  "Meshes": {
    "MeshFiber_0": {
      "nElements": [n_elements],
      "physicalExtent": [n_elements/100.],
      "physicalOffset": [0],
      "inputMeshIsGlobal": True,
    }
  },
  "Solvers": {
    "implicitSolver": {     # solver for the implicit timestepping scheme of the diffusion time step
      "maxIterations":      1e4,
      "relativeTolerance":  1e-10,
      "absoluteTolerance":  0,
      "dumpFormat": "",
      "dumpFilename": "",
      "solverType": "gmres",
      "preconditionerType": "none"
    },
  },
  "RepeatedCall": {
    "timeStepWidth":          dt_splitting,
    "timeStepOutputInterval": 100,
    "endTime":                end_time,
    "MultipleInstances": {
      "ranksAllComputedInstances":  [0],
      "nInstances":                 1,
      "instances": 
      [{
        "ranks": [0],
        "StrangSplitting": {
          #"numberTimeSteps": 1,
          "timeStepWidth":          dt_splitting,
          "timeStepOutputInterval": 100,
          "endTime":                dt_splitting,
          "connectedSlotsTerm1To2": [0],   # transfer slot 0 = state Vm from Term1 (CellML) to Term2 (Diffusion)
          "connectedSlotsTerm2To1": [0],   # transfer the same data back

          "Term1": {      # CellML, i.e. reaction term of Monodomain equation
            "MultipleInstances": {
              "logKey":             "duration_subdomains_z",
              "nInstances":         1,
              "instances": 
              [{
                "ranks":                          [0],    # these rank nos are local nos to the outer instance of MultipleInstances, i.e. from 0 to number of ranks in z direction
                "Heun" : {
                  "timeStepWidth":                dt_0D,  # 5e-5
                  "logTimeStepWidthAsKey":        "dt_0D",
                  "durationLogKey":               "duration_0D",
                  "initialValues":                [],
                  "timeStepOutputInterval":       1e4,
                  "inputMeshIsGlobal":            True,
                  "checkForNanInf":               False,
                  "dirichletBoundaryConditions":  {},
                  "dirichletOutputFilename":      None,                # filename for a vtp file that contains the Dirichlet boundary condition nodes and their values, set to None to disable
                  "nAdditionalFieldVariables":    0,
                  "additionalSlotNames":          [],
                    
                  "CellML" : {
                    "modelFilename":                         "../../../input/hodgkin_huxley_1952.c",                          # input C++ source file or cellml XML file
                    "compilerFlags":                          "-fPIC -O3 -march=native -shared ",
                    "optimizationType":                       "vc",                       # "vc", "simd", "openmp" type of generated optimizated source file
                    "approximateExponentialFunction":         True,                      # if optimizationType is "vc", whether the exponential function exp(x) should be approximate by (1+x/n)^n with n=1024
                    "maximumNumberOfThreads":                 0,                          # if optimizationType is "openmp", the maximum number of threads to use. Default value 0 means no restriction.      
                    "initializeStatesToEquilibrium":          False,
                    
                    "setSpecificStatesFunction":              set_specific_states,                                             # callback function that sets states like Vm, activation can be implemented by using this method and directly setting Vm values, or by using setParameters/setSpecificParameters
                    "setSpecificStatesCallInterval":          0,                                                               # 0 means disabled
                    "setSpecificStatesCallFrequency":         stimulation_frequency,   # set_specific_states should be called variables.stimulation_frequency times per ms
                    "setSpecificStatesFrequencyJitter":       frequency_jitter, # random value to add or substract to setSpecificStatesCallFrequency every stimulation, this is to add random jitter to the frequency
                    "setSpecificStatesRepeatAfterFirstCall":  0.1,                                                            # [ms] simulation time span for which the setSpecificStates callback will be called after a call was triggered
                    "setSpecificStatesCallEnableBegin":       call_enable_begin,# [ms] first time when to call setSpecificStates
                    "additionalArgument":                     0,
                    "stimulationLogFilename":                 "out/stimulation_log.txt",
                    
                    # This uses the old specification of states, algebraics and parameters which needs numbers of the variables. 
                    # The new approach would be to use the option "mappings" which uses the names.
                    "statesForTransfer":                      0,      # which state values to use in further computation, Shorten / Hodgkin Huxley: state 0 = Vm
                    "algebraicsForTransfer":                  [],      # which algebraic values to use in further computation
                    "parametersForTransfer":                  [],
                    "parametersUsedAsAlgebraic":              [],     #[32],       # list of algebraic value indices, that will be set by parameters. Explicitely defined parameters that will be copied to algebraics, this vector contains the indices of the algebraic array. This is ignored if the input is generated from OpenCMISS generated c code.
                    "parametersUsedAsConstant":               [2],    #[65],           # list of constant value indices, that will be set by parameters. This is ignored if the input is generated from OpenCMISS generated c code.
                    "parametersInitialValues":                [0.0],  #[0.0, 1.0],      # initial values for the parameters: I_Stim, l_hs
                    "meshName":                               "MeshFiber_0",
                  },
                },
              }],
            }
          },
          "Term2": {     # Diffusion
            "MultipleInstances": {
              "nInstances": 1,
              "instances": 
              [{
                "ranks":                         [0],    # these rank nos are local nos to the outer instance of MultipleInstances, i.e. from 0 to number of ranks in z direction
                "ImplicitEuler" : {
                  "initialValues":               [],
                  #"numberTimeSteps":            1,
                  "timeStepWidth":               dt_1D,  # 1e-5
                  "timeStepWidthRelativeTolerance": 1e-10,
                  "logTimeStepWidthAsKey":       "dt_1D",
                  "durationLogKey":              "duration_1D",
                  "timeStepOutputInterval":      1e4,
                  "dirichletBoundaryConditions": {}, #{0: -75.0036, -1: -75.0036},
                  "dirichletOutputFilename":     None,                # filename for a vtp file that contains the Dirichlet boundary condition nodes and their values, set to None to disable
                  "inputMeshIsGlobal":           True,
                  "checkForNanInf":              False,
                  "solverName":                  "implicitSolver",
                  "nAdditionalFieldVariables":   0,
                  "additionalSlotNames":         [],
                  
                  "FiniteElementMethod" : {
                    "inputMeshIsGlobal":         True,
                    "meshName":                  "MeshFiber_0",
                    "prefactor":                 0.03,  # resolves to Conductivity / (Am * Cm)
                    "solverName":                "implicitSolver",
                    "slotName":                  "",
                  },
                  "OutputWriter" : [
                    #{"format": "Paraview", "outputInterval": int(1./variables.dt_1D*variables.output_timestep), "filename": "out/fiber_"+str(fiber_no), "binary": True, "fixedFormat": False, "combineFiles": True},
                    #{"format": "Paraview", "outputInterval": 1./variables.dt_1D*variables.output_timestep, "filename": "out/fiber_"+str(i)+"_txt", "binary": False, "fixedFormat": False},
                    #{"format": "ExFile", "filename": "out/fiber_"+str(i), "outputInterval": 1./variables.dt_1D*variables.output_timestep, "sphereSize": "0.02*0.02*0.02"},
                    #{"format": "PythonFile", "filename": "out/fiber_"+str(i), "outputInterval": 1./variables.dt_1D*variables.output_timestep, "binary":True, "onlyNodalValues":True},
                  ]
                },
              }],
              "OutputWriter" : [
                {"format": "PythonFile", "outputInterval": int(1./dt_splitting*output_timestep), "filename": "out/fast/fibers", "binary": True, "fixedFormat": False, "combineFiles": True, "onlyNodalValues": True, "fileNumbering": "incremental"}
              ]
            },
          },
        }
      }]
    },
    "fiberDistributionFile":    fiber_distribution_file,   # for FastMonodomainSolver, e.g. MU_fibre_distribution_3780.txt
    "firingTimesFile":          firing_times_file,         # for FastMonodomainSolver, e.g. MU_firing_times_real.txt
    "onlyComputeIfHasBeenStimulated": False,                          # only compute fibers after they have been stimulated for the first time
    "disableComputationWhenStatesAreCloseToEquilibrium": False,       # optimization where states that are close to their equilibrium will not be computed again
    "valueForStimulatedPoint":  20,                                   # to which value of Vm the stimulated node should be set
  }
}
