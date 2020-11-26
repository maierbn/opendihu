FastMonodomainSolver
======================

This is a very efficient implementation of multiple fibers where the monodomain equation is solved on each.
The fibers all have the same number of elements. The fibers are more or less parallel. 
The parallel partitioning is arbitrary, i.e. fibers can be subdivided to different processes and different fibers do not need to share processes.

All nodes of this geometry form a cuboid in index space :math:`(i,j,k)`, where the `k` index runs over the nodes of a fiber and :math:`(i,j)` specify the fiber in a 2D grid of fibers.
The partitioning is obtained by dividing this cuboid by axis-aligned plane cuts in all three dimensions.

The *FastMonodomainSolver* reuses nested solvers, as given below.

Usage
----------
The following shows the typical usage in the C++ source file. Only the number of states and intermediates, 4 and 9 in this example and the timestepping scheme of the diffusion, `TimeSteppingScheme::ImplicitEuler`, in this example, can be changed.

.. code-block:: c
  :linenos:

  FastMonodomainSolver<                        // a wrapper that improves performance of multidomain
    Control::MultipleInstances<                       // fibers
      OperatorSplitting::Strang<
        Control::MultipleInstances<
          TimeSteppingScheme::Heun<                   // fiber reaction term
            CellmlAdapter<
              4, 9,  // nStates,nIntermediates: 57,1 = Shorten, 4,9 = Hodgkin Huxley
              FunctionSpace::FunctionSpace<
                Mesh::StructuredDeformableOfDimension<1>,
                BasisFunction::LagrangeOfOrder<1>
              >
            >
          >
        >,
        Control::MultipleInstances<
          TimeSteppingScheme::ImplicitEuler<          // fiber diffusion
            SpatialDiscretization::FiniteElementMethod<
              Mesh::StructuredDeformableOfDimension<1>,
              BasisFunction::LagrangeOfOrder<1>,
              Quadrature::Gauss<2>,
              Equation::Dynamic::IsotropicDiffusion
            >
          >
        >
      >
    >
  >

The two template arguments of `CellmlAdapter`, the *number of states* and *number of intermediates* can be adjusted to fit the subcellular CellML model.
Instead of ``TimeSteppingScheme::ImplicitEuler`` for the diffusion problem, ``TimeSteppingScheme::CrankNicholson`` can be used. All other templates must appear exactly as given above.

The *FastMonodomainSolver* solves the same equations as the nested solver would (just as if lines 1 and 27 were not present). The discretization is also the same. A difference is that the diffusion problem is solved in serial using Thomas' algorithm, i.e. in linear time. For this purpose, the data of a single fiber is communicated to a single rank where it gets solved. At he end of the timestep, the results are communicated back.

The improved performance is by roughly a factor of 10. The reason is that the 1D diffusion problem which is a tri-diagonal system gets solved serially and by a Thomas' algorithm which has linear time complexity. All values of a fiber are communicated to a single rank at the beginning of the time span. (Different ranks for different fibers). The fiber is then solved completely on this one rank for all specified timesteps. 
This involves the Strang splitting consisting of solving the subcellular model and the diffusion problem.
At the end, the values are communicated back to the original process. Consequently, the *FastMonodomainSolver* appears to surrounding solvers like its nested solvers with a cubes-like partitioning, but internally the fibers are not split across processors.

For the subcellular model, efficent code of the whole Heun scheme, using `"vc"`, is generated and executed. This is also faster than if `Heun` and `CellMLAdapter` are nested.

Python options
---------------------

The normal python options for the nested solvers are needed. The FastMonodomainSolver extracts most of its settings from the settings of the nested solvers. There are some additional options at the end that are specific to the FastMonodomainSolver.

The following example of python options is taken from the `fibers/fibers_fat_emg` example. Note the additional options at the end. Also note that there is no new top-level key for the FastMonodomainSolver, it simply uses the existing keys.

.. code-block:: python

  {
    "MultipleInstances": {
      "logKey":                     "duration_subdomains_xy",
      "ranksAllComputedInstances":  list(range(n_ranks)),
      "nInstances":                 variables.n_subdomains_xy,
      "instances": 
      [{
        "ranks": list(range(subdomain_coordinate_y*variables.n_subdomains_x + subdomain_coordinate_x, n_ranks, variables.n_subdomains_x*variables.n_subdomains_y)),
        "StrangSplitting": {
          #"numberTimeSteps": 1,
          "timeStepWidth":          variables.dt_splitting,  # 1e-1
          "logTimeStepWidthAsKey":  "dt_splitting",
          "durationLogKey":         "duration_monodomain",
          "timeStepOutputInterval": 100,
          "endTime":                variables.dt_splitting,
          "connectedSlotsTerm1To2": [0],   # transfer slot 0 = state Vm from Term1 (CellML) to Term2 (Diffusion)
          "connectedSlotsTerm2To1": [0],   # transfer the same back, this avoids data copy

          "Term1": {      # CellML, i.e. reaction term of Monodomain equation
            "MultipleInstances": {
              "logKey":             "duration_subdomains_z",
              "nInstances":         n_fibers_in_subdomain_x(subdomain_coordinate_x)*n_fibers_in_subdomain_y(subdomain_coordinate_y),
              "instances": 
              [{
                "ranks":                          list(range(variables.n_subdomains_z)),    # these rank nos are local nos to the outer instance of MultipleInstances, i.e. from 0 to number of ranks in z direction
                "Heun" : {
                  "timeStepWidth":                variables.dt_0D,                         # timestep width of 0D problem
                  "logTimeStepWidthAsKey":        "dt_0D",                                 # key under which the time step width will be written to the log file
                  "durationLogKey":               "duration_0D",                           # log key of duration for this solver
                  "timeStepOutputInterval":       1e4,                                     # how often to print the current timestep
                  "initialValues":                [],                                      # no initial values are specified
                  "dirichletBoundaryConditions":  {},                                      # no Dirichlet boundary conditions are specified
                  "dirichletOutputFilename":      None,                                    # filename for a vtp file that contains the Dirichlet boundary condition nodes and their values, set to None to disable
                  
                  "inputMeshIsGlobal":            True,                                    # the boundary conditions and initial values would be given as global numbers
                  "checkForNanInf":               True,                                    # abort execution if the solution contains nan or inf values
                  "nAdditionalFieldVariables":    0,                                       # number of additional field variables
                  "additionalSlotNames":          [],                                      # names for the additional slots
                    
                  "CellML" : {
                    "modelFilename":                          variables.cellml_file,                          # input C++ source file or cellml XML file
                    #"statesInitialValues":                   [],                                             # if given, the initial values for the the states of one instance
                    "initializeStatesToEquilibrium":          False,                                          # if the equilibrium values of the states should be computed before the simulation starts
                    "initializeStatesToEquilibriumTimestepWidth": 1e-4,                                       # if initializeStatesToEquilibrium is enable, the timestep width to use to solve the equilibrium equation
                    
                    # optimization parameters
                    "optimizationType":                       "vc" if variables.use_vc else "simd",           # "vc", "simd", "openmp" type of generated optimizated source file
                    "approximateExponentialFunction":         True,                                           # if optimizationType is "vc", whether the exponential function exp(x) should be approximate by (1+x/n)^n with n=1024
                    "compilerFlags":                          "-fPIC -O3 -march=native -shared ",             # compiler flags used to compile the optimized model code
                    "maximumNumberOfThreads":                 0,                                              # if optimizationType is "openmp", the maximum number of threads to use. Default value 0 means no restriction.
                    
                    # stimulation callbacks
                    #"libraryFilename":                       "cellml_simd_lib.so",                           # compiled library
                    #"setSpecificParametersFunction":         set_specific_parameters,                        # callback function that sets parameters like stimulation current
                    #"setSpecificParametersCallInterval":     int(1./variables.stimulation_frequency/variables.dt_0D),         # set_specific_parameters should be called every 0.1, 5e-5 * 1e3 = 5e-2 = 0.05
                    "setSpecificStatesFunction":              set_specific_states,                                             # callback function that sets states like Vm, activation can be implemented by using this method and directly setting Vm values, or by using setParameters/setSpecificParameters
                    #"setSpecificStatesCallInterval":         2*int(1./variables.stimulation_frequency/variables.dt_0D),       # set_specific_states should be called variables.stimulation_frequency times per ms, the factor 2 is needed because every Heun step includes two calls to rhs
                    "setSpecificStatesCallInterval":          0,                                                               # 0 means disabled
                    "setSpecificStatesCallFrequency":         variables.get_specific_states_call_frequency(fiber_no, motor_unit_no),   # set_specific_states should be called variables.stimulation_frequency times per ms
                    "setSpecificStatesFrequencyJitter":       variables.get_specific_states_frequency_jitter(fiber_no, motor_unit_no), # random value to add or substract to setSpecificStatesCallFrequency every stimulation, this is to add random jitter to the frequency
                    "setSpecificStatesRepeatAfterFirstCall":  0.01,                                                            # [ms] simulation time span for which the setSpecificStates callback will be called after a call was triggered
                    "setSpecificStatesCallEnableBegin":       variables.get_specific_states_call_enable_begin(fiber_no, motor_unit_no),# [ms] first time when to call setSpecificStates
                    "additionalArgument":                     fiber_no,                                       # last argument that will be passed to the callback functions set_specific_states, set_specific_parameters, etc.
                    
                    # parameters to the cellml model
                    "mappings":                               variables.mappings,                             # mappings between parameters and algebraics/constants and between outputConnectorSlots and states, algebraics or parameters, they are defined in helper.py
                    "parametersInitialValues":                variables.parameters_initial_values,            #[0.0, 1.0],      # initial values for the parameters: I_Stim, l_hs
                    
                    "meshName":                               "MeshFiber_{}".format(fiber_no),                # reference to the fiber mesh
                    "stimulationLogFilename":                 "out/stimulation.log",                          # a file that will contain the times of stimulations
                  },      
                  "OutputWriter" : [
                    {"format": "Paraview", "outputInterval": 1, "filename": "out/" + variables.scenario_name + "/0D_states({},{})".format(fiber_in_subdomain_coordinate_x,fiber_in_subdomain_coordinate_y), "binary": True, "fixedFormat": False, "combineFiles": True}
                  ] if variables.states_output else []
                  
                },
              } for fiber_in_subdomain_coordinate_y in range(n_fibers_in_subdomain_y(subdomain_coordinate_y)) \
                  for fiber_in_subdomain_coordinate_x in range(n_fibers_in_subdomain_x(subdomain_coordinate_x)) \
                    for fiber_no in [get_fiber_no(subdomain_coordinate_x, subdomain_coordinate_y, fiber_in_subdomain_coordinate_x, fiber_in_subdomain_coordinate_y)] \
                      for motor_unit_no in [get_motor_unit_no(fiber_no)]],
            }
          },
          "Term2": {     # Diffusion
            "MultipleInstances": {
              "nInstances": n_fibers_in_subdomain_x(subdomain_coordinate_x)*n_fibers_in_subdomain_y(subdomain_coordinate_y),
              "instances": 
              [{
                "ranks":                         list(range(variables.n_subdomains_z)),   # these rank nos are local nos to the outer instance of MultipleInstances, i.e. from 0 to number of ranks in z direction
                "ImplicitEuler" : {
                  "initialValues":               [],                                      # no initial values are given
                  #"numberTimeSteps":            1,
                  "timeStepWidth":               variables.dt_1D,                         # timestep width for the diffusion problem
                  "timeStepWidthRelativeTolerance": 1e-10,
                  "logTimeStepWidthAsKey":       "dt_1D",                                 # key under which the time step width will be written to the log file
                  "durationLogKey":              "duration_1D",                           # log key of duration for this solver
                  "timeStepOutputInterval":      1e4,                                     # how often to print the current timestep
                  "dirichletBoundaryConditions": {},                                      # old Dirichlet BC that are not used in FastMonodomainSolver: {0: -75.0036, -1: -75.0036},
                  "dirichletOutputFilename":     None,                                    # filename for a vtp file that contains the Dirichlet boundary condition nodes and their values, set to None to disable
                  "inputMeshIsGlobal":           True,                                    # initial values would be given as global numbers
                  "solverName":                  "diffusionTermSolver",                   # reference to the linear solver
                  "nAdditionalFieldVariables":   0,
                  "additionalSlotNames":         [],
                  "checkForNanInf":              False,
                  
                  "FiniteElementMethod" : {
                    "inputMeshIsGlobal":         True,
                    "meshName":                  "MeshFiber_{}".format(fiber_no),
                    "solverName":                "diffusionTermSolver",
                    "prefactor":                 get_diffusion_prefactor(fiber_no, motor_unit_no),  # resolves to Conductivity / (Am * Cm)
                    "slotName":                  "vm",
                  },
                  "OutputWriter" : [
                    #{"format": "Paraview", "outputInterval": int(1./variables.dt_1D*variables.output_timestep), "filename": "out/fiber_"+str(fiber_no), "binary": True, "fixedFormat": False, "combineFiles": True},
                    #{"format": "Paraview", "outputInterval": 1./variables.dt_1D*variables.output_timestep, "filename": "out/fiber_"+str(i)+"_txt", "binary": False, "fixedFormat": False},
                    #{"format": "ExFile", "filename": "out/fiber_"+str(i), "outputInterval": 1./variables.dt_1D*variables.output_timestep, "sphereSize": "0.02*0.02*0.02"},
                    #{"format": "PythonFile", "filename": "out/fiber_"+str(i), "outputInterval": 1./variables.dt_1D*variables.output_timestep, "binary":True, "onlyNodalValues":True},
                  ]
                },
              } for fiber_in_subdomain_coordinate_y in range(n_fibers_in_subdomain_y(subdomain_coordinate_y)) \
                  for fiber_in_subdomain_coordinate_x in range(n_fibers_in_subdomain_x(subdomain_coordinate_x)) \
                    for fiber_no in [get_fiber_no(subdomain_coordinate_x, subdomain_coordinate_y, fiber_in_subdomain_coordinate_x, fiber_in_subdomain_coordinate_y)] \
                      for motor_unit_no in [get_motor_unit_no(fiber_no)]],
              "OutputWriter" : variables.output_writer_fibers,
            },
          },
        }
      } if (subdomain_coordinate_x,subdomain_coordinate_y) == (variables.own_subdomain_coordinate_x,variables.own_subdomain_coordinate_y) else None
      for subdomain_coordinate_y in range(variables.n_subdomains_y)
          for subdomain_coordinate_x in range(variables.n_subdomains_x)]
    },
    "fiberDistributionFile":    variables.fiber_distribution_file,   # for FastMonodomainSolver, e.g. MU_fibre_distribution_3780.txt
    "firingTimesFile":          variables.firing_times_file,         # for FastMonodomainSolver, e.g. MU_firing_times_real.txt
    "onlyComputeIfHasBeenStimulated": variables.fast_monodomain_solver_optimizations,                          # only compute fibers after they have been stimulated for the first time
    "disableComputationWhenStatesAreCloseToEquilibrium": variables.fast_monodomain_solver_optimizations,       # optimization where states that are close to their equilibrium will not be computed again      
    "valueForStimulatedPoint":  variables.vm_value_stimulated,       # to which value of Vm the stimulated node should be set      
    "neuromuscularJunctionRelativeSize": 0.1,                          # range where the neuromuscular junction is located around the center, relative to fiber length. The actual position is draws randomly from the interval [0.5-s/2, 0.5+s/2) with s being this option. 0 means sharply at the center, 0.1 means located approximately at the center, but it can vary 10% in total between all fibers.
  }
  
Instead of the callback function `setSpecificStates` that would normally handle the stimulation, the FastMonodomainSolver does the stimulation differently. Calling the callback functions would be too slow. The same behaviour as with the standard `setSpecificStates` is implemented, respecting the options ``setSpecificStatesCallFrequency``, ``setSpecificStatesFrequencyJitter``, ``setSpecificStatesRepeatAfterFirstCall`` and ``setSpecificStatesCallEnableBegin``. Stimulation is done by setting Vm at a node to the value ``valueForStimulatedPoint``. Which node is determined by ``neuromuscularJunctionRelativeSize``.

.. _stimulation_times_1:
.. figure:: /settings/images/stimulation_times.svg
  :width: 80%
  
  Options that influence the stimulation. A time line is shown from left to right. The red blocks are time spans when there can be a stimulation. If there will be an actual stimulation, i.e. :math:`V_m` will be set to ``valueForStimulatedPoint``, depends on the entry in ``firingTimesFile`` for the current time step. Therefore it makes sense to use the file `"MU_firing_times_always.txt"` with these options as then, the spike trains are completely determined by the options `setSpecificStatesCallEnableBegin`, `setSpecificStatesCallFrequency` and `setSpecificStatesFrequencyJitter`..

fiberDistributionFile
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
This file contains the assignment of fibers to motor units. Such files are located in the `examples/electrophysiology/input` directory and have names ``MU_fibre_distribution*.txt``. 

The file contains a single line with space separated numbers, e.g. ``2 5 4 5 2 2 2 2 2 4 2 5 2 5 2 2 5 5 14 16 2 8 10 7 8``. Each number specifies the motor unit number of the next fiber, i.e., in this example fiber 0 is of MU 2, fiber 1 is of MU 5 etc. If there are more timestep or more motor units than entries in the file, the values wrap around, i.e. after the last column the first is used again. This file format is compatible with the OpenCMISS Iron examples.

firingTimesFile
^^^^^^^^^^^^^^^^^^^^^^^
This file specifies when which motor unit fires. Such files are located in the `examples/electrophysiology/input` directory and have names ``MU_firing_times*.txt``. 
Examples of available files are:

* ``MU_firing_times_always.txt``    Every motor unit fires in every timestep. This file is needed if the `setSpecificStatesCallFrequency` and `setSpecificStatesFrequencyJitter` functionality should be used.
* ``MU_firing_times_immediately.txt``  All motor units fire the first three timesteps at the beginning, then with a certain frequency which is exponentially distributed among the MUs. This can be used for debugging. In small simulation time spans also all MUs fire.
* ``MU_firing_times_once.txt``      All motor units fire only in the first timestep.
* ``MU_firing_times_real.txt``      The motor units fire in a frequency that is exponentially distributed over MU numbers.

It contains multiple lines, one for each time step. Every line consists of indications whether a motor unit fires (1) or not (0). The 0s and 1s are separated by spaces. This means the rows specify timestep numbers and the columns specify motor unit numbers. If there are more timestep or more motor units than entries in the file, the values wrap around, i.e. after the last column the first is used again. This file format is again compatible with the OpenCMISS Iron examples.

onlyComputeIfHasBeenStimulated
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
This option disabled computation of the Monodomain equation as long as the fiber has not been stimulated in therefore is in equilibrium.
This speeds up computation for a ramp scenario where a lot of MUs are inactive at the beginning. Set this option to False for runtime tests of the 0D and 1D problem.

disableComputationWhenStatesAreCloseToEquilibrium
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Similar to `onlyComputeIfHasBeenStimulated`, this checks whether the values have reached the equilibrium and then disables the computation.

valueForStimulatedPoint
^^^^^^^^^^^^^^^^^^^^^^^^^^^
This is the value that will be set for the transmembrane potential :math:`V_m` when it is stimulated.

neuromuscularJunctionRelativeSize
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Relative range of the position of the neuromuscular junction. The neuromuscular junction is the point on a fiber where the nerve innervates and the fiber gets stimulated. This value is a relative number between 0 and 1. It specifies the range around the center of the fiber where this point is located. 

The actual location is draws from a uniform random distribution around the center, 

.. math::
  
  [0.5-\dfrac{s}{2}, 0.5+\dfrac{s}{2}),\\
  \text{with $s=$neuromuscularJunctionRelativeSize.}
  
The interval is multiplied by the number of points on the fiber, i.e. 0.5 indicates the center point. A value of 0 for `neuromuscularJunctionRelativeSize` indicates that the stimulation point is always at the center. A value of 0.1 indicates that the point is randomly at the center range of 10% of the fiber. Thus, for a lot of fibers, the position varies by maximum 10% fiber length.




