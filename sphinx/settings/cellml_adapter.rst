CellMLAdapter
==============

A CellMLAdapter class is used to run a CellML model.
It uses a code generator to produce efficient code for the given model at compile time.

The given CellML model is always computed for all nodes of a mesh, i.e. there are multiple instances being computed.
By specifying a mesh with 0 elements, you get a single instance of the model.

A CellML model is a first-order differential-algebraic system of equations (DAE) of the following form:

.. math::
   \frac{\partial \textbf{y}}{\partial t} = f(t,\textbf{y}(t),\textbf{h}(t),\hat{\textbf{c}},\hat{\textbf{p}}(t)) \\
   \textbf{h}(t) = g(\textbf{y}(t),\hat{\textbf{c}},\hat{\textbf{p}}(t))
   
The values :math:`\textbf{y} \in \mathbb{R}^n` are called *states* and will be integrated in time using a timestepping scheme. 
There are also the algebraic values, :math:`\textbf{h}`, which are not integrated. 
A set of parameters, :math:`\hat{\textbf{p}}`, can be defined in the settings and changed over time.
There are also constants , :math:`\hat{\textbf{c}}`, that are given in the CellML model and cannot be changed.

There exist several different names for the quantities :math:`\textbf{y}, \frac{\partial \textbf{y}}{\partial t}, \textbf{h}` and :math:`\hat{\textbf{p}}`:

=============================================== ================ =========== ========== ================== ==========================
symbol                                          opendihu         OpenCOR     OpenCMISS  computed by model  initial values can be set
=============================================== ================ =========== ========== ================== ==========================
:math:`\textbf{y}`                              states           states      STATES     by timestepping    yes
:math:`\frac{\partial \textbf{y}}{\partial t}`  rates            rates       RATES      yes                no
:math:`\textbf{h}`                              algebraics       algebraic   WANTED     yes                no
:math:`\hat{\textbf{c}}`                        constants        constants   CONSTANTS  no                 no
:math:`\hat{\textbf{p}}`                        parameters       algebraic   KNOWN      no                 yes
=============================================== ================ =========== ========== ================== ==========================
 
Initially the CellML model does not have any 'parameters', all values are given some defined value. 
In opendihu, any *constants* and *algebraics* can be transformed into *parameters* and then have changing values assigned.
This is done by the options ``parametersUsedAsAlgebraic`` and ``parametersUsedAsConstant``.

Usage
----------
The following C++ code shows the typical usage inside a time stepping scheme to solve the model:

.. code-block:: c

  TimeSteppingScheme::ExplicitEuler<
    CellmlAdapter<57,1>  // nStates,nAlgebraics: 57,71 = Shorten, 4,9 = Hodgkin Huxley
  >

The two template arguments of `CellmlAdapter` are the *number of states* and the *number of algebraics*.
This has to match the actual numbers of the CellML model that is to be computed. Consequently, when a specific model should be computed, the CellmlAdapter has be adjusted.

If the numbers are not correct a corresponding error will be shown from which the correct numbers can be determined.
  
Note that only explicit timestepping schemes are possible, which is current ``TimeSteppingScheme::ExplicitEuler`` or ``TimeSteppingScheme::Heun``.

There is an optional third template argument which specifies the function space, on which the CellML instances will be solved. 

.. code-block:: c

  TimeSteppingScheme::ExplicitEuler<
    CellmlAdapter<57,1,FunctionSpace<Mesh::StructuredDeformableOfDimension<1>,BasisFunction::LagrangeOfOrder<1>>>  // nStates,nAlgebraics: 57,71 = Shorten, 4,9 = Hodgkin Huxley
  >

This template argument is required if the mesh should be reused. 
E.g., for Monodomain eq. there is a splitting scheme with CellML and Diffusion and both parts use the same mesh. Then you have to assert that the mesh is the same type in the diffusion and here, e.g. by setting the mesh to structured deformable, as shown above.

The default FunctionSpace is `FunctionSpace::Generic` which is the following typedef:

.. code-block:: c

  typedef FunctionSpace<Mesh::StructuredRegularFixedOfDimension<1>,BasisFunction::LagrangeOfOrder<1>> Generic;

.. code-block:: python

  "CellML": {
    "modelFilename":                          "../../input/hodgkin_huxley_1952.c",    # CellML file (xml) or C++ source file
    #"libraryFilename":                       "cellml_simd_lib.so",                   # (optional) filename of a compiled library, overrides modelFilename
    #"statesInitialValues":                   [],                                     # (optional) initial values of all states, if not set, values from CellML model are used
    "initializeStatesToEquilibrium":          False,                                  # if the equilibrium values of the states should be computed before the simulation starts
    "initializeStatesToEquilibriumTimestepWidth": 1e-4,                               # if initializeStatesToEquilibrium is enable, the timestep width to use to solve the equilibrium equation
   
    # optimization parameters
    "optimizationType":                       "simd",                                 # "vc", "simd", "openmp" or "gpu": type of generated optimizated source file
    "approximateExponentialFunction":         True,                                   # if optimizationType is "vc" or "gpu", whether the exponential function exp(x) should be approximate by (1+x/n)^n with n=1024
    "compilerFlags":                          "-fPIC -O3 -march=native -shared ",     # compiler flags used to compile the optimized model code
    "maximumNumberOfThreads":                 0,                                      # if optimizationType is "openmp", the maximum number of threads to use. Default value 0 means no restriction.
    "useAoVSMemoryLayout":                    use_aovs_memory_layout,                 # if optimizationType is "vc", whether to use the Array-of-Vectorized-Stru    ct (AoVS) memory layout instead of the Struct-of-Vectorized-Array (SoVA) memory layout. Setting to True is faster.
    
    # stimulation callbacks
    #"setSpecificParametersFunction":         set_specific_parameters,                # callback function that sets parameters like stimulation current
    #"setSpecificParametersCallInterval":     int(1./stimulation_frequency/dt_0D),    # set_specific_parameters should be called every 1/stimulation_frequency seconds
    "setSpecificStatesFunction":              set_specific_states,                    # callback function that sets states like Vm, activation can be implemented by using this method and directly setting Vm values, or by using setSpecificParameters
    #"setSpecificStatesCallInterval":         2*int(1./stimulation_frequency/dt_0D),  # set_specific_states should be called stimulation_frequency times per ms, the factor 2 is needed because every Heun step includes two calls to rhs
    "setSpecificStatesCallInterval":          0,                                      # call intervall of the set_specific_states function, 0 means use setSpecificStatesCallFrequency instead
    "setSpecificStatesCallFrequency":         get_specific_states_call_frequency,     # set_specific_states should be called stimulation_frequency times per ms
    "setSpecificStatesFrequencyJitter":       get_specific_states_frequency_jitter,   # list of values to add or substract to setSpecificStatesCallFrequency every stimulation, this is to add random jitter to the frequency
    "setSpecificStatesRepeatAfterFirstCall":  0.01,                                   # [ms] simulation time span for which the setSpecificStates callback will be called after a call was triggered
    "setSpecificStatesCallEnableBegin":       get_specific_states_call_enable_begin,  # [ms] first time when to call setSpecificStates
    "additionalArgument":                     fiber_no,                               # any additional value that will be given to the callback functions
    
    
    "mappings": {                                                                     # mappings between parameters and algebraics/constants and between connectorSlots and states, algebraics or parameters
      ("parameter", 0):           ("constant", "membrane/i_Stim"),                    # parameter 0 is mapped to constant with name "membrane/i_Stim"
      ("connectorSlot", 0):       ("state", "membrane/V"),                            # as output connector slot 0 expose state with name "membrane/V"
    },
    
    #"algebraicsForTransfer":                 [],                                    # alternative way of specifying "mappings": which algebraic values to use in further computation
    #"statesForTransfer":                     [0],                                   # alternative way of specifying "mappings": which state values to use in further computation, Shorten / Hodgkin Huxley: state 0 = Vm
    #"parametersUsedAsAlgebraic":             [32],                                  # alternative way of specifying "mappings": list of algebraic value indices, that will be set by parameters. Explicitely defined parameters that will be copied to algebraics, this vector contains the indices of the algebraic array. This is ignored if the input is generated from OpenCMISS generated c code.
    #"parametersUsedAsConstant":              [65],                                  # alternative way of specifying "mappings": list of constant value indices, that will be set by parameters. This is ignored if the input is generated from OpenCMISS generated c code.
    "parametersInitialValues":                [0.0, 1.0],                            # initial values for the parameters, e.g. I_Stim, l_hs
    "meshName":                               "MeshFiber_{}".format(fiber_no),
    "stimulationLogFilename":                 "out/stimulation.log",
  },      
  
In the following all parameters will be explained.

modelFilename
---------------

This is the filename of the CellML model file. It can either be the XML file or a C/C++ code file. If it is an XML file, *opendihu* will use *OpenCOR* to convert it to a C source code file first.
Afterwards, *opendihu* will generate optimized C code (using the options given by the *optimization parameters*) and will store it as another file in the `src` subdirectory. The code will be compiled to a shared library (extension ’\*.so’) that will get loaded at runtime of the simulation. The shared library will be stored in the `lib` subdirectory.

libraryFilename
---------------

Optional, if given, it should be the filename of a shared object library (*.so) that will be used to compute the model.
This will be used instead of the model given in *modelFilename*. Usually this is only used to reuse library created by opendihu earlier.

statesInitialValues
---------------------
Optional. Default: `"CellML"`

If *statesInitialValues* is a list, it should contain an initial value for each state of the CellML model. 
If there are multiple instances all instances will be initialized by the same values.

If *statesInitialValues* is set to *CellML*, the initial values will be taken from the CellML model file (either XML or C). Usually this is what you want.

If *statesInitialValues* is set to *undefined*, no initial values will be set and the outer time stepping scheme can set initial values by giving `"initialValues"`.

initializeStatesToEquilibrium and initializeStatesToEquilibriumTimestepWidth
--------------------------------------------------------------------------------
If `initializeStatesToEquilibrium` is set to `True`, equilibrum values of the states in the CellML model will be computed before the simulation starts. Then, these values will be used to initialize the states.

Given the CellML model as

.. math::
   \frac{\partial \textbf{y}}{\partial t} = f(t,\textbf{y}(t),\textbf{h}(t),\hat{\textbf{c}},\hat{\textbf{p}}(t)),
   
the equation is solved by a 4th order Runge-Kutta timestepping scheme, until

.. math::
   \Vert\frac{\partial \textbf{y}}{\partial t}\Vert < \epsilon
   
is reached, with :math:`\epsilon = 1e-5`. The timestep width of the Runge-Kutta scheme can be given by `initializeStatesToEquilibriumTimestepWidth`. If an instability with this timestep width is detected (any value gets `inf` or `nan`), the timestep width will be decreased automatically and the computation will be restarted.

The resulting equilibrium values and the residuals are written to a file `<modelfilename>_equilibrium_values.txt`, where `<modelfilename>` is the file name of the model. An example for such a file is given below:

.. code-block:: c++

  // Result of computation of equilibrium values for the states by opendihu on 2020/2/29 10:17:12
  // Number of iterations: 10000000, dt: 0.0015625
  // Maximum ∂u/∂t = 0.0424747 for state 28
  // (If this is a high value, it indicates that the equilibrium was not fully reached.)

  state[0] = -81.0764;      // residuum: 3.15938e-05
  state[1] = -81.0242;      // residuum: 3.15353e-05
  state[2] = 7.25855;       // residuum: 5.68619e-06
  (...more lines follow...)
  state[53] = 0.00249843;   // residuum: 1.95519e-11
  state[54] = 0.213378;     // residuum: -6.67943e-07
  state[55] = 0.228239;     // residuum: -1.38375e-06
  state[56] = 2.8029e-10;   // residuum: -1.57379e-13

    Line to copy for settings:
    "statesInitialValues": [-81.0764, -81.0242, 7.25855, 150.928, 6.13908, 12.6374, 131.485, 132.853, 0.00809159, 0.995921, 0.0312117, 0.546801, 0.784615, 0.0081521, 0.995806, 0.0314177, 0.544509, 0.783771, 1.75163e-06, 5.90311e-06, 7.46021e-06, 4.19024e-06, 8.82585e-07, 0.875814, 0.118062, 0.00596817, 0.000134088, 1.12971e-06, -1580.24, 0.0284811, 53.9751, 0.0284799, 1687.43, 2.98746, 615, 615, 811, 811, 1283.85, 17808.2, 0.107779, 0.107778, 7243.03, 7243.03, 756.867, 756.867, 956.975, 956.975, 0.0343446, 0.0102602, 0.0136077, 0.0314302, 0.00312304, 0.00249843, 0.213378, 0.228239, 2.8029e-10],

The last line can be copy&pasted into the settings file and then specifies the initial values to be used in the next run.

Callbacks
-------------

A CellMLAdapter can have several callback functions. These are python functions that will be called in regular time intervals during the computation and can alter values of the computation.
They can be used, e.g., to stimulate a subcellular model at specific times.

The different callback functions and their time step interval by which the functions will be called are listed below. 
All of them will get the value of the option *additionalArgument* as its last argument. Like this it is possible to distinguish different instances in the functions when *CellMLAdapter* is nested inside *MultipleInstances*. This is the case for multiple fibers, where the *additionalArgument* can be the fiber number.

*setSpecificParametersFunction* and *setSpecificParametersCallInterval*
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Callback function and time step interval by which the function will be called.
This function can change some parameters and has the following signature:

.. code-block:: python

  def set_specific_parameters(n_dofs_global, timestep_no, current_time, global_parameters, additional_argument):
    # n_dofs_global:  (int) global number of dofs in the mesh, i.e. number of CellML instances to be computed
    # timestep_no:    (int)   current time step number, advances by the value of "setSpecificParametersCallInterval"
    # current_time:   (float) the current simulation time
    # global_parameters:  (dict)  initially an empty dict, the parameters to be changed should be indicated in this dict (see below)
    # additional_argument: The value of the option "additionalArgument", can be any Python object.
  
    # set parameters using calls like the following
    
    global_parameters{([x,y,z], nodal_dof_index, parameter_no)} = value
    # [x,y,z] are the global coordinates of the node to set the parameter
    # nodal_dof_index is the dof number of the node, usually 0. Only for Hermite ansatz functions it can be higher.
    # parameter_no is the parameter number to set 
    # value is the new parameter value

.. _callbackmesh:
.. figure:: images/callback_mesh.svg
  :width: 50%
  :align: center
  
  Example mesh with two subdomains and global natural ordering of the nodes.

For example, consider a mesh as in :numref:`callbackmesh` where a CellML model is computed on each node. The mesh is partitioned to two subdomains.
Rank 0 computes the grey nodes, rank 1 computes the blue nodes. The global natural ordering is given in the figure.

Then, on rank 0, ``dof_nos_global_natural`` will contain the list ``[0,1,4,5,8,9]`` and on rank 1, the list will be  ``[2,3,6,7,10,11]``. 
This shows to which global nodes the values in the `parameters` list correspond. With this information, the callback function could decide which parameters to update.

*setSpecificStatesFunction* and *setSpecificStatesCallInterval*
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Callback function and time step interval by which the function will be called.
This function can change some states and has the following signature:

.. code-block:: python

  def set_specific_states(n_dofs_global, timestep_no, current_time, global_states, additional_argument):
    # n_dofs_global:  (int) global number of dofs in the mesh, i.e. number of CellML instances to be computed
    # timestep_no:    (int)   current time step number, advances by the value of "setSpecificParametersCallInterval"
    # current_time:   (float) the current simulation time
    # global_states:  (dict)  initially an empty dict, the states to be changed should be indicated in this dict (see below)
    # additional_argument: The value of the option "additionalArgument", can be any Python object.
  
    # set states using calls like the following
    
    global_states{([x,y,z], nodal_dof_index, state_no)} = value
    # [x,y,z] are the global coordinates of the node for which to set the state
    # nodal_dof_index is the dof number of the node, usually 0. Only for Hermite ansatz functions it can be higher.
    # state_no is the state number to set 
    # value is the new state value

If ``setSpecificStatesFunction`` will be called, this happens during the time step update just before each evaluation of the right hand side / the CellML model.
I.e. for Heun's method it will be called up to twice per time step (depending on the other `setSpecificStates*` settings).
    
*setSpecificStatesCallEnableBegin*, *setSpecificStatesCallFrequency* and *setSpecificStatesFrequencyJitter*
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
If *setSpecificStatesCallInterval* is set to 0, the times when to call *setSpecificStatesFunction* are given by *setSpecificStatesCallEnableBegin*, *setSpecificStatesCallFrequency* and *setSpecificStatesFrequencyJitter*.

With these options, it is possible to efficiently specify a repeating pattern of calling the callback function. This is the recommended way to model a frequency encoded stimulation.

The first call of the callback is at simulation time *setSpecificStatesCallEnableBegin*. Using this parameter, a "ramp" can be modelled.
The callback is then called according to the frequency in *setSpecificStatesCallFrequency*. The frequency is :math:`1/T` and thus does not count timesteps, as with *setSpecificStatesCallInterval*, but uses the simulation time directly.

The frequency is modulated by applying a relative jitter, given in a list by *setSpecificStatesFrequencyJitter*. The jitter values are taken from the list and repeated. A value of 0 indicates no jitter, i.e. the frequency is met exactly. E.g., a value of 1.1 means a 10% longer time between subsequent calls to the function.

After the callback was called it will be repeated in the next timesteps *setSpecificStatesRepeatAfterFirstCall* times. Using this setting, a "square" signal can be modelled.

A visualization of the options is shown in :numref:`stimulation_times_2`.

.. _stimulation_times_2:
.. figure:: /settings/images/stimulation_times.svg
  :width: 80%
  
  Options that influence the stimulation. A time line is shown from left to right. The red blocks are time spans when `setSpecificStates` will be called. Because setSpecificStates usually checks a `firing times file` whether or not to activate the fiber, it can make sense to use the file `"MU_firing_times_always.txt"`. This file always indicates stimulation. Thus, the spike trains are completely determined by the options `setSpecificStatesCallEnableBegin`, `setSpecificStatesCallFrequency` and `setSpecificStatesFrequencyJitter`.
    
*handleResultFunction* and *handleResultCallInterval*
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Callback function and time step interval by which the function will be called.
This function can be used to postprocess the result and has the following signature:

.. code-block:: python

  def handle_result(n_instances, time_step_no, current_time, states_list, algebraics_list, name_information, additional_argument):
    # n_instances:         (int) local number of CellML instances to be computed
    # time_step_no:        (int)   current time step number, advances by the value of "setSpecificParametersCallInterval"
    # current_time:        (float) the current simulation time
    # states_list:         (list of floats) all local state values in struct-of-array memory layout,
    #                       i.e. [instance0state0, instance1state0, ... instanceNstate0, instance0state1, instance1state1, ...]
    # algebraics_list:  (list of floats) all local algebraic values in struct-of-array memory layout, 
    #                       i.e. [instance0algebraic0, instance1algebraic0, ... instanceNalgebraic0, instance0algebraic1, instance1algebraic1, ...]
    # name_information:    a map with the keys "stateNames" and "algebraicNames", contains lists of all CellML names of the states and algebraics
    # additional_argument: The value of the option "additionalArgument", can be any Python object.

    
    # asign some states to variables
    Vm = states[name_information["stateNames"].index("membrane/V") * n_instances]
    Ca_1 = states[name_information["stateNames"].index("razumova/Ca_1") * n_instances]
    A_1 = states[name_information["stateNames"].index("razumova/A_1") * n_instances]
    A_2 = states[name_information["stateNames"].index("razumova/A_2") * n_instances]
    x_1 = states[name_information["stateNames"].index("razumova/x_1") * n_instances]
    x_2 = states[name_information["stateNames"].index("razumova/x_2") * n_instances]
    
    # assign some algebraics to variables
    active_stress = algebraics[name_information["algebraicNames"].index("razumova/activestress") * n_instances]
    activation = algebraics[name_information["algebraicNames"].index("razumova/activation") * n_instances]
      
The example shows how one can access the state and algebraic variables by their name. The call to 

.. code-block:: python

  name_information["stateNames"].index("razumova/Ca_1")
  
gives the index of the state with the given name. Because the data for all locally computed instances is contained in the states array, we need to multiply this index with ``n_instances`` to get the first entry of the given state. This is now the index in ``states`` for the first instance. If the problem is monodomain on a fiber, in order to get the value at the center, use

.. code-block:: python

    Ca_1 = states[name_information["stateNames"].index("razumova/Ca_1") * n_instances + int(n_instances/2)]
      
How to specify mappings of states, algebraics and parameters
--------------------------------------------------------------------

The algebraics and constants in the CellML model can be replaced by so-called `parameters`. It is possible to define an arbitrary number of parameters (not completely arbitrary - the number has to be lower than the number of algebraics). These parameters act like constants during computation of the model. After each computation, their values can be changed either by callback functions or, if they are connected via an output slot to another solver, the values are set by the other solver.

The model to be computed appears as if the specified `algebraics` and `constants` had been replaced by the respective parameters.
This replacement relation is called `mapping` and can be defined in two different ways: the older way is by setting `parametersUsedAsAlgebraic` and `parametersUsedAsConstant`. The newer and recommended way is by using `mappings`.

Furthermore, some of the `states` and `algebraics` as well as some `parameters` can be connected to an output slot of the timestepping scheme and thereby reused by a different solver within a coupling or operator splitting scheme. Which `states`, `algebraics` and `parameters` to connect can again be specified in two ways: either by `algebraicsForTransfer` and `statesForTransfer` and `parametersForTransfer` or by `mappings`.

These settings will be explained in the following.

parametersUsedAsAlgebraic
^^^^^^^^^^^^^^^^^^^^^^^^^^^^
(list of int) List of algebraic numbers that will be replaced by parameters.
There are explicitely defined parameter values that will be copied to these algebraics. 
This vector contains the indices of the algebraic array. 
Note, that these values can also be set by the ``mappings`` option, which is more clear.

parametersUsedAsConstant
^^^^^^^^^^^^^^^^^^^^^^^^^^^^
(list of int) List of indices, which constants in the computation will be replaced by parameters.
Note, that these values can also be set by the ``mappings`` option, which is more clear.

*algebraicsForTransfer* and *statesForTransfer*
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
(list of ints) Which algebraics and states should be transferred to the other solver in either a `Coupling`, `GodunovSplitting` or `StrangSplitting`.

The total number of field variables to be transferred is the sum of the length of these two settings (+number of parameters if specified).

Note, that these values can also be set by the ``mappings`` option, which is more clear.

parametersInitialValues
---------------------------
(list of float) List of values of the parameters. This also defines the number of parameters.

Example:

.. code-block:: python

  parametersInitialValues = [1.0, 2.0, 3.0]
  parametersUsedAsAlgebraic = [5, 2]
  parametersUsedAsConstant[10]
  
This example will compute the given CellML model with the following modifications: The algebraic/algebraic values ``algebraics[5]`` and ``algebraics[2]`` will not be computed by the model, but get the values ``1.0`` and ``2.0``. These values may be changed later using one of the callback functions.
The variable ``constants[10]`` will be set to ``3.0`` and not changed.
  
mappings
-------------
(dict)
Under ``mapping`` it is possible to specify the connection of `parameters` to `algebraics` and `constants`, 
as well as the connection of `connectorSlots` to `states`, `algebraics` and `parameters`. An example is given below (the actual names are only dummies and make no sense):
  
.. code-block:: python

  "mappings" : {
      ("parameter", 0):           ("algebraic", "wal_environment/I_HH"),
      ("parameter", 1):           ("constant", "razumova/L_x"),
      
      ("connectorSlot", 0):       ("state", "wal_environment/vS"),
      ("connectorSlot", 1):       ("state", 5),  
      ("connectorSlot", 2):       ("state", "potassium_channel_n_gate/n"),
      ("connectorSlot", 2):       "potassium_channel_n_gate/n",             # alternative
      ("connectorSlot", 3, "A"):  ("algebraic", "leakage_current/i_L"),
      ("connectorSlot", 3, "A"):  "leakage_current/i_L",                    # alternative
      ("connectorSlot", "slotB"): ("parameter", 0),
      ("connectorSlot", "lambda"):("constant", "razumova/L_S"),         # expose fiber stretch to get the current fiber stretch from the mechanics solver
    }
    
The value of `mappings` is a Python Dict. 
Each key (left hand side) has one of the following formats:

* ``("parameter", 0)`` to specify a parameter with given number. The number is needed to identify the initial values for the parameters.
* ``("connectorSlot", 0)`` where ``0`` can be any integer number, to specify a connector slot, the number is arbitrary and is only used to order multiple slots in relation to each other.
* ``("connectorSlot", "slotA")`` here with a slot name, the name has to be maximum 10 characters long.
* ``("connectorSlot", 0, "slotA")`` This is a combination of the two formats above, it specifies a slot name and also a number for ordering the slots.

The value that corresponds to the key (right hand side) of one `mappings` item is a two-element tuple or string of the form 

* ``("name", "cellml name")``
* or ``("name", 0)``
* or ``"cellml name"``,

where ``"name"`` has to be either ``"constant"``, ``"state"``, ``"algebraic"`` or ``"parameter"``. The ``"cellml name"`` is the name of the variable in the CellML model in the form ``"componentName/variableName"`` and ``0`` can be any valid index. This means, it is possible to identify, e.g. a state by its name as well as by its index in the C code file.
If there is no tuple but only the "cellml name", it will be determine automatically if it is a `state`, `algebraic` or `constant` by searching among all available cellml names.

For the parameters, the index must start with `0` and increase by one for all further parameters. As already mentioned, the mapped variable for a parameter can be an `"algebraic"` or a `"constant"`. The beginning of the parameters list must all map to algebraics and the rest must map to constants. I.e., every constant must be mapped to a parameter with lower index than all the parameters that are mapped to algebraics. The specified mappings will internally be transferred to the ``parametersUsedAsAlgebraic`` and ``parametersUsedAsConstant`` lists that can otherwise also be set directly by these options.

Also for the `"connectorSlots"` there is a required order. At first, all mapped `"states"` have to be given, then all `"algebraics"` and then all `"parameters"`. 

Note that the values of parameters will not be changed by the CellML model. If you need to reuse values computed within the CellML model, use states or algebraics. The purpose of connecting parameters to output slots is to allow the initial parameter value to be set by a different solver.

Typical mappings and initial values of parameters by commonly used cellml models (in variable ``cellml_file``) are given below. Note that these do not set slot names. But for more complex examples it would be good to add slot names.

.. code-block:: python

  # set variable mappings for cellml model
  if "hodgkin_huxley" in cellml_file:
    # parameters: I_stim
    mappings = {
      ("parameter", 0):     ("constant", "membrane/i_Stim"),      # parameter 0 is constant 2 = I_stim
      ("connectorSlot", 0): ("state", "membrane/V"),              # expose state 0 = Vm to the operator splitting
    }
    parameters_initial_values = [0.0]                         # initial value for stimulation current
    
  elif "shorten" in cellml_file:
    # parameters: stimulation current I_stim, fiber stretch λ
    mappings = {
      ("parameter", 0):     ("algebraic", "wal_environment/I_HH"), # parameter is algebraic 32
      ("parameter", 1):     ("constant", "razumova/L_x"),             # parameter is constant 65, fiber stretch λ, this indicates how much the fiber has stretched, 1 means no extension
      ("connectorSlot", 0): ("state", "wal_environment/vS"),          # expose state 0 = Vm to the operator splitting
    }
    parameters_initial_values = [0.0, 1.0]                        # stimulation current I_stim, fiber stretch λ
    
  elif "slow_TK_2014" in cellml_file:   # this is (3a, "MultiPhysStrain", old tomo mechanics) in OpenCMISS
    # parameters: I_stim, fiber stretch λ
    mappings = {
      ("parameter", 0):     ("constant", "wal_environment/I_HH"), # parameter 0 is constant 54 = I_stim
      ("parameter", 1):     ("constant", "razumova/L_S"),         # parameter 1 is constant 67 = fiber stretch λ
      ("connectorSlot", 0): ("state", "wal_environment/vS"),      # expose state 0 = Vm to the operator splitting
      ("connectorSlot", 1): ("algebraic", "razumova/stress"),     # expose algebraic 12 = γ to the operator splitting
      ("connectorSlot", "lambda"):("constant", "razumova/L_S"),   # expose fiber stretch to get the current fiber stretch from the mechanics solver
    }
    parameters_initial_values = [0.0, 1.0]                    # wal_environment/I_HH = I_stim, razumova/L_S = λ
    
  elif "Aliev_Panfilov_Razumova_2016_08_22" in cellml_file :   # this is (3, "MultiPhysStrain", numerically more stable) in OpenCMISS, this only computes A1,A2,x1,x2 not the stress
    # parameters: I_stim, fiber stretch λ, fiber contraction velocity \dot{λ}
    mappings = {
      ("parameter", 0):     ("constant", "Aliev_Panfilov/I_HH"),  # parameter 0 is constant 0 = I_stim
      ("parameter", 1):     ("constant", "Razumova/l_hs"),        # parameter 1 is constant 8 = fiber stretch λ
      ("parameter", 2):     ("constant", "Razumova/velo"),        # parameter 2 is constant 9 = fiber contraction velocity \dot{λ}
      ("connectorSlot", 0): ("state", "Aliev_Panfilov/V_m"),      # expose state 0 = Vm to the operator splitting
      ("connectorSlot", 1): ("algebraic", "Razumova/sigma"),      # expose algebraic 0 = γ to the operator splitting
    }
    parameters_initial_values = [0, 1, 0]                     # Aliev_Panfilov/I_HH = I_stim, Razumova/l_hs = λ, Razumova/velo = \dot{λ}
    
  elif "Aliev_Panfilov_Razumova_Titin" in cellml_file:   # this is (4, "Titin") in OpenCMISS
    # parameters: I_stim, fiber stretch λ, fiber contraction velocity \dot{λ}
    mappings = {
      ("parameter", 0):     ("constant", "Aliev_Panfilov/I_HH"),  # parameter 0 is constant 0 = I_stim
      ("parameter", 1):     ("constant", "Razumova/l_hs"),        # parameter 1 is constant 11 = fiber stretch λ
      ("parameter", 2):     ("constant", "Razumova/rel_velo"),    # parameter 2 is constant 12 = fiber contraction velocity \dot{λ}
      ("connectorSlot", 0): ("state", "Aliev_Panfilov/V_m"),      # expose state 0 = Vm to the operator splitting
      ("connectorSlot", 1): ("algebraic", "Razumova/ActiveStress"),   # expose algebraic 4 = γ to the operator splitting
      ("connectorSlot", 2): ("algebraic", "Razumova/Activation"),     # expose algebraic 5 = α to the operator splitting
    }
    parameters_initial_values = [0, 1, 0]                     # Aliev_Panfilov/I_HH = I_stim, Razumova/l_hs = λ, Razumova/rel_velo = \dot{λ}
    

meshName
------------------------------------------------
The mesh to use, to be defined under "Meshes". For details, see :ref:`define_meshes`. You can instead also just specify ``nElements`` to directly set the number of instances to be computed.

If no mesh is specified at all, the standard is ``"nElements": 0``. This corresponds to 1 node, i.e. one instance of the CellML problem. There will be the warning about the missing *nElements* though.

stimulationLogFilename
------------------------------------------------
Default: "out/stimulation.log"

A file name of an output file that will contain all firing times.

optimizationType
--------------------
Possible values: ``simd``, ``vc``, ``openmp`` or ``gpu``. Which type of code to generate. ``openmp`` produces code for shared-memory parallelization, using OpenMP. ``simd`` produces auto-vectorizable code. ``vc`` produces explicitly vectorized code (fastest). ``gpu`` is only available if the :doc:`fast_monodomain_solver` is used.

See also the notes on ``vc`` about AVX-512 on the page of :doc:`fast_monodomain_solver`.

compilerFlags
-----------------
Additional compiler flags for the compilation of the source file. Default: ``-fPIC -finstrument-functions -ftree-vectorize -fopt-info-vec-optimized=vectorizer_optimized.log -shared``

When compiled in release target, ``-O3`` is added. In debug target, ``-O0 -ggdb`` is added. If *optimizationType* is ``openmp``, ``-fopenmp`` is added.

