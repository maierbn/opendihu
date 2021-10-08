MapDofs
===============

This is a class that copies and transforms values from one mesh to another in a configurable way. 
With parallel computation this also involves communicating the values between the processes.

It is a wrapper to another solver. It adds a number of field variables which can be defined on an own function space. The slots of `MapDofs` consists of the slots of the nested solver, as well as the additional field variables.
The value transfer can be done between any two of these slots, the `from` or `source` slot and the `to` or `target` slot. With advanced callbacks it is also possible to map to multiple output slots at once.

There are several modes, what exactly is done in this transfer:

  * There is a simple local transfer, where values are copied from one slot to another (mode `"copyLocal"`).
  * There is a conditional transfer, where only positive values are copied (mode `"copyLocalIfPositive"`).
  * Instead of copying from the first slot, one can define the value to be set in the second slot if the value of the first slot is above a defined threshold (mode `"localSetIfAboveThreshold"`).
  * The previous modes are simple actions that may occur frequently and thus have been implemented in the C++ core. However, other use cases can be programmed in python.
    A Python callback function can handle the transfer and involve any computation. This allows to e.g. delaying a signal or to map from :math:`M` inputs to :math:`N` outputs (mode `"callback"`).
  * Whereas the previous modes all work entirely on the local domain, there is the option to communicate between processes. Local or global degree of freedom (dof) numbers can be defined for specifying the input and output values.

C++ code
^^^^^^^^^^

The first template argument is a function space. The additional field variables will be defined with this function space.
The second template argument is the nested solver.
  
.. code-block:: c

  Control::MapDofs<
    // function space for additional field variables
    ,
    // nested solver
  >
  
For example:

.. code-block:: c

  Control::MapDofs<
    FunctionSpace::FunctionSpace<Mesh::StructuredRegularFixedOfDimension<1>,BasisFunction::LagrangeOfOrder<1>>,
    OperatorSplitting::Strang<
      ...
    >
  >

Python settings
^^^^^^^^^^^^^^^^^

The follows shows all python settings. 

.. code-block:: python

  "MapDofs": {
    "nAdditionalFieldVariables":  1,                              # number of additional field variables that are defined by this object. They have 1 component, use the templated function space and mesh given by meshName.
    "additionalSlotNames":        ["mn"],                         # list of names of the slots for the additional field variables
    "meshName":                   "motoneuronMesh",               # the mesh on which the additional field variables will be defined
    
    # mapping from motoneuronMesh which contains on every rank as many nodes as there are motoneurons to the 3D domain
    # map from motoneuronMesh (algebraics) to 3Dmesh (solution)
    "beforeComputation": [                                        # transfer/mapping of dofs that will be performed before the computation of the nested solver, can be None if not needed
      {                                                 
        # alternative forms:
        "fromConnectorSlot":                2,                    # The slot no or slot name from which the data will be taken.
        "fromConnectorSlot":                "slotA",              # The slot no or slot name from which the data will be taken.
        
        # alternative forms:
        "toConnectorSlots":                 0,                    # The slots to which the data will be written. This can be a list if multiple slots are needed. Either as slot no.s or names.
        "toConnectorSlots":                 [0,1],                # The slots to which the data will be written. This can be a list if multiple slots are needed. Either as slot no.s or names.
        "toConnectorSlots":                 "slotB",              # The slots to which the data will be written. This can be a list if multiple slots are needed. Either as slot no.s or names.
        "toConnectorSlots":                 ["slotB","slotC"],    # The slots to which the data will be written. This can be a list if multiple slots are needed. Either as slot no.s or names.
        
        "fromOutputConnectorArrayIndex":    0,                    # which fiber/compartment for the input slot, if there are multiple
        "toOutputConnectorArrayIndex":      0,                    # which fiber/compartment for the output slot, if there are multiple      
        "mode":                             "callback",           # "copyLocal", "copyLocalIfPositive", "localSetIfAboveThreshold", "callback" or "communicate"
        "fromDofNosNumbering":              "local",              # if the input dofs are interpreted as local or global numbers, refers to "dofsMapping" and "inputDofs" 
        "toDofNosNumbering":                "global",             # if the output dofs are interpreted as local or global numbers, refers to "dofsMapping" and "outputDofs"
        "dofsMapping":                      {0: [n_elements/2-1, n_elements/2, n_elements/2+1]},    # a mapping from input to output dofs, for all modes except "callback"
        "inputDofs":                        0,                    # the input dofs for the mode "callback"
        "outputDofs":                       [n_elements/2-1, n_elements/2, n_elements/2+1],       # the output dofs for the mode "callback", if multiple output dofs are used, this can be a list of lists of the dofs
        "callback":                         callback_motoneuron,  # the callback function if the mode is "callback"
        #"thresholdValue":                   20,                  # if mode is "localSetIfAboveThreshold", this is the threshold, if the value is above it, set the value `valueToSet`
        #"valueToSet":                       20,                  # if mode is "localSetIfAboveThreshold", this is the value to set the target dof to, if the source dof is above thresholdValue.
      }
    ],
    "afterComputation":             None,                         # transfer/mapping of dofs that will be performed after the computation of the nested solver, can be None if not needed
    
    # Nested solver, e.g.
    "StrangSplitting": {
      ...
    }
  }

`nAdditionalFieldVariables`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Number of field variables that will be defined on the additional function space that was given as first template argument in the c++ code. These can be useful to map from the nested solvers' slots to some other mesh.

`additionalSlotNames`
^^^^^^^^^^^^^^^^^^^^^^^
A list of strings with as many entries as `nAdditionalFieldVariables`. The list should contains the names for the slots that will be created for the additional field variables. Each name has to consist of 6 or less characters.

`meshName`
^^^^^^^^^^^^^^^^^
Specification or reference of the mesh (see :doc:`mesh` for details how to specify meshes inline and under the `"Meshes"` key) to be used for the additional field variables. The type of the mesh is as given by the first template argument in the c++ code.
    
`beforeComputation` and `afterComputation`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Under these keys, a list of mapping actions can be defined that will be performed directly before and after the solution step of the nested solver.
If no actions for one of them is needed, `None` can be used. Each action is given as a dict containing more specific options to the action.

`fromConnectorSlot` and `toConnectorSlots`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
These two values specify the slots for the value transfer. The slots can either be identified by the slot number or by the slot name.

* If the slots should be used to identify the slots, consult the `solver_structure.txt` representation to find out the slot numbers. Then, option `fromConnectorSlot` is simply a single integer. The option `toConnectorSlots` is either a single number or a list of multiple numbers for multiple output slots. Multiple output slots are only supported for the mode "callback".
* If the slots should be identified by slot names, the same holds, but instead of numbers, simply provide the slot names as strings.

`fromOutputConnectorArrayIndex` and `toOutputConnectorArrayIndex`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
If the slots contain multiple instances of the actual slot, these two options specify which of the instance to use for the mapping. 
This occurs, if a nested solver contains `MultipleInstances` somewhere. For example for multiple fibers or multiple compartments for multidomain.

Note that `toOutputConnectorArrayIndex` is always a single number, even if multiple output slots are used for `toConnectorSlotNo`.

`mode`
^^^^^^^^
One of "copyLocal", "copyLocalIfPositive", "localSetIfAboveThreshold", "callback" or "communicate", specifies what to do for the transfer.

* `copyLocal`: Copy dofs within the local domain as specified in `"dofsMapping"`, dofs on other processes are ignored.
* `copyLocalIfPositive`: Same as `copyLocal`, but the target value is only set if the source value is positive.
* `localSetIfAboveThreshold`: Similar to `copyLocalIfPositive`, but the threshold value can be customized by the option `"thresholdValue"`. Instead of copying the source dof, a fixed value given by `"valueToSet"` is used.
* `communicate`: Perform the mapping specified in `"dofsMapping"` and also consider dofs on remote processes. The source dofs can be given as either local or global numbers. The target dofs have to be given as global numbers, i.e. `"toDofNosNumbering"` has to be `"global"`.
* `callback`: Do not use the `"dofsMapping"`, instead specify what to map by a custom callback function. The function is provided in `"callback"`, see below for the signature. The input dofs and output dofs are given by `"inputDofs"` and `"outputDofs"` and can both be specified in local or global numbering. Again, only the locally present dofs are considered. If you need the callback plus global communication, use two actions, one with mode "communicate" and one with "callback".

  The `callback` option is the only one to allow to map to multiple output slots. If this is needed, the option `"toConnectorSlotNo"` is a list (of lists) with entries for the different slots.

Depending on the mode, other options have to be given.
All modes need the options `"fromDofNosNumbering"` and `"toDofNosNumbering"`. These specify if dof numbers for the source and target slots are specified in *local numbering* or *global numbering*.
For the modes *"copyLocal"*, *"copyLocalIfPositive"*, *"localSetIfAboveThreshold"* and *"communicate"*, the additional option `"inputDofs"` is needed.

For mode *"localSetIfAboveThreshold"*, additionally, `"thresholdValue"` and `"valueToSet"` are needed.

For the mode *"callback"*, the options *"inputDofs"*, *"outputDofs"* and *"callback"* need to be given, instead of `dofsMapping`.

`"fromDofNosNumbering"` and `"toDofNosNumbering"`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
One of "local", "global". Specifies if the dof numbers given as `key:value` pairs in the dict `dofsMapping` are interpreted as *local numbering* or *global numbering*. (`fromDofNosNumbering` refers to key, `toDofNosNumbering` refers to value dofs).
For the `callback` mode, it refers to the numbers in `inputDofs` (`fromDofNosNumbering`) and `outputDofs  (`fromDofNosNumbering`).

`dofsMapping`
^^^^^^^^^^^^^^^^^

Specification of which dofs values at the "from" slot will be transferred to which dofs at the "to" slot. It is a dict such as,
e.g., ``{0: 1, 2: [5,6,8], 3: 0}``. This would copy the value at dof 0 to the other slot at dof 1,
the value at dof 2 to the other slot at three values at once (5,6,8) and dof 3 to dof 0.

The dof numbers are interpreted either as local or global numbers, depending on the valeu of `"fromDofNosNumbering"` and `"toDofNosNumbering"`. 
Global numbers that are not present on the own process are ignored, for both the source and the target dofs.

`inputDofs`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If `"mode"` is `"callback"`, this is a list of dofs for which the values are taken from the field variable at the input slot.

`outputDofs`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If `"mode"` is `"callback"`, this is a list of lists of dofs for which the values are taken from the field variable for each output slot.
If only one output slot is needed, it has to be of the form ```[[0,1,2,...]]```, i.e. also as list of lists.

`callback`
^^^^^^^^^^^^^

A python function that performs the mapping between a potentially different number of input and output dofs, used for mode "callback".
 The function has the following form. An example is given that delays the input signal number 0, and writes a gaussian stimulus with maximum value of 20 to all output dofs.

.. code-block:: python

  def callback_motoneuron(input_values, output_values, current_time, slot_nos, buffer):
    """
    Callback function that transform a number of input_values to a number of output_values.
    This function gets called by a MapDofs object.
    :param input_values: (list of float values) The input values from the slot as defined in the MapDofs settings.
    :param output_values: (list of list of float values) output_values[slotIndex][valueIndex]
                          The output values buffer, potentially for multiple slots.
                          Initially, this is a list of the form [[None, None, ..., None]] with the size matching 
                          the number of required output values. The function should set some of the entries to a computed value.
                          The entries that are not None will be set in the output slot at the dofs defined by MapDofs.
    :param current_time:  Current simulation time.
    :param slot_nos:      List of [fromSlotNo, toSlotNo, fromArrayIndex, toArrayIndex].
    :param buffer:        A persistent helper buffer. This variable can be set to anything and will be provided back to 
                          this function every time. Using this buffer, it is possible to implement a time delay of signals.
    """
      
    # get number of input and output values
    n_input_values = len(input_values)      
    n_output_values = len(output_values[0])     # number of output values for the first output slot
    
    # initialize buffer the first time
    if 0 not in buffer:
      buffer[0] = None
    
    # determine spike by threshold
    if input_values[0] > 20:
      buffer[0] = current_time    # store time of last activation in buffer[0]
      
    # if there has been a stimulation so far
    if buffer[0] is not None:
      
      # convolute Dirac delta, kernel is a shifted and scaled gaussian
      t_delay = 10              # [ms] delay of the signal
      gaussian_std_dev = 0.1    # [ms] width of the gaussian curve
      convolution_kernel = lambda t: np.exp(-0.5 * ((t - t_delay) / gaussian_std_dev)**2)
      delayed_signal = convolution_kernel(current_time - buffer[0]) * 20
        
      # loop over output values and set all to the computed signal, cut off at 1e-5
      if delayed_signal > 1e-5:
        print("motoneuron t: {}, last_activation: {}, computed delayed_signal: {}".format(current_time, buffer[0], delayed_signal))
        for i in range(n_output_values):
          output_values[0][i] = delayed_signal
      else:
        for i in range(n_output_values):
          output_values[0][i] = None     # do not set any values
      
Exemplary solver structure
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The following is a solver structure that uses a MapDofs. The actions are indicated by the two arrows with double tips.

.. code-block:: bash

  Solver structure: 

  ├── Coupling                                              
  │  output slots:                                          
  │  [a] solution.membrane/V                     +────── ¤0 x
  │  [a] solution                                :+───── ¤1 x
  │  [b] additionalFieldVariable0                ::+──── ¤2 x
  │  [b] additionalFieldVariable1                :::+─── ¤3 x
  │                                              ::::       
  │  slot connections:                           ::::       
  │  1¤ <═> ¤2                                   ::::       
  │  2¤ <─> ¤3                                   ::::       
  │                                              ::::       
  │ ├── Heun                                     ::::       
  │ │   ("Term1")                                ::::       
  │ │  output slots:                             ::::       
  │ │  [b] solution.membrane/V                   +÷÷÷─── ¤0 x
  │ │  [b] firing_threshold/V_extern_out (in var  +÷÷─── ¤1════╗
  │ │  [b] (P)firing_threshold/V_extern_in (in v   +÷─── ¤2<─┐ ║
  │ │                                               :        │ ║
  │ │ └── CellmlAdapter                             :        │ ║
  │ └                                               :        │ ║
  │                                                 :        │ ║
  │ ├── MapDofs                                     :        │ ║
  │ │   ("Term2")                                   :        │ ║
  │ │  output slots:                                :        │ ║
  │ │  [a] solution.membrane/V              ┌»┌     +─── ¤0 x│ ║
  │ │  [a] solution                         │ │     :+── ¤1 x│ ║
  │ │  [b] additionalFieldVariable0         └ │     ::   ¤2══┼═╝
  │ │  [b] additionalFieldVariable1           └»    ::   ¤3<─┘
  │ │                                               ::      
  │ │ ├── StrangSplitting                           ::      
  │ │ │  output slots:                              ::      
  │ │ │  [a] solution.membrane/V                    +÷── ¤0 x
  │ │ │  [a] solution                               :+── ¤1 x
  │ │ │                                             ::      
  │ │ │  slot connections:                          ::      
  │ │ │  0¤ <═> ¤0                                  ::      
  │ │ │  1¤ <═> ¤1                                  ::      
  │ │ │  2¤ <═> ¤2                                  ::      
  │ │ │                                             ::      
  │ │ │ ├── Heun                                    ::      
  │ │ │ │   ("Term1")                               ::      
  │ │ │ │  output slots:                            ::      
  │ │ │ │  [a] solution.membrane/V                  +÷── ¤0══╗
  │ │ │ │                                            :       ║
  │ │ │ │ └── CellmlAdapter                          :       ║
  │ │ │ └                                            :       ║
  │ │ │                                              :       ║
  │ │ │ ├── CrankNicolson                            :       ║
  │ │ │ │   ("Term2")                                :       ║
  │ │ │ │  output slots:                             :       ║
  │ │ │ │  [a] solution                              +── ¤0══╝
  │ │ │ │                                                   
  │ │ │ │ ├── FiniteElementMethod                           
  │ │ │ │ │  output slots:                                  
  │ │ │ │ │  [a] solution                                ¤0 x
  │ │ │ │ │                                                 
  │ │ │ └                                                   
  │ │ └                                                     
  │ └                                                       
  └                                                         
                                                            
  Connection Types:
    +··+   Internal connection, no copy
    ════   Reuse variable, no copy
    ───>   Copy data in direction of arrow
    ─m──   Mapping between different meshes

  Referenced Meshes:
    [a] "MeshFiber", 1D regular fixed, linear Lagrange basis
    [b] "motoneuronMesh", 1D regular fixed, linear Lagrange basis



  
