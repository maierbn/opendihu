PrescribedValues
=================

The `PrescribedValues` allows to directly set the values of a field variable over time from a python callback function.
It is a normal timestepping scheme. The "solution" step consists of calling a python callback for every field variable. 
The field variable can be 1D, 2D or 3D, with any mesh and ansatz function.

Two different types of field variables with different numbers of components can be specified. Then, of each type, any number of field variables can be created. E.g. if scalar and two-component vector-valued field variables should be used, there can be any number of scalar field variables plus any number of vector-valued field variables.

The first component of every field variable is automatically connected to the output slot of the `PrescribedValues` class (see :doc:`output_connector_slots`).

Usage
----------
The PrescribedValues class has three template parameters, of which the last two are optional (with default value 1).

.. code-block:: c

  PrescribedValues</*FunctionSpace*/, /*nComponents1*/, /*nComponents2*/>
  
for example: 

.. code-block:: c

  PrescribedValues<
    FunctionSpace::FunctionSpace<
      Mesh::StructuredRegularFixedOfDimension<2>,
      BasisFunction::LagrangeOfOrder<2>
    >,
    1,
    3
  >
  
The above example uses a 2D structured regular fixed mesh with quadratic Lagrange ansatz functions. 
It allows to define different field variables with one and three components, i.e. scalar field variables and field variables with 3 components. If the last two numbers are not specified, they are both 1.

The python settings are as given below:

.. code-block:: python

  "PrescribedValues": {
    "meshName":            "mesh",        # the mesh to use for the field variables
    "numberTimeSteps":     1,             # number of timesteps to call the callback functions subsequently, this is usually 1 for prescribed values, because it is enough to set the reaction term only once per time step
    "timeStepOutputInterval": 10,         # if the time step should be written to console, a value > 1 produces no output
    "slotNames":           ["reac"],      # list of names of the connector slots, maximum length is 10 characters per name
    
    # a list of field variables that will get values assigned in every timestep, by the provided callback function
    "fieldVariables1": [
      {"name": "reaction_term", "callback": set_reaction_term_values},
    ],
    "fieldVariables2":     [],
    "additionalArgument":  None,          # a custom argument to the fieldVariables callback functions, this will be passed on as the last argument
    
    "OutputWriter" : [
      #{"format": "Paraview", "interval": 1, "filename": "out", "binaryOutput": "false", "fixedFormat": False, "onlyNodalValues":True, "fileNumbering": "incremental"},
      {"format": "PythonFile", "outputInterval": 100, "filename": "out/reaction", "binary": True, "onlyNodalValues":True, "fileNumbering": "incremental"}
    ]
  },      

In the following all parameters will be explained.

meshName
----------
The :doc:`mesh` on which the field variable will live. It is also possible to specify the mesh in-line, as explained for :ref:`Finite Element Method <femesh>`.

slotNames
----------
A list of strings, names for the connector slots. Each name should be smaller or equal than 10 characters. 
In general, named slots are used to connect the slots from a global setting "connectedSlots". See :doc:`output_connector_slots` for details.

additionalArgument
----------------------
This is a custom argument that will be passed to the fieldVariables callback functions as the last argument.


Callback Function
-------------------

The callback functions that are specified under "callback" have the following signatures. The function will provide the current local values of the field variable in the ``values`` parameter, as well as some other information.
The values in this list can be changed inside the callback function and are then updated in the field variable. Only the local values can be accessed inside the callback function. But there is enough information passed to the callback to conveniently access the correct values.

The following example also shows, how one can iterate over the local values.

.. code-block:: python

  def set_reaction_term_values(n_dofs_global, n_nodes_global_per_coordinate_direction, time_step_no, current_time, values, global_natural_dofs, additional_argument):
    # n_dofs_global:       (int) global number of dofs in the mesh where to set the values
    # n_nodes_global_per_coordinate_direction (list of ints)   [mx, my, mz] number of global nodes in each coordinate direction. 
    #                       For composite meshes, the values are only for the first submesh, for other meshes sum(...) equals n_dofs_global
    # time_step_no:        (int)   current time step number
    # current_time:        (float) the current simulation time
    # values:              (list of floats) all current local values of the field variable, if there are multiple components, they are stored in struct-of-array memory layout 
    #                       i.e. [point0_component0, point0_component1, ... point0_componentN, point1_component0, point1_component1, ...]
    #                       After the call, these values will be assigned to the field variable.
    # global_natural_dofs  (list of ints) for every local dof no. the dof no. in global natural ordering
    # additional_argument: The value of the option "additionalArgument", can be any Python object.
    
    for local_dof_no in range(len(values)):
      # get the global no. of the current dof
      global_dof_no = global_natural_dofs[local_dof_no]
      
      i = global_dof_no % n_nodes_global_per_coordinate_direction[0]        # index in x direction
      j = int(global_dof_no / n_nodes_global_per_coordinate_direction[0])   # index in y direction
      
      # now we know that values[local_dof_no] is the value at node (i,j) in the global mesh
      
      # e.g. set gaussian
      center = np.array((2,4))
      x = np.linalg.norm(np.array((i,j)) - center)
      values[local_dof_no] = scipy.stats.norm.pdf(x)
      
Dummy
=================

The `Dummy` class can be used as a placeholder, e.g. in a Coupling scheme when the second Term is not yet implemented but the first one should already be tested.

It is used without any template argument:

.. code-block:: c

  Dummy
  
It does not have any Python settings. Also there are no field variables and no connector slots.
