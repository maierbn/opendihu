MultipleInstances
==================

MultipleInstances is a class that creates multiple instances out of the underlying solver. A use case is, for example, when we want to simulate multiple muscle fibers. We only need to setup the solver structure for a single fiber and can wrap it by a MultipleInstances class. In the settings we can specify how many fibers we want to have and then for each instance or fiber the individual settings. Like this we can specify e.g. different geometry for the fibers in the Python settings and use the same C++ solver code at the same time.

The MultipleInstances class simply creates the specified number of instances and runs them one after another in a for loop.
In the mentioned example with the fibers, this would mean that all fibers will get computed one ofter another. But, usually we want to compute the equations of the fibers synchronously, from timestep to timestep for all fibers at once. This is also possible, for this we have to put the MultipleInstances class into the inner solver. For details, study the C++ file of the `fibers_emg` example.

The C++ code is simply:

.. code-block:: c

  Control::MultipleInstances<
    // nested solvers
  >
  
The settings have the following form:

.. code-block:: python

  "MultipleInstances": {
    "nInstances": nInstances,
    "instances": [
      {
        "ranks": [0,1,2],
        
        # further solver property objects, e.g. for GodunovSplitting
        "GodunovSplitting": {
        }
      },
      {
      #... further instance
      }
    ]
    "OutputWriter" : [...],
  }
  
nInstances
------------
The number of instance to create and run in total. 

instances
------------
A list of settings for the instances. The size has to be at least `nInstances`. If the list is larger, exceeding entries will be ignored.
Each instance takes the next settings from the list. 

One settings instance is a Python dictionary with an item ``"ranks"`` (see below) and the settings for the instance.

ranks
--------
A list of MPI ranks on which the instance should be computed. Since this property is part of the settings for a specific instance, it specifies the processes that will be used for the particular instance. This list can be different for different instances.

A process only takes part in the computation of those instances which have the rank number of the process included in `ranks`. Consequently, only those settings in the `instances` list are significant. The others can be set to `None`. This helps to reduce the parsing time and memory footprint. It is required for large settings on more than 10.000 cores, for smaller, medium-scale problems it is not required.

When specifying the `instances` settings, it can be helpful to make use of `List Comprehensions <https://docs.python.org/3/tutorial/datastructures.html#list-comprehensions>`_, e.g. for 4 processes, like

.. code-block:: python

    "instances": [
      {
        "ranks": list(range(4)),   # use all 4 processes for all fibers
        # further settings, depending on variable i, e.g.:
        ...
            "meshName": "MeshFiber_{}".format(i)
        ...
      }
    ] for i in range(n_fibers)     # iterate over settings with i=0..n_fibers-1
