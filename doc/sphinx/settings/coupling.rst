Coupling
=========

The coupling class takes two timestepping schemes and executes them alternatingly. It sets divides the total time span in smaller time spans. For each of the
smaller time spans :math:`[t_i,t_{i+1}]`, it does the following:

- set time span :math:`[t_i,t_{i+1}]` in time stepping scheme 1
- solve time stepping scheme 1 (call ``advanceTimeSpan()``)
- transfer data from time stepping scheme 1 to 2
- set time span :math:`[t_i,t_{i+1}]` in time stepping scheme 2
- solve time stepping scheme 2 (call ``advanceTimeSpan()``)
- transfer data from time stepping scheme 2 to 1

The Coupling is itself a timestepping scheme and has all options listed on the (:doc:`timestepping_schemes_ode`) page.

It is equivalent to :doc:`splitting`, but with the keyword "Coupling" instead of "GodunovSplitting".

Note that usually data should be transferred between the two solvers. This is done with the settings ``connectedSlotsTerm1To2`` 
and ``connectedSlotsTerm2To1``. See the notes on :doc:`/settings/output_connector_slots`.

C++ code:

.. code-block:: c

  Control::Coupling<
    /*first timestepping scheme*/,
    /*second timestepping scheme*/
  >

Python settings:

.. code-block:: python

  "Coupling": {
    "numberTimeSteps": 1,
    "timeStepWidth": 1e-1,
    "endTime": 10.0,
    "logTimeStepWidthAsKey": "dt_3D",
    "durationLogKey": "duration_total",
    "timeStepOutputInterval": 10,
    
    "connectedSlotsTerm1To2": [0],  # list of slots of term 2 that are connected to the slots of term 1
    "connectedSlotsTerm2To1": [0],  # list of slots of term 1 that are connected to the slots of term 2
    
    "Term1": {
      "ExplicitEuler" : {  # these are the properties for the first timestepping
      }
    },
    "Term2": {
      "ImplicitEuler" : {  # these are the properties for the second timestepping
      }
    }
  }

MultipleCoupling
=================

If you need more than two child solvers that should be executed in a row, you can use `MultipleCoupling`. 

Without this, a nested structure of `Coupling` schemes would be needed. For example, for 3 solvers `A`, `B` and `C`, you would do something like:

.. code-block:: c

  Control::Coupling<
    A,
    Control::Coupling<
      B,
      C
    >
  >

with settings 

.. code-block:: python

  "Coupling": {
    "Term1": {
      # settings for A
    },
    "Term2": {
      "Coupling": {
        "Term1": {
          # settings for B
        },
        "Term2": {
          # settings for C
        }
      }
    }
  }

The `MultipleCoupling` abstracts this construction for any number of subsolvers. The usage is

.. code-block:: c

  Control::MultipleCoupling<
    A,B,C, //...   any number of subsolvers (at least 2)
  >

This is more convenient than the construction with the nested Coupling schemes.
Also the Python settings are simplified:

.. code-block:: python

  "MultipleCoupling": {
    "endTime":                variables.end_time,           # end time of the simulation
    "timeStepWidth":          variables.dt_neuron_transfer, # time step width of the data transfer between the sub solvers
    "logTimeStepWidthAsKey":  "dt_neuron_transfer",         # string under which the timestep width will be stored in the log file
    "durationLogKey":         "duration_total",             # string under which the total duration will be stored in the log file
    "timeStepOutputInterval": 500,                          # how often to print the current time step
    "connectedSlotsTerm1To2": None,                         # this would be the connected slots for a normal Coupling scheme, but here it should be set to None, use global option "connectedSlots" instead
    "connectedSlotsTerm2To1": None,                         # this would be the connected slots for a normal Coupling scheme, but here it should be set to None, use global option "connectedSlots" instead

    "Term1": {
      # settings for A
    },
    "Term2": {
      # settings for B
    },
    "Term3": {
      # settings for C
    },
    // etc.
  }

Internally, `MultipleCoupling` constructs the tree of nested `Coupling` schemes. It is implemented using variadic templates and therefore works for any number of nested solvers. 

The settings for `MultipleCoupling` are shared for all constructed `Coupling` schemes. This means that the set timestep width is reused for all terms of the coupling. The options ``"connectedSlotsTerm1To2"`` and ``"connectedSlotsTerm2To1"`` should be set to ``None``. The slots should be connected using the `"connectedSlots"` option, for details see :doc:`output_connector_slots`.
The settings for the nested solvers can be given under `"Term1"`, `"Term2"`, `"Term3"`, etc., as with the normal Coupling scheme.


