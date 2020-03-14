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
