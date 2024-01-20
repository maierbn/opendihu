TimeSteppingSchemesOde
======================
Several timestepping schemes to solve ordinary differential equations are implemented. The C++ classes are listed in the following

C++ code:

.. code-block:: c
  
  TimeSteppingScheme::ExplicitEuler</* inner object, DiscretizableInTime*/>
  TimeSteppingScheme::ImplicitEuler</* inner object, DiscretizableInTime*/>
  TimeSteppingScheme::Heun</* inner object, DiscretizableInTime*/>
  TimeSteppingScheme::HeunAdaptive</* inner object, DiscretizableInTime*/>
  TimeSteppingScheme::CrankNicolson</* inner object, DiscretizableInTime*/>

They all have the following properties in common.

The following example is a configuration of the ExplicitEuler scheme:
 
.. code-block:: python

  "ExplicitEuler" : {
    "endTime":                     0.1,
    "numberTimeSteps":             5,
    "timeStepWidth":               1e-5,
    "logTimeStepWidthAsKey":       "dt_1D",
    "durationLogKey":              "duration_1D",
    "timeStepOutputInterval":      1e4,
    "initialValues":               [],
    "dirichletBoundaryConditions": {0: -75.0036, -1: -75.0036},
    "inputMeshIsGlobal":           True,
    "nAdditionalFieldVariables":   1,
    "additionalSlotNames":         ["a"],  
    "OutputWriter": ...
  }

Common properties
-------------------

In the following, the properties that can be specified for all time stepping schemes to solve ODEs (Heun, Euler, CrankNicolson) are listed.

endTime, numberTimeSteps and timeStepWidth
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
*The default values of the settings are: ``endTime``: 1, ``timestepWidth``: 0.001, ``numberTimeSteps``: 10.*

Every timestepping scheme has a  ``run()`` method and a ``advanceTimeSpan()`` method.
The ``run()`` method is called from the C++ main source file, if the timestepping scheme is the outermost class of the composite class templates.
It performs the simulation for :math:`t \in [0, \text{endTime}]`. The other method, ``advanceTimeSpan()`` is only called internally, 
if the time stepping scheme is wrapped by another solver. 
Then the wrapping class sets the values for startTime and endTime and ``advanceTimeSpan()`` simulates the time span :math:`[\text{startTime}, \text{endTime}]`.
In this case, it is not necessary to specify ``endTime`` in the python settings.

The scheme solves the equations in potentially multiple timesteps, with a time step width. The timestep width is constant within the time span, i.e. all time steps have the same width. It is determined either by ``numberTimeSteps`` or by ``timeStepWidth``. 
Only one of the two options needs to be specified. If both are given, ``numberTimeSteps`` has precedence if it is nonzero.

In the case of ``numberTimeSteps``, the time step width is computed as the time span divided by the number of time steps. 
If ``timeStepWidth`` is specified directly, this determines the time step width. 

In the case of ``timestepWidth``, the specified time step width will be modified to produce equally large time steps in the whole time span, 
such that there are no small "remainder" time steps at the end of the time span. 

The rule is that the time step width will be increased by maximum 10% or decreased otherwise.
Increasing time step width is potentially dangerous, because it can make the timestepping scheme unstable.
But ideally, changing the timestep width should not be possible at all, if it is given correctly in the settings.

Examples:

- endTime = 1.09, timeStepWidth = 1.0 -> increase dt to 1.09, one timestep
- endTime = 1.11, timeStepWidth = 1.0 -> decrease dt to 0.555, two timesteps

logTimeStepWidthAsKey, logNumberTimeStepsAsKey and durationLogKey
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
These options are optional.
If they are set, information will be stored and written to the log file  ``logs/log.csv`` at the end of the execution, under the specified key.
This log file contains one line per run of the program, i.e. at the end of execution, one line will be added. The file is in ``csv`` format, all data fields are separated by semicolons (``;``).
There will be a header that specifies the names of the logged variables (so called keys).

If ``logTimeStepWidthAsKey`` is set, the time step width of this scheme will be stored under the name, that is given by ``logTimeStepWidthAsKey``. 
If ``logNumberTimeStepsAsKey`` is set, the number of time steps of this scheme will be stored under the name, that is given by ``logNumberTimeStepsAsKey``.
This is useful to distinguish multiple runs with different time step widths.

If ``durationLogKey`` is set, there will be a duration measurement of the walltime for the timestepping scheme. The duration is measured by the ``MPI_Wtime`` call, 
and is therefore the total time that has passed for the computation. This duration will be logged under the name given by ``durationLogKey``.

If ``durationLogKey`` is not specified, the duration measurement will not take place.

timeStepOutputInterval
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
*Default: 100*

When the simulation runs, it is useful to get some information on the console, at which point it currently is. 
But, if there were output every timestep, this would slow down the simulation. In case of a nested solver structure, it may, however, be useful to output every timestep, if it takes some time to compute it.

The option ``timeStepOutputInterval`` is a positive integer value that specifies the interval in which timesteps are printed to standard output.

initialValues
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
A list of double values to use as initial values. The solution is set to these values upon initialization.

dirichletBoundaryConditions and inputMeshIsGlobal
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The Dirichlet-type boundary conditions that hold for the complete time span.
This is a dictionary with degrees of freedom as key and the value as value (i.e. ``{"dof": value, ...}``.
Negative values count from the end of possible degrees of freedom, i.e. -1 means the last degree of freedom, -2 the second last and so on.

The degrees of freedom are interpreted in global numbering, if ``inputMeshIsGlobal`` is set to ``True``, or in local numbering of the process, if ``inputMeshIsGlobal`` is ``False``.

OutputWriter
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The output writers for this time stepping scheme, see :doc:`output_writer`.

nAdditionalFieldVariables
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
(integer) A number of additional field variables that will be created. The purpose is to allow to write additional values with the output writers of this timestepping scheme.
The additional field variables can be set by connecting them via their connector slot. Their values will also be written by the output writers.

additionalSlotNames
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
A list of strings, names for of connector slots for the additional field variables. Each name should be smaller or equal than 10 characters. 
In general, named slots are used to connect the slots from a global setting "connectedSlots". See :doc:`output_connector_slots` for details.

ExplicitEuler
----------------
The explicit Euler or *forward integration* is a 1st order consistent scheme for integration of ordinary differential equations. 
The keyword for the settings is ``"ExplicitEuler"``. It only uses the common properties.

ImplicitEuler
----------------
The implicit Euler or *backward integration* is a 1st order consistent implicit scheme for integration of ordinary differential equations.
The keyword for the settings is ``"ImplicitEuler"``.
In addition to the common properties, is has the options:

.. code-block:: python
  
  "solverName" : "solver",
  "timeStepWidthRelativeTolerance" : 1e-10,
  "timeStepWidthRelativeToleranceAsKey" : "some_key",
  "durationInitTimeStepLogKey": "duration_init_1D",

``solverName`` is the name of the :doc:`solver` to use for the linear system of equations that results from the implicit scheme. 
Alternatively, the solver options can be specified directly under "ImplicitEuler", for details see the :doc:`solver` page.

The ``timeStepWidthRelativeTolerance`` is the tolerance for the time step width which controls when the system matrix has to be recomputed. 
This value will also stored in the log file if ``timeStepWidthRelativeToleranceAsKey`` is given.

If ``durationInitTimeStepLogKey`` is set, there will be a duration measurement of the walltime for the time step initialization. 
This includes both, the initial setup and potential re-initializations if the time step size changed. 

Heun
----------------
Heun integration is a 2st order consistent scheme. The keyword for the settings is ``"Heun"``.

HeunAdaptive
----------------
The HeunAdaptive class also implements the Heun method but with a time-adaptive step width. It was implemented 2019 in the Bachelor thesis by Sebastian Kreuder.

The solution is computed with different time step widths, the error is estimated and compared to a tolerance. 
If the estimated error is too high, the time step width gets dynamically reduced, down to a lower bound (``minTimeStepWidth``). 

The keyword for the settings is ``"HeunAdaptive"``.
In addition to the standard options, it has the following options:

.. code-block:: python
  
  "tolerance": 1e7,
  "lowestMultiplier": 1000,
  "minTimeStepWidth": 1e-5,
  "timeStepAdaptOption": "regular",

tolerance
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
*Default: 0.1*
The tolerance for the estimated error. It is guaranteed, that the error is always smaller than this value.

minTimeStepWidth
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
*Default: 1e-6*

The minimum timestepwidth to use. The timestep witdh will not be decreased below this value, even if the estimated error is still above the tolerance.
This avoids starving of the computation and allows to step over very badly conditioned parts of the solution process.

timeStepAdaptOption
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
*Default: regular*

Method for the adaptive time step width computation.
Possible values are ``"regular"`` and ``"modified"``. The regular method is to compute the new timestep width, :math:`dt_\text{new}`, by

.. math::

  \alpha = \left(\dfrac{\text{tolerance}}{\text{estimator}}\right)^{1/3}
  
  dt_\text{new} = 0.9 \cdot \alpha \cdot dt_\text{old}

Care is taken to not produce too small remainder timesteps (< 0.1*dt) at the end of the time span.

The modified version only allows equally sized time step widths for the total time span. The time step width is fixed at the beginning of the time span. This is only usefull in the inner timestepping of a splitting scheme.
For details on the effects see the Bachelor thesis document.

lowestMultiplier
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
*Default: 1000*
This is the minimum number of timesteps to perform in the time span for the "modified" method. E.g. by default there will be at least 1000 time steps in the time span.


CrankNicolson
-------------------
The Crank Nicolson scheme is implicit and 2nd order consistent. 
The keyword for the settings is ``"CrankNicolson"``.
In addition to the common properties, it has the options:

.. code-block:: python
  
  "solverName" : "solver"
  "timeStepWidthRelativeTolerance" : 1e-10,
  "timeStepWidthRelativeToleranceAsKey" : "some_key",
  "durationInitTimeStepLogKey": "duration_init_1D",

``solverName`` is the name of the :doc:`solver` to use for the linear system of equations that results from the implicit scheme. 
Alternatively, the solver options can be specified directly under "CrankNicolson", for details see the :doc:`solver` page. 

The ``timeStepWidthRelativeTolerance`` is the tolerance for the time step width which controls when the system matrix has to be recomputed. 
This value will also stored in the log file if ``timeStepWidthRelativeToleranceAsKey`` is given.

If ``durationInitTimeStepLogKey`` is set, there will be a duration measurement of the walltime for the time step initialization. 
This includes both, the initial setup and potential re-initializations if the time step size changed. 
