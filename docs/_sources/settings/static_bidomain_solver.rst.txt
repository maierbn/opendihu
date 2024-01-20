StaticBidomainSolver
=====================

This class solves the static bidomain equation. It can be used to compute the electromyography (EMG) signal in a muscle from :math:`V_m` values on the fibers.

.. math::

  \mathrm{div}\big((\sigma_i+\sigma_e)\,\mathrm{grad}(\phi_e)\big) + \mathrm{div}\big(\sigma_i\,\mathrm{grad}(V_m)\big) = 0

The value to be solved for is the extra-cellular potential, :math:`\phi_e` which corresponds to the EMG values. The input is the transmembrane potential, :math:`V_m`. It should be provided to the `StaticBidomainSolver` by the connector slot ``Vm``.

C++ instantiation
-----------------

The `StaticBidomainSolver` class needs two `FiniteElementMethod` classes, the first for a potential flow to estimate fiber directions (for the conductivity tensors), the second one for the actual gradient operator of the equation.

.. code-block:: c

  TimeSteppingScheme::StaticBidomainSolver<
    SpatialDiscretization::FiniteElementMethod<       //FEM for initial potential flow, fiber directions
      Mesh::StructuredDeformableOfDimension<3>,
      BasisFunction::LagrangeOfOrder<1>,
      Quadrature::Gauss<3>,
      Equation::Static::Laplace
    >,
    SpatialDiscretization::FiniteElementMethod<       // anisotropic diffusion
      Mesh::StructuredDeformableOfDimension<3>,
      BasisFunction::LagrangeOfOrder<1>,
      Quadrature::Gauss<5>,
      Equation::Dynamic::DirectionalDiffusion
    >
  >

Python Settings
^^^^^^^^^^^^^^^^^^^

The python settings are as given below. The values for ``"PotentialFlow"`` and ``"Activation"`` are the same as for a :doc:`finite_element_method`.

.. code-block:: python

  "StaticBidomainSolver": {
    "timeStepWidth":          variables.dt_3D,
    "timeStepOutputInterval": 50,
    "durationLogKey":         "duration_bidomain",
    "solverName":             "activationSolver",
    "initialGuessNonzero":    variables.emg_initial_guess_nonzero,
    "slotNames:"              [],
    "PotentialFlow": {
      "FiniteElementMethod" : {
        "meshName":           "3Dmesh",
        "solverName":         "potentialFlowSolver",
        "prefactor":          1.0,
        "dirichletBoundaryConditions": variables.potential_flow_dirichlet_bc,
        "neumannBoundaryConditions":   [],
        "inputMeshIsGlobal":  True,
      },
    },
    "Activation": {
      "FiniteElementMethod" : {
        "meshName":           "3Dmesh",
        "solverName":         "activationSolver",
        "prefactor":          1.0,
        "inputMeshIsGlobal":  True,
        "dirichletBoundaryConditions": {},
        "neumannBoundaryConditions":   [],
        "diffusionTensor": [[      # sigma_i, fiber direction is (1,0,0), one list item = same tensor for all elements, multiple list items = a different tensor for each element
          8.93, 0, 0,
          0, 0.893, 0,
          0, 0, 0.893
        ]],
        "extracellularDiffusionTensor": [[      # sigma_e, one list item = same tensor for all elements, multiple list items = a different tensor for each element
          6.7, 0, 0,
          0, 6.7, 0,
          0, 0, 6.7,
        ]],
      },
    },
    "OutputWriter" : variables.output_writer_emg,   #  list of output writers
  }
  
slotNames
----------
A list of strings, names for the connector slots. Each name should be smaller or equal than 10 characters. 
In general, named slots are used to connect the slots from a global setting "connectedSlots". See :doc:`output_connector_slots` for details.
