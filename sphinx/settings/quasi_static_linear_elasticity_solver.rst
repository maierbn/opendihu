QuasiStaticLinearElasticitySolver
===================================

This solves quasi static linear elasticity.

C++ code:

.. code-block:: c

  TimeSteppingScheme::QuasiStaticLinearElasticitySolver<
    SpatialDiscretization::FiniteElementMethod<       // finite element method for static linear elasticity, has the normal options:
      Mesh::StructuredDeformableOfDimension<3>,       // Mesh
      BasisFunction::LagrangeOfOrder<1>,              // BasisFunction
      Quadrature::Gauss<3>,                           // Quadrature
      Equation::Static::LinearElasticityActiveStress
    >
  >

dirichletBoundaryConditions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The value to set is a vector because displacements in all spatial directions can be prescribed. 
Specify a list of the components for each prescribed dof, e.g. ``[1.0, 2.0, 3.0]`` to set a Dirichlet boundary condition of :math:`\bar{u} = (1,2,3)^\top`. When not all components should be prescribed, replace the entry by ``None``, e.g. ``[None, 2.0, None]`` to only prescribe the y component.

Python Settings
^^^^^^^^^^^^^^^^^^^

The python settings are as given below:

.. code-block:: python

  # timestepping solver that solves quasi-static elasticity formulation
    "QuasiStaticLinearElasticitySolver": {
      "fiberDirection":           [0, 0, 1],      # direction for anisotropy of elasticity formulation
      "PotentialFlow":            None,           # if fiberDirection is not set to a constant direction, a potential flow simulation can be used where the fiber direction is set to the streamlines of the flow through the volume. In this case, set "PotentialFlow" to the settings for the FEM for the potential flow.
      "maximumActiveStress":      1.0,            # scaling factor to the active stress, σ_active = activation * anisotropyTensor * maximumActiveStress
      "strainScalingCurveWidth":  1.0,            # parameter for strain-stress curve of active stress, has no effect, because strain-stress curve is commented out in the code
      "scalingFactor":            1.0,            # scaling factor for displacements, to overrate them, if != 0 it is only for visualization purposes and not physical
      "inputMeshIsGlobal":        True,           # if boundary conditions are specified in global numbering
      "slotNames":                [], 
        
      # Anisotropy for active stress.
      # The tensor is given in a local basis where the fiber direction is (1,0,0), one list item = same tensor for all elements, multiple list items = a different tensor for each element.
      # The tensor has to be symmetric.
      "anisotropyTensor": [              
        1, 0, 0,
        0, 0, 0,
        0, 0, 0
      ],
      
      # linear elasticity finite element method
      "FiniteElementMethod" : {   
        "meshName":             "elasticityMesh",       # reference to the mesh, defined under "Meshes"
        "solverName":           "linearElasticitySolver",   # reference to the linear system solver, defined under "Solvers"
        "prefactor":            1.0,                    # prefactor c of the equation c*Δu = f
        "inputMeshIsGlobal":    True,                   # if boundary conditions are specified in global numbering
        "dirichletBoundaryConditions":  dirichlet_bc,   # dirichlet boundary conditions
        "neumannBoundaryConditions":    neumann_bc,     # neumann boundary conditions
        "divideNeumannBoundaryConditionValuesByTotalArea": False,  # if the neumann boundary condition vectors should be divided by the total surface area where surface loads are applied, this allows to specify the total force that acts on the surface. If set to False (default), the given traction is a per-surface quantity.
              
        # material parameters
        "bulkModulus":        150,                      # bulk modulus, K, material parameter for compressibility
        "shearModulus":       20,                       # shear modulus, μ
      },
      
      "OutputWriter": [
        {"format": "Paraview", "outputInterval": 1, "filename": "out/elasticity", "binary": False, "fixedFormat": False, "onlyNodalValues":True, "combineFiles":True, "fileNumbering": "incremental"},
        #{"format": "PythonFile", "filename": "out/elasticity", "outputInterval": 1, "binary":False, "onlyNodalValues":True, "fileNumbering": "incremental"},
      ]
    }
  }

slotNames
----------
A list of strings, names for the connector slots. Each name should be smaller or equal than 10 characters. 
In general, named slots are used to connect the slots from a global setting "connectedSlots". See :doc:`output_connector_slots` for details.
