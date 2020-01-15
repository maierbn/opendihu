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
