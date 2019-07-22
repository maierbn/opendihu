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
