MultidomainSolver
===================

This solves the multidomain equation.

C++ code of usage:

.. code-block:: c
  
  OperatorSplitting::Strang<
    Control::MultipleInstances<
      TimeSteppingScheme::Heun<
        CellmlAdapter<
          4,   // 57 for Hodgkin-Huxley
          FunctionSpace::FunctionSpace<MeshType,BasisFunction::LagrangeOfOrder<1>>  // same function space as for anisotropic diffusion
        >  
      >
    >,
    TimeSteppingScheme::MultidomainSolver<              // multidomain
      SpatialDiscretization::FiniteElementMethod<       //FEM for initial potential flow, fibre directions
        MeshType,
        BasisFunction::LagrangeOfOrder<1>,
        Quadrature::Gauss<3>,
        Equation::Static::Laplace
      >,
      SpatialDiscretization::FiniteElementMethod<   // anisotropic diffusion
        MeshType,
        BasisFunction::LagrangeOfOrder<1>,
        Quadrature::Gauss<5>,
        Equation::Dynamic::DirectionalDiffusion
      >
    >
  >
  
See the multidomain3d example for details.
