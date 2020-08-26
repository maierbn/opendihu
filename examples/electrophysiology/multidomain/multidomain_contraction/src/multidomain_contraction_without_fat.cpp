#include <iostream>
#include <cstdlib>

#include "opendihu.h"

int main(int argc, char *argv[])
{
  // 3D multidomain coupled with contraction
  
  // initialize everything, handle arguments and parse settings from input file
  DihuContext settings(argc, argv);
  
  typedef Mesh::StructuredDeformableOfDimension<3> MeshType;

  Control::Coupling
  <
    OperatorSplitting::Strang<
      Control::MultipleInstances<
        TimeSteppingScheme::Heun<
          CellmlAdapter<
            57,71,  // nStates,nAlgebraics: 57,1 = Shorten, 4,9 = Hodgkin Huxley
            FunctionSpace::FunctionSpace<MeshType,BasisFunction::LagrangeOfOrder<1>>  // same function space as for anisotropic diffusion
          >  
        >
      >,
      TimeSteppingScheme::MultidomainWithFatSolver<       // multidomain
        SpatialDiscretization::FiniteElementMethod<       //FEM for initial potential flow, fibre directions
          MeshType,
          BasisFunction::LagrangeOfOrder<1>,
          Quadrature::Gauss<3>,
          Equation::Static::Laplace
        >,
        SpatialDiscretization::FiniteElementMethod<       // anisotropic diffusion
          MeshType,
          BasisFunction::LagrangeOfOrder<1>,
          Quadrature::Gauss<5>,
          Equation::Dynamic::DirectionalDiffusion
        >,
        SpatialDiscretization::FiniteElementMethod<       // isotropic diffusion in fat layer
          MeshType,
          BasisFunction::LagrangeOfOrder<1>,
          Quadrature::Gauss<5>,
          Equation::Dynamic::IsotropicDiffusion
        >
      >
    >,
    MuscleContractionSolver<>
  > problem(settings);
  
  problem.run();
  
  return EXIT_SUCCESS;
}
