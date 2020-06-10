#include <iostream>
#include <cstdlib>

#include "opendihu.h"

int main(int argc, char *argv[])
{
  // 3D multidomain: implicit Euler, FEM
  
  // initialize everything, handle arguments and parse settings from input file
  DihuContext settings(argc, argv);
  
  typedef Mesh::StructuredDeformableOfDimension<3> MeshType;

  MultidomainWrapper<
    OperatorSplitting::Strang<
      Control::MultipleInstances<
        TimeSteppingScheme::Heun<
          CellmlAdapter<
            4,9,  // nStates,nIntermediates: 57,1 = Shorten, 4,9 = Hodgkin Huxley
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
  > problem(settings);
  
  problem.run();
  
  
  // test methods of wrapper
  // get the number of local solution values
  int nSolutionValuesLocal = problem.nSolutionValuesLocal();
  
  // allocate a vector to hold all solution values
  std::vector<double> values(nSolutionValuesLocal);
  
  // get the current solution values
  problem.getSolution(values.data(), 0, 0);
  LOG(INFO) << "number of values: " << nSolutionValuesLocal << ", values: " << values;
  
  // set the current solution values
  problem.setSolution(values.data());
  
  return EXIT_SUCCESS;
}
