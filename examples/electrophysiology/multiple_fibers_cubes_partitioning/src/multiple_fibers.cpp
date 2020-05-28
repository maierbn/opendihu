#include <Python.h>
#include <iostream>
#include <cstdlib>

#include <iostream>
#include "easylogging++.h"

#include "opendihu.h"

int main(int argc, char *argv[])
{
  // 1D reaction-diffusion equation du/dt = c du^2/dx^2 + R(t), R is from cellml file
  
  // initialize everything, handle arguments and parse settings from input file
  DihuContext settings(argc, argv);
  
  // define problem
  Control::MultipleInstances<
    OperatorSplitting::Strang<
      Control::MultipleInstances<
        TimeSteppingScheme::Heun<
          CellmlAdapter<
            4,9,  // nStates,nAlgebraics: 57,1 = Shorten, 4,9 = Hodgkin Huxley, 57,71 = slow_TK_2014_12_08
            FunctionSpace::FunctionSpace<
              Mesh::StructuredDeformableOfDimension<1>,
              BasisFunction::LagrangeOfOrder<1>
            >
          >  
        >
      >,
      Control::MultipleInstances<
        TimeSteppingScheme::CrankNicolson<
          SpatialDiscretization::FiniteElementMethod<
            Mesh::StructuredDeformableOfDimension<1>,
            BasisFunction::LagrangeOfOrder<1>,
            Quadrature::Gauss<2>,
            Equation::Dynamic::IsotropicDiffusion
          >
        >
      >
    >
  > problem(settings);
  
  // run problem
  problem.run();
  
  return EXIT_SUCCESS;
}




