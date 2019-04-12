#include <iostream>
#include <cstdlib>

#include "opendihu.h"

int main(int argc, char *argv[])
{
  // 1D diffusion equation du/dt = c du^2/dx^2
  
  // initialize everything, handle arguments and parse settings from input file
  DihuContext settings(argc, argv);
  
  PythonConfig topLevelSettings = settings.getPythonConfig(); 
  
  if(topLevelSettings.hasKey("ExplicitEuler"))
  {
    LOG(INFO) << "ExplicitEuler";
    
    TimeSteppingScheme::ExplicitEuler<
    SpatialDiscretization::FiniteElementMethod<
      Mesh::StructuredRegularFixedOfDimension<1>,
      BasisFunction::LagrangeOfOrder<>,
      Quadrature::None,
      Equation::Dynamic::IsotropicDiffusion
    >
    > problem(settings);
  
    problem.run();
  
    return EXIT_SUCCESS;
  } 
  else if(topLevelSettings.hasKey("ImplicitEuler"))
  {
    LOG(INFO) << "ImplicitEuler";
    
    TimeSteppingScheme::ImplicitEuler<
    SpatialDiscretization::FiniteElementMethod<
      Mesh::StructuredRegularFixedOfDimension<1>,
      BasisFunction::LagrangeOfOrder<>,
      Quadrature::None,
      Equation::Dynamic::IsotropicDiffusion
    >
    > problem(settings);
    
    problem.run();
    
    return EXIT_SUCCESS;
  }
  else if(topLevelSettings.hasKey("CrankNicolson"))
  {
    LOG(INFO) << "CrankNicolson";
    
    TimeSteppingScheme::CrankNicolson<
    SpatialDiscretization::FiniteElementMethod<
    Mesh::StructuredRegularFixedOfDimension<1>,
    BasisFunction::LagrangeOfOrder<>,
    Quadrature::None,
    Equation::Dynamic::IsotropicDiffusion
    >
    > problem(settings);
    
    problem.run();
    
    return EXIT_SUCCESS;
  }
  else
    LOG(ERROR) << "No valid time integration scheme in settings.py";
   
}
