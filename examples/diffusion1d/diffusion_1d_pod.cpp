#include <iostream>
#include <cstdlib>

#include "opendihu.h"

int main(int argc, char *argv[])
{
  // 1D diffusion equation du/dt = c du^2/dx^2
  
  // initialize everything, handle arguments and parse settings from input file
  DihuContext settings(argc, argv);
  
  PythonConfig topLevelSettings = settings.getPythonConfig();
  
  if(topLevelSettings.hasKey("ModelOrderReduction"))
  {
    DihuContext settings_timestepping=settings["ModelOrderReduction"];
    PythonConfig topLevelSettings_timeStepping = settings_timestepping.getPythonConfig();
    if(topLevelSettings_timeStepping.hasKey("ExplicitEuler"))
  {
    LOG(INFO) << "Reduced order ExplicitEuler";
    ModelOrderReduction::TimeSteppingSchemeOdeReducedExplicit<
      TimeSteppingScheme::ExplicitEuler<
        SpatialDiscretization::FiniteElementMethod<
          Mesh::StructuredRegularFixedOfDimension<1>,
          BasisFunction::LagrangeOfOrder<>,
          Quadrature::None,
          Equation::Dynamic::IsotropicDiffusion
        >
      >
    > problem(settings);
  
    problem.run();
  
    return EXIT_SUCCESS;
  } 
  else if(topLevelSettings_timeStepping.hasKey("ImplicitEuler"))
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
  else if(topLevelSettings_timeStepping.hasKey("CrankNicolson"))
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
else
  LOG(ERROR) << "No valid Model order reduction technique";
} 