#include <iostream>
#include <cstdlib>

#include "opendihu.h"

int main(int argc, char *argv[])
{
  // 1D diffusion equation du/dt = c du^2/dx^2
  
  // initialize everything, handle arguments and parse settings from input file
  DihuContext settings(argc, argv);
  
  PyObject *topLevelSettings = settings.getPythonConfig();
  
  if(PythonUtility::hasKey(topLevelSettings, "ModelOrderReduction"))
  {
    DihuContext settings_timestepping=settings["ModelOrderReduction"];
    PyObject *topLevelSettings_timeStepping = settings_timestepping.getPythonConfig();
  
    if(PythonUtility::hasKey(topLevelSettings_timeStepping, "ExplicitEuler"))
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
  else if(PythonUtility::hasKey(topLevelSettings, "ImplicitEuler"))
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
  else if(PythonUtility::hasKey(topLevelSettings, "CrankNicolson"))
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