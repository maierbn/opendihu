#include <iostream>
#include <cstdlib>

#include "opendihu.h"

int main(int argc, char *argv[])
{
  // linear elasticity
  
  // initialize everything, handle arguments and parse settings from input file
  DihuContext settings(argc, argv);
  
  // coupling scheme that repeatedly calls the two nested solvers
  Control::Coupling<
  
    // class that prescribes values on a 3D mesh
    PrescribedValues<
      FunctionSpace::FunctionSpace<
        Mesh::StructuredDeformableOfDimension<3>,
        BasisFunction::LagrangeOfOrder<1>
      >
    >,
    
    // timestepping scheme that wraps the linear elasticity and computes the quasi-static problem, this also add the active stress term to the right hand side
    TimeSteppingScheme::QuasiStaticLinearElasticitySolver<              
    
      // finite element method that solves Δu = f, in this case with u ∈ ℝ^3
      SpatialDiscretization::FiniteElementMethod< 
        Mesh::StructuredDeformableOfDimension<3>,
        BasisFunction::LagrangeOfOrder<1>,
        Quadrature::Gauss<3>,
        Equation::Static::LinearElasticityActiveStress
      >
    >
  > problem(settings);
  
  problem.run();
  
  return EXIT_SUCCESS;
}
