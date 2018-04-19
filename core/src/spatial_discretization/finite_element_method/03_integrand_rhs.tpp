#include "spatial_discretization/finite_element_method/03_integrand_rhs.h"

#include <array>

#include "utility/math_utility.h"

namespace SpatialDiscretization
{

template<int D,typename EvaluationsType,typename BasisOnMeshType,typename Term>
EvaluationsType IntegrandRightHandSide<D,EvaluationsType,BasisOnMeshType,Term,Equation::hasRhs<Term>>::
evaluateIntegrand(const std::array<Vec3,D> &jacobian, const std::array<double,D> xi)
{
  EvaluationsType evaluations;
  
  // get the factor in the integral that arises from the change in integration domain from world to coordinate space
  double integrationFactor = MathUtility::computeIntegrationFactor<D>(jacobian);

  // loop over pairs of basis functions and evaluation integrand at xi
  for (int i=0; i<BasisOnMeshType::nDofsPerElement(); i++)
  {
    for (int j=0; j<BasisOnMeshType::nDofsPerElement(); j++)
    {
      //VLOG(2) << "    integrationFactor " << integrationFactor << ", jacobian: " << jacobian;
      double integrand = BasisOnMeshType::phi(i,xi) * BasisOnMeshType::phi(j,xi) * integrationFactor;
      evaluations(i,j) = integrand;
    }
  }
  
  return evaluations;
};

};  // namespace