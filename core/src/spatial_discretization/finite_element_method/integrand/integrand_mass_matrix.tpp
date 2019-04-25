#include "spatial_discretization/finite_element_method/integrand/integrand_mass_matrix.h"

#include <array>

#include "utility/math_utility.h"

namespace SpatialDiscretization
{

template<int D,typename EvaluationsType,typename FunctionSpaceType,typename Term,typename Dummy>
EvaluationsType IntegrandMassMatrix<D,EvaluationsType,FunctionSpaceType,Term,Dummy>::
evaluateIntegrand(const std::array<Vec3,D> &jacobian, const std::array<double,D> xi)
{
  EvaluationsType evaluations;

  // get the factor in the integral that arises from the change in integration domain from world to coordinate space
  double integrationFactor = MathUtility::computeIntegrationFactor<D>(jacobian);

  // loop over pairs of basis functions and evaluation integrand at xi
  for (int i=0; i<FunctionSpaceType::nDofsPerElement(); i++)
  {
    for (int j=0; j<FunctionSpaceType::nDofsPerElement(); j++)
    {
      //VLOG(2) << "    integrationFactor " << integrationFactor << ", jacobian: " << jacobian;
      double integrand = FunctionSpaceType::phi(i,xi) * FunctionSpaceType::phi(j,xi) * integrationFactor;
      evaluations(i,j) = integrand;
    }
  }

  return evaluations;
}

} // namespace
