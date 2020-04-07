#include "spatial_discretization/finite_element_method/integrand/integrand_mass_matrix.h"

#include <array>

#include "utility/math_utility.h"

namespace SpatialDiscretization
{

template<int D,typename EvaluationsType,typename FunctionSpaceType,int nComponents,typename double_v_t,typename Term,typename Dummy>
EvaluationsType IntegrandMassMatrix<D,EvaluationsType,FunctionSpaceType,nComponents,double_v_t,Term,Dummy>::
evaluateIntegrand(const std::array<VecD<3,double_v_t>,D> &jacobian, const std::array<double,D> xi)
{
  EvaluationsType evaluations;

  // get the factor in the integral that arises from the change in integration domain from world to coordinate space
  double_v_t integrationFactor = MathUtility::computeIntegrationFactor(jacobian);

  // loop over pairs of basis functions and evaluation integrand at xi
  for (int i = 0; i < FunctionSpaceType::nDofsPerElement(); i++)
  {
    for (int j = 0; j < FunctionSpaceType::nDofsPerElement(); j++)
    {
      // loop over components, for not solid mechanics, nComponents is 1
      for (int componentNo = 0; componentNo < nComponents; componentNo++)
      {
        double_v_t integrand = FunctionSpaceType::phi(i,xi) * FunctionSpaceType::phi(j,xi) * integrationFactor;
        LOG(DEBUG) << "    integrationFactor " << integrationFactor << ", jacobian: " << jacobian << ", integrationFactor: " << integrationFactor << " -> " << integrand;
        evaluations(i*nComponents + componentNo, j*nComponents + componentNo) = integrand;
      }
    }
  }

  return evaluations;
}

} // namespace
