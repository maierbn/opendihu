#pragma once

#include "equation/type_traits.h"
#include "control/types.h"

namespace SpatialDiscretization
{

/** class of integrand for rhs
 */
template<int D,typename EvaluationsType,typename FunctionSpaceType,int nComponents,typename double_v_t,typename Term,typename=Term>
class IntegrandMassMatrix
{
public:
  static EvaluationsType evaluateIntegrand(const std::array<VecD<3,double_v_t>,D> &jacobian, const std::array<double,D> xi,
                                           std::shared_ptr<FunctionSpaceType> functionSpace, element_no_t elementNoLocal);
};



} // namespace

#include "spatial_discretization/finite_element_method/integrand/integrand_mass_matrix.tpp"
