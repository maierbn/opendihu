#pragma once

#include "equation/type_traits.h"
#include "control/types.h"

namespace SpatialDiscretization
{

/** class of integrand for rhs
 */
template<int D,typename EvaluationsType,typename FunctionSpaceType,typename Term,typename=Term>
class IntegrandMassMatrix
{
public:
  static EvaluationsType evaluateIntegrand(const std::array<Vec3,D> &jacobian, const std::array<double,D> xi);
};



};  // namespace

#include "spatial_discretization/finite_element_method/integrand/integrand_mass_matrix.tpp"
