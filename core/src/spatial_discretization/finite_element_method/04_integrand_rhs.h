#pragma once

#include "equation/type_traits.h"
#include "control/types.h"

namespace SpatialDiscretization
{

/** class of integrand for rhs
 */
template<int D,typename EvaluationsType,typename BasisOnMeshType,typename Term,typename=Term>
class IntegrandRightHandSide
{
public:
};

/** partial specialization for equations that have a rhs
 */
template<int D,typename EvaluationsType,typename BasisOnMeshType,typename Term>
class IntegrandRightHandSide<D,EvaluationsType,BasisOnMeshType,Term,Equation::hasRhs<Term>>
{
public:
  static EvaluationsType evaluateIntegrand(const std::array<Vec3,D> &jacobian, const std::array<double,D> xi);
};



};  // namespace

#include "spatial_discretization/finite_element_method/04_integrand_rhs.tpp"
