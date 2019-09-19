#pragma once

#include "equation/type_traits.h"
#include "control/types.h"

namespace SpatialDiscretization
{

/** partial specialization for finite elasticity
 */
template<int D, typename EvaluationsType,typename FunctionSpaceType,typename Term>
class IntegrandStiffnessMatrix<D,EvaluationsType,FunctionSpaceType,D,Term,Equation::isLinearElasticity<Term>>
{
public:
  static EvaluationsType evaluateIntegrand(const Data::FiniteElements<FunctionSpaceType,D,Term> &data,
                                           const std::array<Vec3,D> &jacobian, element_no_t elementNoLocal,
                                           const std::array<double,D> xi);
protected:
  //! evaluate the stiffness tensor C_abcd
  static double stiffness(const Data::FiniteElements<FunctionSpaceType,D,Term> &data, int a, int b, int c, int d);
};

} // namespace

#include "spatial_discretization/finite_element_method/integrand/integrand_stiffness_matrix_linear_elasticity.tpp"
