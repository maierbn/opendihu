#pragma once

#include "control/types.h"
#include "equation/type_traits.h"

/*
namespace Equation
{
template<typename Term>
using hasLaplaceOperator = std::enable_if_t<Term::hasLaplaceOperator, Term>;
}
*/

namespace SpatialDiscretization
{

/** base class of the integrand that produces the stiffness matrix
 * D = dimension (e.g. 2D, 3D), EvaluationsType = a composite data structure of multiple values, nComponents=1 for classical FE, D for solid mechanics
 */
template<int D,typename EvaluationsType,typename FunctionSpaceType,int nComponents,typename double_v_t,typename element_no_v_t,typename Term,typename=Term>
class IntegrandStiffnessMatrix
{
public:
  static EvaluationsType evaluateIntegrand(const Data::FiniteElements<FunctionSpaceType,nComponents,Term> &data,
                                           const std::array<VecD<3,double_v_t>,D> &jacobian, element_no_v_t elementNoLocal, std::array<double,D> xi);
};

/** partial specialization for laplace operator, dimension 1
 */
template<typename EvaluationsType,typename FunctionSpaceType,typename double_v_t,typename element_no_v_t,typename Term>
class IntegrandStiffnessMatrix<1,EvaluationsType,FunctionSpaceType,1,double_v_t,element_no_v_t,Term,Equation::hasLaplaceOperator<Term>>
{
public:
  static EvaluationsType evaluateIntegrand(const Data::FiniteElements<FunctionSpaceType,1,Term> &data,
                                           const std::array<VecD<3,double_v_t>,1> &jacobian, element_no_v_t elementNoLocal, const std::array<double,1> xi);
};


/** partial specialization for laplace operator, dimension 2
 */
template<typename EvaluationsType,typename FunctionSpaceType,typename double_v_t,typename element_no_v_t,typename Term>
class IntegrandStiffnessMatrix<2,EvaluationsType,FunctionSpaceType,1,double_v_t,element_no_v_t,Term,Equation::hasLaplaceOperator<Term>>
{
public:
  static EvaluationsType evaluateIntegrand(const Data::FiniteElements<FunctionSpaceType,1,Term> &data,
                                           const std::array<VecD<3,double_v_t>,2> &jacobian, element_no_v_t elementNoLocal, const std::array<double,2> xi);
};


/** partial specialization for laplace operator, dimension 3
 */
template<typename EvaluationsType,typename FunctionSpaceType,typename double_v_t,typename element_no_v_t,typename Term>
class IntegrandStiffnessMatrix<3,EvaluationsType,FunctionSpaceType,1,double_v_t,element_no_v_t,Term,Equation::hasLaplaceOperator<Term>>
{
public:
  static EvaluationsType evaluateIntegrand(const Data::FiniteElements<FunctionSpaceType,1,Term> &data,
                                           const std::array<VecD<3,double_v_t>,3> &jacobian, element_no_v_t elementNoLocal, const std::array<double,3> xi);
};


} // namespace

#include "spatial_discretization/finite_element_method/integrand/integrand_stiffness_matrix_generalized_laplace.h"
#include "spatial_discretization/finite_element_method/integrand/integrand_stiffness_matrix_laplace.tpp"
