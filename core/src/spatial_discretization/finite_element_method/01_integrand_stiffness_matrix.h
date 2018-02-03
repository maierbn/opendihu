#pragma once

#include "equation/type_traits.h"
#include "control/types.h"

namespace SpatialDiscretization
{

/** base class of the integrand that produces the stiffness matrix 
 */
template<int D,typename EvaluationsType,typename BasisOnMeshType,typename Term,typename=Term>
class IntegrandStiffnessMatrix
{
public:
  static EvaluationsType evaluateIntegrand(const std::shared_ptr<BasisOnMeshType> mesh, 
                                           const std::array<Vec3,D> &jacobian,std::array<double,D> xi);
};
 
/** partial specialization for laplace operator, dimension 1
 */
template<typename EvaluationsType,typename BasisOnMeshType,typename Term>
class IntegrandStiffnessMatrix<1,EvaluationsType,BasisOnMeshType,Term,Equation::hasLaplaceOperator<Term>>
{
public:
  static EvaluationsType evaluateIntegrand(const std::shared_ptr<BasisOnMeshType> mesh,
                                           const std::array<Vec3,1> &jacobian, const std::array<double,1> xi);
};

 
/** partial specialization for laplace operator, dimension 2
 */
template<typename EvaluationsType,typename BasisOnMeshType,typename Term>
class IntegrandStiffnessMatrix<2,EvaluationsType,BasisOnMeshType,Term,Equation::hasLaplaceOperator<Term>>
{
public:
  static EvaluationsType evaluateIntegrand(const std::shared_ptr<BasisOnMeshType> mesh,
                                           const std::array<Vec3,2> &jacobian, const std::array<double,2> xi);
};

 
/** partial specialization for laplace operator, dimension 3
 */
template<typename EvaluationsType,typename BasisOnMeshType,typename Term>
class IntegrandStiffnessMatrix<3,EvaluationsType,BasisOnMeshType,Term,Equation::hasLaplaceOperator<Term>>
{
public:
  static EvaluationsType evaluateIntegrand(const std::shared_ptr<BasisOnMeshType> mesh,
                                           const std::array<Vec3,3> &jacobian, const std::array<double,3> xi);
};

};  // namespace

#include "spatial_discretization/finite_element_method/01_integrand_stiffness_matrix.tpp"
