#pragma once
#include "basis_function/basis_function.h"
#include "control/types.h"

namespace BasisFunction
{

/** Complete polynomials of a dimension and order contain all terms up to the given polynomial degree.
 *  The function space is called P_k for order (polynomial degree) k.
 *  2D Examples: order=0: p=c0, order=1: p=c0 + c1*x + c2*y, order=2: p=c0 + c1*x^2 + c2*x*y + c3*y^2
 *  Coefficients of the basis functions in FE computation are later element-based, not node-based like Lagrange basis functions.
 *  The class CompletePolynomialNDofs computes the numbers of dofs for D and order.
  */
template<int D, int order>
class CompletePolynomialNDofs
{};

/** class that computes number of dofs for 1D
 */
template<int order>
class CompletePolynomialNDofs<1,order>
{
public:
  //! number of degrees of freedom of this basis
  static constexpr int nDofsPerBasis();
};

/** class that computes number of dofs for 2D
 */
template<int order>
class CompletePolynomialNDofs<2,order>
{
public:
  //! number of degrees of freedom of this basis
  static constexpr int nDofsPerBasis();
};

/** class that computes number of dofs for 3D
 */
template<int order>
class CompletePolynomialNDofs<3,order>
{
public:
  //! number of degrees of freedom of this basis
  static constexpr int nDofsPerBasis();
};
 
/**  Defines complete polynomials of the order, space P_k in  interval xi in [0,1] of given order.
 *   2D Examples: order=0: p=c0, order=1: p=c0 + c1*x + c2*y, order=2: p=c0 + c1*x^2 + c2*x*y + c3*y^2
  */
template<int D, int order>
class CompletePolynomialOfDimensionAndOrder : 
  public BasisFunction,
  public CompletePolynomialNDofs<D,order>
{
public:
  
  //! number of degrees of freedom associated with a node in world space. This is 0 because the degrees of freedom are associated with elements.
  static constexpr int nDofsPerNode();
  
  //! evaluate the basis function corresponding to element-local dof dofIndex at xi, xi lives in [0,1]^D
  static double phi(int dofIndex, std::array<double,D> xi);
  
  //! evaluate the derivative of phi(xi) w.r.t xi_i, where i is given by derivativeIdx, i.e. Phi_{dofIndex,derivativeIdx}(xi)
  static double dphi_dxi(int dofIndex, int derivativeIdx, std::array<double,D> xi);
  
  //! return the basis order value as used in python files and callbacks, e.g. 2
  static constexpr int getBasisOrder();
  
  //! return a basis function type string as used in python files and callbacks, e.g. "Lagrange"
  static std::string getBasisFunctionString();
private:
};



}  // namespace

#include "basis_function/complete_polynomial.tpp"
