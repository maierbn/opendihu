#pragma once

#include "basis_function/basis_function.h"
#include "basis_function/lagrange.h"
#include "control/types.h"

namespace BasisFunction
{

class Hermite : public BasisFunction
{
public:

  typedef LagrangeOfOrder<1> BasisFunctionUsingOnlyNodalValues;  //< the same basis function, but only using nodal values, for Hermite this is the linear Lagrange functions. This construct is used for Neumann BC, where the integration over the surface/line is done using Lagrange ansatz functions.

  //! number of degrees of freedom of this basis
  static constexpr int nDofsPerBasis();

  //! number of degrees of freedom associated with a node in world space
  static constexpr int nDofsPerNode();

  //! evaluate the 1D basis function corresponding to element-local dof i at xi, interval for xi is [0,1]
  static double phi(int i, double xi);

  //! evaluate the first derivative of the 1D basis function corresponding to element-local dof i at xi, interval for xi is [0,1]
  static double dphi_dxi(int i, double xi);

  //! return the basis order value as used in python files and callbacks, i.e. 3
  static int getBasisOrder();

  //! return a basis function type string as used in python files and callbacks, i.e. "Hermite"
  static std::string getBasisFunctionString();

  static constexpr bool isNodalBased = true;  //< specify that this basis function is nodal based
};

}  // namespace

#include "basis_function/hermite.tpp"
