#pragma once

#include "basis_function/basis_function.h"
#include "control/types.h"

namespace BasisFunction
{

template<int order=1>
class Lagrange : public BasisFunction
{
public:
  
  //! number of degrees of freedom of this basis
  static constexpr int nDofsPerBasis();
  
  //! number of degrees of freedom associated with a node in world space
  static constexpr int nDofsPerNode();
  
  //! average number of dofs per element, i.e. in a big mesh how many dofs there are per element without considering the border elements
  //! example: 1D L1-element: 1, L2-element: 2 (the 3rd is already at the next elements)
  static constexpr int averageNDofsPerElement();
  
  //! evaluate the 1D basis function corresponding to element-local dof i at xi, interval for xi is [0,1]
  static double phi(int i, double xi);
  
  //! evaluate the first derivative of the 1D basis function corresponding to element-local dof i at xi, interval for xi is [0,1]
  static double dphi_dxi(int i, double xi);
  
private:
};

}  // namespace

#include "basis_function/lagrange.tpp"