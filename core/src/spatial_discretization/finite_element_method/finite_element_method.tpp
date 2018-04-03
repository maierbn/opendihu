#include "spatial_discretization/finite_element_method/finite_element_method.h"

#include <iostream>
#include <petscksp.h>
#include <vector>
#include <numeric>

#include <Python.h>
#include "easylogging++.h"

#include "control/types.h"
#include "utility/python_utility.h"
#include "equation/type_traits.h"


namespace SpatialDiscretization
{

// set right hand side for laplace equation (set to 0)
template<typename MeshType, typename BasisFunctionType, typename QuadratureType, typename Term>
void FiniteElementMethod<MeshType, BasisFunctionType, QuadratureType, Term, Equation::hasNoRhs<Term>, BasisFunction::notMixed<BasisFunctionType>>::
setRightHandSide()
{
  LOG(TRACE)<<"setRightHandSide for Laplace equation";

  int nDegreesOfFreedom = this->data_.nDegreesOfFreedom();
  
  // fill rhs vector
  PetscErrorCode ierr;
  
  Vec &rightHandSide = this->data_.rightHandSide().values();
  
  for (node_no_t nodeNo = 0; nodeNo < nDegreesOfFreedom; nodeNo++)
  {
   
    //                 vector         row     value
    ierr = VecSetValue(rightHandSide, nodeNo, 0.0, INSERT_VALUES); CHKERRV(ierr);
  }
}

};