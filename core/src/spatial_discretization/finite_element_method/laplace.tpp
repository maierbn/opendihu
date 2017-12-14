#include "spatial_discretization/finite_element_method/finite_element_method.h"

#include <iostream>
#include <petscksp.h>
#include <vector>
#include <numeric>

#include <Python.h>
#include "easylogging++.h"

#include "control/types.h"
#include "utility/python_utility.h"


namespace SpatialDiscretization
{

template<typename MeshType, typename BasisFunctionType, typename IntegratorType>
void FiniteElementMethod<MeshType, BasisFunctionType, IntegratorType, Equation::Static::Laplace>::
setRightHandSide()
{
 LOG(DEBUG)<<"setRightHandSide for Laplace equation";

 int nDegreesOfFreedom = this->data_.nDegreesOfFreedom();
  
 // fill rhs vector
 PetscErrorCode ierr;
 
 Vec &rightHandSide = this->data_.rightHandSide();

 for (node_idx_t nodeNo = 0; nodeNo < nDegreesOfFreedom; nodeNo++)
 {
  //                 vector         row     value
  ierr = VecSetValue(rightHandSide, nodeNo, 0.0, INSERT_VALUES); CHKERRV(ierr);
 }
}

};