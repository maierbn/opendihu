#include "spatial_discretization/finite_element_method.h"

#include <iostream>
#include <petscksp.h>
#include <vector>
#include <numeric>

#include <Python.h>
#include "easylogging++.h"

#include "control/types.h"
#include "control/python_utility.h"


namespace SpatialDiscretization
{

template<typename Mesh, typename BasisFunction>
FiniteElementMethod<Mesh, BasisFunction, Equation::Static::Laplace>::
FiniteElementMethod(const DihuContext &context)
  : FiniteElementMethodBase<Mesh, BasisFunction>(context)
{
}

template<class MeshType, class BasisFunctionType>
void FiniteElementMethod<MeshType, BasisFunctionType, Equation::Static::Laplace>::
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