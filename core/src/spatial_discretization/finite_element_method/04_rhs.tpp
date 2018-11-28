#include "spatial_discretization/finite_element_method/04_rhs.h"

#include <iostream>
#include <petscksp.h>
#include <vector>
#include <numeric>

#include <Python.h>
#include "easylogging++.h"

#include "control/types.h"
#include "utility/python_utility.h"
#include "utility/petsc_utility.h"

namespace SpatialDiscretization
{

template<typename FunctionSpaceType, typename QuadratureType, typename Term>
void FiniteElementMethodRhs<FunctionSpaceType, QuadratureType, Term>::
setRightHandSide()
{
  LOG(TRACE) << "setRightHandSide";

  dof_no_t nUnknownsLocal = this->data_.nUnknownsLocalWithoutGhosts();     // local unknows without ghosts
  std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,1>> rightHandSide = this->data_.rightHandSide();

  std::vector<double> localValues;
  
  // parse values from config
  bool inputMeshIsGlobal = this->specificSettings_.getOptionBool("inputMeshIsGlobal", true);
  if (inputMeshIsGlobal)
  {
    global_no_t nUnknownsGlobal = this->data_.nUnknownsGlobal();
    this->specificSettings_.getOptionVector("rightHandSide", (int)nUnknownsGlobal, localValues);

    std::shared_ptr<Mesh::Mesh> functionSpace = this->data_.functionSpace();
    functionSpace->meshPartitionBase()->extractLocalDofsWithoutGhosts(localValues);
  }
  else 
  {
    this->specificSettings_.getOptionVector("rightHandSide", nUnknownsLocal, localValues);
  }

  // assign read values to rightHandSide variable. They are now stored in "strong form" and need to be transformed to "weak form" by multiplying with mass matrix
  rightHandSide->setValuesWithoutGhosts(localValues);

  // transform the entries from strong form to weak form
  this->multiplyRightHandSideWithMassMatrix();

  // if implemented for the current equation, further manipulate the rhs values that are now in weak form
  // this method is empty by default
  this->manipulateWeakRhs();
}

};
