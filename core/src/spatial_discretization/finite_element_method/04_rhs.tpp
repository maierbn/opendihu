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

template<typename BasisOnMeshType, typename QuadratureType, typename Term>
void FiniteElementMethodRhs<BasisOnMeshType, QuadratureType, Term>::
setRightHandSide()
{
  LOG(TRACE) << "setRightHandSide";

  dof_no_t nUnknownsLocal = this->data_.nUnknownsLocalWithoutGhosts();     // local unknows without ghosts
  FieldVariable::FieldVariable<BasisOnMeshType,1> &rightHandSide = this->data_.rightHandSide();
  rightHandSide.startVectorManipulation();

  std::vector<double> localValues;
  
  // parse values from config
  bool inputMeshIsGlobal = PythonUtility::getOptionBool(this->specificSettings_, "inputMeshIsGlobal", true);
  if (inputMeshIsGlobal)
  {
    global_no_t nUnknownsGlobal = this->data_.nUnknownsGlobal();
    PythonUtility::getOptionVector(this->specificSettings_, "rightHandSide", (int)nUnknownsGlobal, localValues);

    std::shared_ptr<Mesh::Mesh> mesh = this->data_.mesh();
    mesh->meshPartitionBase()->extractLocalDofsWithoutGhosts(localValues);
  }
  else 
  {
    PythonUtility::getOptionVector(this->specificSettings_, "rightHandSide", nUnknownsLocal, localValues);
  }

  // assign read values to rightHandSide variable. They are now stored in "strong form" and need to be transformed to "weak form" by multiplying with mass matrix
  rightHandSide.setValuesWithoutGhosts(localValues);
  rightHandSide.finishVectorManipulation();

  // transform the entries from strong form to weak form
  this->multiplyRightHandSideWithMassMatrix();

  // if implemented for the current equation, further manipulate the rhs values that are now in weak form
  // this method is empty by default
  this->manipulateWeakRhs();
}

};
