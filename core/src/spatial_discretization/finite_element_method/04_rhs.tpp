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
  LOG(TRACE)<<"setRightHandSide";

  dof_no_t nLocalUnknowns = this->data_.nLocalUnknowns();
  FieldVariable::FieldVariable<BasisOnMeshType,1> &rightHandSide = this->data_.rightHandSide();

  std::vector<double> localValues;
  
  bool inputMeshIsGlobal = PythonUtility::getOptionBool(this->specificSettings_, "inputMeshIsGlobal", true);
  if (inputMeshIsGlobal)
  {
    global_no_t nGlobalUnknowns = this->data_.nGlobalUnknowns();
    PythonUtility::getOptionVector(this->specificSettings_, "rightHandSide", (int)nGlobalUnknowns, localValues);

    std::shared_ptr<Mesh::Mesh> mesh = this->data_.mesh();
    mesh->meshPartitionBase()->extractLocalDofs(localValues);
  }
  else 
  {
    PythonUtility::getOptionVector(this->specificSettings_, "rightHandSide", nLocalUnknowns, localValues);
  }
  
#ifndef NDEBUG
  LOG(DEBUG) << "Read in rhs localValues from config:";
  std::stringstream s;
  for(auto value : localValues)
  {
    s << " " << value;
  }
  LOG(DEBUG) << s.str();
#endif

  rightHandSide.setValues(localValues);

  // transform the entries from strong form to weak form
  this->transferRhsToWeakForm();

  // if implemented for the current equation, further manipulate the rhs values that are now in weak form
  // this method is empty by default
  this->manipulateWeakRhs();

#ifndef NDEBUG
  LOG(DEBUG) << "Transferred to weak form:";
  s.str("");
  PetscUtility::getVectorEntries(rightHandSide.values(), localValues);
  for(auto value : localValues)
  {
    s << " " << value;
  }
  LOG(DEBUG) << s.str();
#endif
}

};