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

  dof_no_t nUnknowns = this->data_.nUnknowns();
  Vec &rightHandSide = this->data_.rightHandSide().values();
  
  std::vector<double> values;
  PythonUtility::getOptionVector(this->specificSettings_, "rightHandSide", nUnknowns, values);
  
#ifndef NDEBUG  
  LOG(DEBUG) << "Read in rhs values from config:";
  std::stringstream s;
  for(auto value : values)
  {
    s << " " << value;
  }
  LOG(DEBUG) << s.str();
#endif
  
  PetscUtility::setVector(values, rightHandSide);
  
  // transform the entries from strong form to weak form
  this->transferRhsToWeakForm();
  
  // if implemented for the current equation, further manipulate the rhs values that are now in weak form
  // this method is empty by default
  this->manipulateWeakRhs();

#ifndef NDEBUG  
  LOG(DEBUG) << "Transferred to weak form:";
  s.str("");
  PetscUtility::getVectorEntries(rightHandSide, values);
  for(auto value : values)
  {
    s << " " << value;
  }
  LOG(DEBUG) << s.str();
#endif
}

};