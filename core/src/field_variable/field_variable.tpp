#include "field_variable/field_variable.h"

namespace FieldVariable
{
 
//! this has to be called before the vector is manipulated (i.e. VecSetValues or vecZeroEntries is called), to ensure that the current state of the vector is fetched from the global vector
template<typename FunctionSpaceType,int nComponents>
void FieldVariable<FunctionSpaceType,nComponents>::
startGhostManipulation()
{
  this->values_->startGhostManipulation();
}
  
//! this has to be called after the vector is manipulated (i.e. VecSetValues or vecZeroEntries is called), to ensure that operations on different partitions are merged by Petsc
template<typename FunctionSpaceType,int nComponents>
void FieldVariable<FunctionSpaceType,nComponents>::
finishGhostManipulation()
{
  this->values_->finishGhostManipulation();
}

};  // namespace
