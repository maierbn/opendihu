#include "field_variable/field_variable.h"

#include <cmath>

namespace FieldVariable
{
 
//! this has to be called before the vector is manipulated (i.e. VecSetValues or vecZeroEntries is called), to ensure that the current state of the vector is fetched from the global vector
template<typename FunctionSpaceType,int nComponents>
void FieldVariable<FunctionSpaceType,nComponents>::
startGhostManipulation()
{
  // if there is an internal values_ vector (this is not the case for geometry fields of stencil-type settings)
  if (this->values_)
    this->values_->startGhostManipulation();
}

//! this has to be called after the vector is manipulated (i.e. VecSetValues or vecZeroEntries is called), to ensure that operations on different partitions are merged by Petsc
template<typename FunctionSpaceType,int nComponents>
void FieldVariable<FunctionSpaceType,nComponents>::
zeroGhostBuffer()
{
  if (this->values_) // if there is an internal values_ vector (this is not the case for geometry fields of stencil-type settings)
    this->values_->zeroGhostBuffer();
}

//! this has to be called after the vector is manipulated (i.e. VecSetValues or vecZeroEntries is called), to ensure that operations on different partitions are merged by Petsc
template<typename FunctionSpaceType,int nComponents>
void FieldVariable<FunctionSpaceType,nComponents>::
finishGhostManipulation()
{
  if (this->values_) // if there is an internal values_ vector (this is not the case for geometry fields of stencil-type settings)
    this->values_->finishGhostManipulation();
}

//! this has to be called after the vector is manipulated (i.e. VecSetValues or vecZeroEntries is called), to ensure that operations on different partitions are merged by Petsc
template<typename FunctionSpaceType,int nComponents>
void FieldVariable<FunctionSpaceType,nComponents>::
setRepresentationGlobal()
{
  if (this->values_) // if there is an internal values_ vector (this is not the case for geometry fields of stencil-type settings)
    this->values_->setRepresentationGlobal();
}

//! this has to be called after the vector is manipulated (i.e. VecSetValues or vecZeroEntries is called), to ensure that operations on different partitions are merged by Petsc
template<typename FunctionSpaceType,int nComponents>
void FieldVariable<FunctionSpaceType,nComponents>::
setRepresentationLocal()
{
  if (this->values_) // if there is an internal values_ vector (this is not the case for geometry fields of stencil-type settings)
    this->values_->setRepresentationLocal();
}

//! this has to be called after the vector is manipulated (i.e. VecSetValues or vecZeroEntries is called), to ensure that operations on different partitions are merged by Petsc
template<typename FunctionSpaceType,int nComponents>
void FieldVariable<FunctionSpaceType,nComponents>::
setRepresentationContiguous()
{
  if (this->values_) // if there is an internal values_ vector (this is not the case for geometry fields of stencil-type settings)
    this->values_->setRepresentationContiguous();
}

template<typename FunctionSpaceType,int nComponents>
bool FieldVariable<FunctionSpaceType,nComponents>::
containsNanOrInf()
{
  if (this->values_)
  {
    // get all local values
    std::vector<VecD<nComponents>> values;
    this->getValuesWithoutGhosts(values);

    // loop over values and check if they are neither nan nor inf
    for (int i = 0; i < values.size(); i++)
    {
      for (int componentNo = 0; componentNo < nComponents; componentNo++)
      {
        if (!std::isfinite(values[i][componentNo]))
        {
          LOG(ERROR) << "containsNanOrInf(): value " << i << "/" << values.size() << ", component " << componentNo << "/" << nComponents << ": " << values[componentNo][i];
          return true;
        }
      }
    }
  }
  return false;
}

} // namespace
