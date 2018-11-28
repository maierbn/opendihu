#include "field_variable/00_field_variable_base.h"

namespace FieldVariable
{

template<typename FunctionSpaceType>
FieldVariableBase<FunctionSpaceType>::
FieldVariableBase() : functionSpace_(nullptr)
{
}

template<typename FunctionSpaceType>
std::shared_ptr<FunctionSpaceType> FieldVariableBase<FunctionSpaceType>::
functionSpace()
{
  return functionSpace_;
}

template<typename FunctionSpaceType>
std::string FieldVariableBase<FunctionSpaceType>::
name() const
{
  return this->name_;
}

template<typename FunctionSpaceType>
bool FieldVariableBase<FunctionSpaceType>::
isGeometryField() const
{
  return this->isGeometryField_;
}

template<typename FunctionSpaceType>
void FieldVariableBase<FunctionSpaceType>::
checkNansInfs(int componentNo) const
{
  // get all local values without ghosts for the given componentNo
  std::vector<double> values;
  getValuesWithoutGhosts(componentNo, values);

  // determine if there are nans or high values
  int nNans = 0;
  int nHighValues = 0;
  for (int i = 0; i < values.size(); i++)
  {
    if (std::isnan(values[i]))
      nNans++;
    else if (fabs(values[i]) > 1e100)
      nHighValues++;
  }

  if (nNans > 0)
  {
    LOG(WARNING) << "Solution contains " << nNans << " Nans, out of " << values.size() << " total values";
  }

  if (nHighValues > 0)
  {
    LOG(WARNING) << "Solution contains " << nHighValues << " high values with absolute value > 1e100, out of " << values.size() << " total values";
  }

  if (nNans == values.size())
  {
    LOG(FATAL) << "There are only Nans, abort computation.";
  }
}

//! get the number of dofs
template<typename FunctionSpaceType>
dof_no_t FieldVariableBase<FunctionSpaceType>::
nDofsLocalWithoutGhosts() const
{
  return this->functionSpace_->meshPartition()->nDofsLocalWithoutGhosts();
}

//! get the number of dofs
template<typename FunctionSpaceType>
dof_no_t FieldVariableBase<FunctionSpaceType>::
nDofsGlobal() const
{
  return this->functionSpace_->meshPartition()->nDofsGlobal();
}

} // namespace
