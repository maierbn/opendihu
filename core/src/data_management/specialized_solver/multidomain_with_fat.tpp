#include "data_management/specialized_solver/multidomain_with_fat.h"

#include <iomanip>
#include <iostream>
#include <sstream>
#include <vector>
#include <numeric>
#include <memory>

#include "easylogging++.h"

#include "utility/python_utility.h"
#include "control/dihu_context.h"
#include "utility/petsc_utility.h"

namespace Data
{

template<typename FunctionSpaceType, typename FunctionSpaceFatType>
MultidomainWithFat<FunctionSpaceType,FunctionSpaceFatType>::
MultidomainWithFat(DihuContext context) :
  Data<FunctionSpaceType>::Data(context)
{
}

template<typename FunctionSpaceType, typename FunctionSpaceFatType>
void MultidomainWithFat<FunctionSpaceType,FunctionSpaceFatType>::
createPetscObjects()
{
  LOG(DEBUG) << "MultidomainWithFat::createPetscObjects";

  this->extraCellularPotentialFat_ = this->functionSpace_->template createFieldVariable<1>("phi_b");
}

//! initialize the function space
template<typename FunctionSpaceType, typename FunctionSpaceFatType>
void MultidomainWithFat<FunctionSpaceType,FunctionSpaceFatType>::
setDataMultidomain(std::shared_ptr<Multidomain<FunctionSpaceType>> dataMultidomain)
{
  this->dataMultidomain_ = dataMultidomain;
}

template<typename FunctionSpaceType, typename FunctionSpaceFatType>
std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceFatType,1>> MultidomainWithFat<FunctionSpaceType,FunctionSpaceFatType>::
extraCellularPotentialFat()
{
  return this->extraCellularPotentialFat_;
}

template<typename FunctionSpaceType, typename FunctionSpaceFatType>
typename MultidomainWithFat<FunctionSpaceType,FunctionSpaceFatType>::FieldVariablesForOutputWriter MultidomainWithFat<FunctionSpaceType,FunctionSpaceFatType>::
getFieldVariablesForOutputWriter()
{
  // these field variables will be written to output files
  std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceFatType,3>> geometryFieldFat
    = std::make_shared<FieldVariable::FieldVariable<FunctionSpaceFatType,3>>(this->functionSpace_->geometryField());

  LOG(DEBUG) << "getFieldVariablesForOutputWriter: this->extraCellularPotentialFat_: " << *this->extraCellularPotentialFat_;

  return std::tuple_cat(dataMultidomain_->getFieldVariablesForOutputWriter(), std::make_tuple(geometryFieldFat, this->extraCellularPotentialFat_));
}

} // namespace
