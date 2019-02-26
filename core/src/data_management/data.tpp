#include "data_management/data.h"

#include <iomanip>
#include <iostream>
#include <sstream>
#include <vector>
#include <numeric>

#include "easylogging++.h"

#include "utility/python_utility.h"
#include "control/dihu_context.h"

namespace Data
{

template<typename FunctionSpaceType>
Data<FunctionSpaceType>::
Data(DihuContext context) :
  context_(context), rankSubset_(nullptr)
{
}

template<typename FunctionSpaceType>
Data<FunctionSpaceType>::
~Data()
{
}

template<typename FunctionSpaceType>
dof_no_t Data<FunctionSpaceType>::
nUnknownsLocalWithGhosts()
{
  return this->functionSpace_->nDofsLocalWithGhosts() * 1;  // value for 1 component, can be overloaded to also have the factor nComponents in it
}

template<typename FunctionSpaceType>
dof_no_t Data<FunctionSpaceType>::
nUnknownsLocalWithoutGhosts()
{
  return this->functionSpace_->nDofsLocalWithoutGhosts() * 1;  // value for 1 component, can be overloaded to also have the factor nComponents in it
}

template<typename FunctionSpaceType>
global_no_t Data<FunctionSpaceType>::
nUnknownsGlobal()
{
  return this->functionSpace_->nDofsGlobal() * 1;  // value for 1 component, can be overloaded to also have the factor nComponents in it
}

template<typename FunctionSpaceType>
void Data<FunctionSpaceType>::
setFunctionSpace(std::shared_ptr<FunctionSpaceType> functionSpace)
{
  this->functionSpace_ = functionSpace;
}

template<typename FunctionSpaceType>
void Data<FunctionSpaceType>::
initialize()
{
  if (this->initialized_)
    return;

  this->createPetscObjects();
  this->initialized_ = true;
}

template<typename FunctionSpaceType>
void Data<FunctionSpaceType>::
reset()
{
  this->initialized_ = false;
}

template<typename FunctionSpaceType>
void Data<FunctionSpaceType>::
setRankSubset(Partition::RankSubset rankSubset)
{
  rankSubset_ = std::make_shared<Partition::RankSubset>(rankSubset);
}
  
template<typename FunctionSpaceType>
const std::shared_ptr<FunctionSpaceType> Data<FunctionSpaceType>::
functionSpace() const
{
  return this->functionSpace_;
}

} // namespace Data
