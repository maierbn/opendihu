#include "data_management/neumann_boundary_conditions.h"

#include "easylogging++.h"

namespace Data
{

template<typename FunctionSpaceType,int nComponents>
NeumannBoundaryConditions<FunctionSpaceType,nComponents>::
~NeumannBoundaryConditions()
{
}

template<typename FunctionSpaceType,int nComponents>
void NeumannBoundaryConditions<FunctionSpaceType,nComponents>::
initialize()
{
  Data<FunctionSpaceType>::initialize();
}

template<typename FunctionSpaceType,int nComponents>
void NeumannBoundaryConditions<FunctionSpaceType,nComponents>::
reset()
{
  // set initalize_ to false
  Data<FunctionSpaceType>::reset();

  // deallocate Petsc vectors
  this->rhs_ = nullptr;
}

template<typename FunctionSpaceType,int nComponents>
void NeumannBoundaryConditions<FunctionSpaceType,nComponents>::
createPetscObjects()
{
  LOG(TRACE) << "NeumannBoundaryConditions::createPetscObjects";

  assert(this->functionSpace_);
  this->rhs_ = this->functionSpace_->template createFieldVariable<nComponents>("-rhs_NeumannBC");
}

template<typename FunctionSpaceType,int nComponents>
std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,nComponents>> NeumannBoundaryConditions<FunctionSpaceType,nComponents>::
rhs()
{
  return this->rhs_;
}

} // namespace Data
