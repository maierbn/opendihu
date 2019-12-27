#include "data_management/time_stepping/time_stepping.h"

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
#include "partition/mesh_partition/01_mesh_partition.h"
#include "partition/partitioned_petsc_mat/partitioned_petsc_mat.h"

namespace Data
{

template<typename FunctionSpaceType,int nComponents>
TimeStepping<FunctionSpaceType,nComponents>::
TimeStepping(DihuContext context) : Data<FunctionSpaceType>(context)
{
  this->debuggingName_ = "timestepping";
}

template<typename FunctionSpaceType,int nComponents>
TimeStepping<FunctionSpaceType,nComponents>::
~TimeStepping()
{
  // free PETSc objects
  if (this->initialized_)
  {
    //PetscErrorCode ierr;
    //ierr = VecDestroy(&solution_); CHKERRV(ierr);
  }
}

template<typename FunctionSpaceType,int nComponents>
void TimeStepping<FunctionSpaceType,nComponents>::
createPetscObjects()
{
  LOG(DEBUG) << "TimeStepping<FunctionSpaceType,nComponents>::createPetscObjects(" <<nComponents << ")";
  assert(this->functionSpace_);
  assert(this->functionSpace_->meshPartition());

  if (componentNames_.empty())
  {
    this->solution_ = this->functionSpace_->template createFieldVariable<nComponents>("solution");
    this->increment_ = this->functionSpace_->template createFieldVariable<nComponents>("increment");
  }
  else 
  {
    // if there are component names stored, use them for construction of the field variables 
    this->solution_ = this->functionSpace_->template createFieldVariable<nComponents>("solution", componentNames_);
    this->increment_ = this->functionSpace_->template createFieldVariable<nComponents>("increment", componentNames_);
  }

  outputConnectorData_.addFieldVariable(this->solution_);

  // create additional field variables that appear as output connector slots and can be connected to discretizableInTime_ and surrounding solvers
  int nAdditionalFieldVariables = this->context_.getPythonConfig().getOptionInt("nAdditionalFieldVariables", 0, PythonUtility::NonNegative);
  additionalFieldVariables_.resize(nAdditionalFieldVariables);

  for (int i = 0; i < nAdditionalFieldVariables; i++)
  {
    std::stringstream name;
    name << "additionalFieldVariable" << i;
    additionalFieldVariables_[i] = this->functionSpace_->template createFieldVariable<1>(name.str());

    outputConnectorData_.addFieldVariable2(additionalFieldVariables_[i]);
  }
}

template<typename FunctionSpaceType,int nComponents>
std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,nComponents>> TimeStepping<FunctionSpaceType,nComponents>::
solution()
{
  return this->solution_;
}

template<typename FunctionSpaceType,int nComponents>
std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,nComponents>> TimeStepping<FunctionSpaceType,nComponents>::
increment()
{
  return this->increment_;
}

template<typename FunctionSpaceType,int nComponents>
dof_no_t TimeStepping<FunctionSpaceType,nComponents>::
nUnknownsLocalWithGhosts()
{
  return this->functionSpace_->nNodesLocalWithGhosts() * nComponents;
}

template<typename FunctionSpaceType,int nComponents>
dof_no_t TimeStepping<FunctionSpaceType,nComponents>::
nUnknownsLocalWithoutGhosts()
{
  return this->functionSpace_->nNodesLocalWithoutGhosts() * nComponents;
}

template<typename FunctionSpaceType,int nComponents>
constexpr int TimeStepping<FunctionSpaceType,nComponents>::
getNDofsPerNode()
{
  return nComponents;
}

template<typename FunctionSpaceType,int nComponents>
void TimeStepping<FunctionSpaceType,nComponents>::
print()
{
  if (!VLOG_IS_ON(4))
    return;

  VLOG(4) << "======================";
  VLOG(4) << *this->increment_;
  VLOG(4) << *this->solution_;
  VLOG(4) << "======================";
}

template<typename FunctionSpaceType,int nComponents>
void TimeStepping<FunctionSpaceType,nComponents>::
setComponentNames(std::vector<std::string> componentNames)
{
  componentNames_ = componentNames;
}
/*
template<typename FunctionSpaceType,int nComponents>
void TimeStepping<FunctionSpaceType,nComponents>::
setOutputComponentNo(int outputComponentNo)
{
  outputComponentNo_ = outputComponentNo;
}

template<typename FunctionSpaceType,int nComponents>
void TimeStepping<FunctionSpaceType,nComponents>::
setPrefactor(double prefactor)
{
  prefactor_ = prefactor;
}*/

template<typename FunctionSpaceType,int nComponents>
typename TimeStepping<FunctionSpaceType,nComponents>::OutputConnectorDataType &TimeStepping<FunctionSpaceType,nComponents>::
getOutputConnectorData()
{
  return outputConnectorData_;
}

template<typename FunctionSpaceType,int nComponents>
typename TimeStepping<FunctionSpaceType,nComponents>::FieldVariablesForOutputWriter TimeStepping<FunctionSpaceType,nComponents>::
getFieldVariablesForOutputWriter()
{
  // these field variables will be written to output files
  return FieldVariablesForOutputWriter(
    std::make_shared<FieldVariable::FieldVariable<FunctionSpaceType,3>>(this->functionSpace_->geometryField()),
    solution_,
    additionalFieldVariables_
  );
}

//! output the given data for debugging
template<typename FunctionSpaceType,int nComponents>
std::string TimeStepping<FunctionSpaceType,nComponents>::
getString(typename TimeStepping<FunctionSpaceType,nComponents>::OutputConnectorDataType &data)
{
  std::stringstream s;
  s << "<" << debuggingName_ << ":" << data << ">";

  return s.str();
}



} // namespace
