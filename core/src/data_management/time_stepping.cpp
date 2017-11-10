#include "data_management/time_stepping.h"

#include <iomanip>
#include <iostream>
#include <sstream>
#include <vector>
#include <numeric>

#include "easylogging++.h"

#include "control/python_utility.h"
#include "control/dihu_context.h"
#include "control/petsc_utility.h"

namespace Data 
{

TimeStepping::TimeStepping(DihuContext &context) : Data(context)
{
  disablePrinting_ = PythonUtility::getOptionBool(context_.getPythonConfig(), "disablePrinting", false);
}

TimeStepping::~TimeStepping()
{
  PetscErrorCode ierr;
  // free PETSc objects
  if (initialized_)
  {
    ierr = VecDestroy(&solution_); CHKERRV(ierr);
  }
}

void TimeStepping::createPetscObjects()
{
  element_idx_t n = mesh_->nDegreesOfFreedom();
  
  LOG(DEBUG)<<"TimeStepping::createPetscObjects("<<n<<")"<<std::endl;
  
  PetscErrorCode ierr;
  // create PETSc vector object
  ierr = VecCreate(PETSC_COMM_WORLD, &solution_);  CHKERRV(ierr);
  ierr = PetscObjectSetName((PetscObject) solution_, "solution"); CHKERRV(ierr);
  
  // initialize size of vector
  ierr = VecSetSizes(solution_, PETSC_DECIDE, n); CHKERRV(ierr);
  
  // set sparsity type and other options
  ierr = VecSetFromOptions(solution_);  CHKERRV(ierr);
  
  // create the other vector
  ierr = VecDuplicate(solution_, &increment_);  CHKERRV(ierr);
}

Vec& TimeStepping::solution()
{
  return solution_;
}

Vec& TimeStepping::increment()
{
  return increment_;
}

void TimeStepping::print()
{
  if (disablePrinting_)
    return;
  
  LOG(INFO)<<"======================";
  
  int nEntries;
  VecGetSize(increment_, &nEntries);
  LOG(INFO)<<"increment ("<<nEntries<<" entries):";
  LOG(INFO)<<PetscUtility::getStringVector(increment_);
  LOG(INFO)<<"======================";
  
  VecGetSize(solution_, &nEntries);
  LOG(INFO)<<"solution ("<<nEntries<<" entries):";
  LOG(INFO)<<PetscUtility::getStringVector(solution_);
  LOG(INFO)<<"======================";
}



} // namespace