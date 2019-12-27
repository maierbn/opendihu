#include "time_stepping_scheme/03_time_stepping_implicit.h"

#include <Python.h>  // has to be the first included header

#include "utility/python_utility.h"
#include "utility/petsc_utility.h"
#include <petscksp.h>
#include "solver/solver_manager.h"
#include "solver/linear.h"
#include "data_management/time_stepping/time_stepping_implicit.h"

namespace TimeSteppingScheme
{

template<typename DiscretizableInTimeType>
TimeSteppingImplicit<DiscretizableInTimeType>::TimeSteppingImplicit(DihuContext context, std::string name) :
TimeSteppingSchemeOde<DiscretizableInTimeType>(context, name)
{
}

template<typename DiscretizableInTimeType>
void TimeSteppingImplicit<DiscretizableInTimeType>::
initialize()
{
  if (this->initialized_)
    return;

  // initialize data objects that are needed for TimeSteppingSchemeOde<DiscretizableInTimeType>::initialize();
  this->data_ = std::make_shared<Data::TimeSteppingImplicit<typename DiscretizableInTimeType::FunctionSpace, DiscretizableInTimeType::nComponents()>>(this->context_); // create data object for implicit euler
  this->dataImplicit_ = std::static_pointer_cast<Data::TimeSteppingImplicit<typename DiscretizableInTimeType::FunctionSpace, DiscretizableInTimeType::nComponents()>>(this->data_);

  TimeSteppingSchemeOde<DiscretizableInTimeType>::initialize();
  LOG(TRACE) << "TimeSteppingImplicit::initialize";

  // compute the system matrix
  this->setSystemMatrix(this->timeStepWidth_);

  LOG(DEBUG) << "time_stepping_implicit applyInSystemMatrix, from TimeSteppingImplicit::initialize";
  // set the boundary conditions to system matrix, i.e. zero rows and columns of Dirichlet BC dofs and set diagonal to 1
  this->dirichletBoundaryConditions_->applyInSystemMatrix(this->dataImplicit_->systemMatrix(), this->dataImplicit_->boundaryConditionsRightHandSideSummand());

  // initialize the linear solver that is used for solving the implicit system
  initializeLinearSolver();
  
  // set matrix used for linear system and preconditioner to ksp context
  Mat &systemMatrix = this->dataImplicit_->systemMatrix()->valuesGlobal();
  assert(this->ksp_);
  PetscErrorCode ierr;
  ierr = KSPSetOperators(*ksp_, systemMatrix, systemMatrix); CHKERRV(ierr);
  
  this->initialized_ = true;
}

template<typename DiscretizableInTimeType>
void TimeSteppingImplicit<DiscretizableInTimeType>::
reset()
{
  TimeSteppingSchemeOdeBaseDiscretizable<DiscretizableInTimeType>::reset();

  LOG(DEBUG) << "set linearSolver_ to nullptr";
  if (linearSolver_)
  {
    LOG(DEBUG) << "delete linear solver";
    this->context_.solverManager()->deleteSolver(linearSolver_->name());
  }

  linearSolver_ = nullptr;
}

template<typename DiscretizableInTimeType>
Data::TimeSteppingImplicit<typename DiscretizableInTimeType::FunctionSpace, DiscretizableInTimeType::nComponents()>
&TimeSteppingImplicit<DiscretizableInTimeType>::
dataImplicit()
{
  return *dataImplicit_;
}

template<typename DiscretizableInTimeType>
void TimeSteppingImplicit<DiscretizableInTimeType>::
solveLinearSystem(Vec &input, Vec &output)
{
  // solve systemMatrix*output = input for output
  Mat &systemMatrix = this->dataImplicit_->systemMatrix()->valuesGlobal();
  
  PetscUtility::checkDimensionsMatrixVector(systemMatrix, input);
  
  if (VLOG_IS_ON(1))
  {
    linearSolver_->solve(input, output, "Linear system of implicit time stepping solved");
  }
  else
  {
    linearSolver_->solve(input, output);
  }
}

template<typename DiscretizableInTimeType>
void TimeSteppingImplicit<DiscretizableInTimeType>::
initializeLinearSolver()
{ 
  LOG(DEBUG) << "initializeLinearSolver, linearSolver_ == nullptr: " << (linearSolver_ == nullptr);
  if (linearSolver_ == nullptr)
  {
    LOG(DEBUG) << "Implicit time stepping: initialize linearSolver";
    
    // retrieve linear solver
    linearSolver_ = this->context_.solverManager()->template solver<Solver::Linear>(
      this->specificSettings_, this->data_->functionSpace()->meshPartition()->mpiCommunicator());
    ksp_ = linearSolver_->ksp();
  }
  else 
  {
    VLOG(2) << ": linearSolver_ already set";
  }
}

/*
//! output the given data for debugging
template<typename DiscretizableInTimeType>
std::string TimeSteppingImplicit<DiscretizableInTimeType>::
getString(typename TimeSteppingSchemeOde<DiscretizableInTimeType>::OutputConnectorDataType &data)
{
  return dataImplicit_->getString(data);
}*/


} // namespace TimeSteppingScheme
