#include "spatial_discretization/finite_element_method/05_time_stepping.h"

#include <Python.h>
#include <iostream>
#include <sstream>
#include <petscksp.h>
#include <vector>
#include <numeric>
#include <omp.h>

#include "easylogging++.h"

#include "control/types.h"
#include "utility/python_utility.h"
#include "utility/petsc_utility.h"
#include "solver/solver_manager.h"
#include "solver/linear.h"


namespace SpatialDiscretization
{

template<typename FunctionSpaceType, typename QuadratureType, typename Term>
FiniteElementMethodTimeStepping<FunctionSpaceType, QuadratureType, Term>::
FiniteElementMethodTimeStepping(DihuContext context, std::shared_ptr<FunctionSpaceType> functionSpace)
  : AssembleRightHandSide<FunctionSpaceType, QuadratureType, Term>(context, functionSpace),
  DiscretizableInTime(), Splittable(), linearSolver_(nullptr), ksp_(nullptr)
{
}

template<typename FunctionSpaceType, typename QuadratureType, typename Term>
void FiniteElementMethodTimeStepping<FunctionSpaceType, QuadratureType, Term>::
setBoundaryConditionHandlingEnabled(bool boundaryConditionHandlingEnabled)
{
  BoundaryConditions<FunctionSpaceType,QuadratureType,Term,Term>::setBoundaryConditionHandlingEnabled(boundaryConditionHandlingEnabled);
}

template<typename FunctionSpaceType, typename QuadratureType, typename Term>
void FiniteElementMethodTimeStepping<FunctionSpaceType, QuadratureType, Term>::
initialize()
{
  LOG(DEBUG) << "FiniteElementMethodTimeStepping::initialize";
  
  // call initialize of the parent class
  FiniteElementMethodBase<FunctionSpaceType,QuadratureType,Term>::initialize();

  // initialize the linear solver
  this->initializeLinearSolver();

  // print a warning if this finite element class has output writers, because we do not have solution data to write
  if (this->outputWriterManager_.hasOutputWriters())
  {
    LOG(WARNING) << "You have specified output writers for a FiniteElementMethod which is used for a time stepping problem. "
      "The output will not contain any solution data, only geometry. Probably you want to get output from the time stepping scheme, then define the output writers there.";
  }
}

template<typename FunctionSpaceType, typename QuadratureType, typename Term>
void FiniteElementMethodTimeStepping<FunctionSpaceType, QuadratureType, Term>::
initializeForImplicitTimeStepping()
{
  LOG(DEBUG) << "FiniteElementMethodTimeStepping::initializeForImplicitTimeStepping()";

  // initialize everything needed for implicit time stepping
  // currently this is executed regardless of explicit or implicit time stepping scheme

  // initialize matrices
  this->data_.initializeMassMatrix();
  this->data_.initializeInverseLumpedMassMatrix();

  // compute the mass matrix
  this->setMassMatrix();

  // compute inverse lumped mass matrix
  this->setInverseLumpedMassMatrix();
}

template<typename FunctionSpaceType, typename QuadratureType, typename Term>
void FiniteElementMethodTimeStepping<FunctionSpaceType, QuadratureType, Term>::
reset()
{
  linearSolver_ = nullptr;
  ksp_ = nullptr;
}

//! hook to set initial values for a time stepping from this FiniteElement context, return true if it has set the values or don't do anything and return false
template<typename FunctionSpaceType, typename QuadratureType, typename Term>
bool FiniteElementMethodTimeStepping<FunctionSpaceType, QuadratureType, Term>::
setInitialValues(std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,1>> initialValues)
{
  // Do not set initial values from within the "FiniteElements" section of the config. (therefore return false)
  // The initial values are set by the time stepping scheme under its section.
  return false;
}

template<typename FunctionSpaceType, typename QuadratureType, typename Term>
void FiniteElementMethodTimeStepping<FunctionSpaceType, QuadratureType, Term>::
setRankSubset(Partition::RankSubset rankSubset)
{
  FiniteElementMethodBase<FunctionSpaceType, QuadratureType, Term>::setRankSubset(rankSubset);
}

template<typename FunctionSpaceType, typename QuadratureType, typename Term>
void FiniteElementMethodTimeStepping<FunctionSpaceType, QuadratureType, Term>::
initializeLinearSolver()
{ 
  if (linearSolver_ == nullptr)
  {
    std::stringstream s;
    s << "[" << omp_get_thread_num() << "/" << omp_get_num_threads() << "]";
    LOG(DEBUG) << s.str() << ": FiniteElementMethodTimeStepping: initialize linearSolver";
    
    // retrieve linear solver
    // The mpiCommunicator is needed such that the solver knowns which ranks to use (it could be a subset of all ranks).
    linearSolver_ = this->context_.solverManager()->template solver<Solver::Linear>(
      this->specificSettings_, this->data_.functionSpace()->meshPartition()->mpiCommunicator());
    ksp_ = linearSolver_->ksp();
  }
  else 
  {
    std::stringstream s;
    s << "[" << omp_get_thread_num() << "/" << omp_get_num_threads() << "]";
    VLOG(2) << s.str() << ": linearSolver_ already set";
    
  }
}

template<typename FunctionSpaceType, typename QuadratureType, typename Term>
constexpr int FiniteElementMethodTimeStepping<FunctionSpaceType, QuadratureType, Term>::
nComponents()
{
  return 1;   // this may be different for structural mechanics
}

template<typename FunctionSpaceType, typename QuadratureType, typename Term>
bool FiniteElementMethodTimeStepping<FunctionSpaceType, QuadratureType, Term>::
knowsMeshType()
{
  return true;
}

template<typename FunctionSpaceType, typename QuadratureType, typename Term>
typename FiniteElementMethodTimeStepping<FunctionSpaceType, QuadratureType, Term>::TransferableSolutionDataType
FiniteElementMethodTimeStepping<FunctionSpaceType, QuadratureType, Term>::
getSolutionForTransferInOperatorSplitting()
{
  // check for nans or infs
  //this->data_->solution()->checkNanInf();
  return this->data_->getSolutionForTransferInOperatorSplitting();
}

template<typename FunctionSpaceType, typename QuadratureType, typename Term>
std::shared_ptr<FunctionSpaceType> FiniteElementMethodTimeStepping<FunctionSpaceType, QuadratureType, Term>::
functionSpace()
{
  return FiniteElementMethodBase<FunctionSpaceType, QuadratureType, Term>::functionSpace();
}

} // namespace SpatialDiscretization
