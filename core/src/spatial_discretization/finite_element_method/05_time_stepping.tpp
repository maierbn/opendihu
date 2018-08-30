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
FiniteElementMethodTimeStepping(DihuContext context)
  : AssembleRightHandSide<FunctionSpaceType, QuadratureType, Term>(context),
  DiscretizableInTime(SolutionVectorMapping(true)), linearSolver_(nullptr), ksp_(nullptr)
{
  // The solutionVectorMapping_ object stores the information which range of values of the solution will be further used
  // in methods that use the result of this method, e.g. in operator splittings. Since there are no internal values
  // in this FEM, set the range to all values.
  solutionVectorMapping_.setOutputRange(0, this->data_.functionSpace()->nNodesLocalWithoutGhosts());   // without ghosts because CellML vectors do not have ghost nodes
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
initialize(double timeStepWidth)
{
  LOG(DEBUG) << "FiniteElementMethodTimeStepping::initialize(timeStepWidth=" << timeStepWidth << ")";

  // initialize everything needed for implicit time stepping
  // currently this is executed regardless of explicit or implicit time stepping scheme

  // initialize matrices
  this->data_.initializeMassMatrix();
  this->data_.initializeInverseLumpedMassMatrix();

  // compute the mass matrix
  this->setMassMatrix();

  // compute inverse lumped mass matrix
  setInverseLumpedMassMatrix();

  // initialize and compute the system matrix
  setSystemMatrix(timeStepWidth);
}

template<typename FunctionSpaceType, typename QuadratureType, typename Term>
void FiniteElementMethodTimeStepping<FunctionSpaceType, QuadratureType, Term>::
reset()
{
  linearSolver_ = nullptr;
  ksp_ = nullptr;
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
    linearSolver_ = this->context_.solverManager()->template solver<Solver::Linear>(this->specificSettings_);
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
std::shared_ptr<Mesh::Mesh> FiniteElementMethodTimeStepping<FunctionSpaceType, QuadratureType, Term>::
mesh()
{
  return FiniteElementMethodBase<FunctionSpaceType, QuadratureType, Term>::mesh();
}

} // namespace SpatialDiscretization
