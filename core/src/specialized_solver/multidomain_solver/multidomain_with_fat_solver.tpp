#include "specialized_solver/multidomain_solver/multidomain_with_fat_solver.h"

#include <Python.h>  // has to be the first included header

namespace TimeSteppingScheme
{

template<typename FiniteElementMethodPotentialFlow,typename FiniteElementMethodDiffusionMuscle,typename FiniteElementMethodDiffusionFat>
MultidomainWithFatSolver<FiniteElementMethodPotentialFlow,FiniteElementMethodDiffusionMuscle,FiniteElementMethodDiffusionFat>::
MultidomainWithFatSolver(DihuContext context) :
  MultidomainSolver<FiniteElementMethodPotentialFlow,FiniteElementMethodDiffusionMuscle>(context),
  dataFat_(this->context_),
  finiteElementMethodFat_(this->context_["Fat"])
{
}


template<typename FiniteElementMethodPotentialFlow,typename FiniteElementMethodDiffusionMuscle,typename FiniteElementMethodDiffusionFat>
void MultidomainWithFatSolver<FiniteElementMethodPotentialFlow,FiniteElementMethodDiffusionMuscle,FiniteElementMethodDiffusionFat>::
initialize()
{
  MultidomainSolver<FiniteElementMethodPotentialFlow,FiniteElementMethodDiffusionMuscle>::initialize();

  // indicate in solverStructureVisualizer that now a child solver will be initialized
  DihuContext::solverStructureVisualizer()->beginChild("Fat");

  LOG(DEBUG) << "initialize fat FEM";
  // initialize the potential flow finite element method, this also creates the function space
  finiteElementMethodFat_.initialize();

  // indicate in solverStructureVisualizer that the child solver initialization is done
  DihuContext::solverStructureVisualizer()->endChild();

  LOG(DEBUG) << "function space \"" << finiteElementMethodFat_.functionSpace()->meshName() << "\".";

  // initialize the data object
  dataFat_.setFunctionSpace(finiteElementMethodFat_.functionSpace());
  dataFat_.setDataMultidomain(
    std::make_shared<typename MultidomainSolver<FiniteElementMethodPotentialFlow,FiniteElementMethodDiffusionMuscle>::Data>(
    this->dataMultidomain_));
  
  dataFat_.initialize();

  DihuContext::solverStructureVisualizer()->setOutputConnectorData(this->getOutputConnectorData());
}

template<typename FiniteElementMethodPotentialFlow,typename FiniteElementMethodDiffusionMuscle,typename FiniteElementMethodDiffusionFat>
void MultidomainWithFatSolver<FiniteElementMethodPotentialFlow,FiniteElementMethodDiffusionMuscle,FiniteElementMethodDiffusionFat>::
callOutputWriter(int timeStepNo, double currentTime)
{
  // write current output values
  this->outputWriterManager_.writeOutput(this->dataFat_, timeStepNo, currentTime);
}

template<typename FiniteElementMethodPotentialFlow,typename FiniteElementMethodDiffusionMuscle,typename FiniteElementMethodDiffusionFat>
void MultidomainWithFatSolver<FiniteElementMethodPotentialFlow,FiniteElementMethodDiffusionMuscle,FiniteElementMethodDiffusionFat>::
setSystemMatrix(double timeStepWidth)
{
  MultidomainSolver<FiniteElementMethodPotentialFlow,FiniteElementMethodDiffusionMuscle>::setSystemMatrix(timeStepWidth);
}

template<typename FiniteElementMethodPotentialFlow,typename FiniteElementMethodDiffusionMuscle,typename FiniteElementMethodDiffusionFat>
void MultidomainWithFatSolver<FiniteElementMethodPotentialFlow,FiniteElementMethodDiffusionMuscle,FiniteElementMethodDiffusionFat>::
solveLinearSystem()
{
  VLOG(1) << "in solveLinearSystem";

  // configure that the initial value for the iterative solver is the value in solution, not zero
  PetscErrorCode ierr;
  if (this->initialGuessNonzero_)
  {
    LOG(DEBUG) << "set initial guess nonzero";
    ierr = KSPSetInitialGuessNonzero(*this->linearSolver_->ksp(), PETSC_TRUE); CHKERRV(ierr);
  }

  // copy the values from a nested Petsc Vec to a single Vec that contains all entries
  NestedMatVecUtility::createVecFromNestedVec(this->nestedRightHandSide_, this->singleRightHandSide_, data().functionSpace()->meshPartition()->rankSubset());

  // copy the values from a nested Petsc Vec to a single Vec that contains all entries
  NestedMatVecUtility::createVecFromNestedVec(this->nestedSolution_, this->singleSolution_, data().functionSpace()->meshPartition()->rankSubset());

  // solve the linear system
  // this can be done using the nested Vecs and nested Mat (nestedSolution_, nestedRightHandSide_, nestedSystemMatrix_),
  // or the single Vecs and Mats that contain all values directly  (singleSolution_, singleRightHandSide_, singleSystemMatrix_) 


  if (this->showLinearSolverOutput_)
  {
    this->linearSolver_->solve(this->singleRightHandSide_, this->singleSolution_, "Linear system of multidomain problem solved");
  }
  else
  {
    // solve without showing output
    this->linearSolver_->solve(this->singleRightHandSide_, this->singleSolution_);
  }
  
  // copy the values back from a single Vec that contains all entries to a nested Petsc Vec
  NestedMatVecUtility::fillNestedVec(this->singleSolution_, this->nestedSolution_);
}

template<typename FiniteElementMethodPotentialFlow,typename FiniteElementMethodDiffusionMuscle,typename FiniteElementMethodDiffusionFat>
typename MultidomainWithFatSolver<FiniteElementMethodPotentialFlow,FiniteElementMethodDiffusionMuscle,FiniteElementMethodDiffusionFat>::DataFat &MultidomainWithFatSolver<FiniteElementMethodPotentialFlow,FiniteElementMethodDiffusionMuscle,FiniteElementMethodDiffusionFat>::
data()
{
  return dataFat_;
}

} // namespace TimeSteppingScheme
