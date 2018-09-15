#include "time_stepping_scheme/multidomain_solver.h"

#include <Python.h>  // has to be the first included header

#include "utility/python_utility.h"
#include "utility/petsc_utility.h"
#include "data_management/multidomain.h"

namespace TimeSteppingScheme
{

template<typename FiniteElementMethodPotentialFlow,typename CellMLAdapter,typename FiniteElementMethodDiffusion>
MultidomainSolver<FiniteElementMethodPotentialFlow,CellMLAdapter,FiniteElementMethodDiffusion>::
MultidomainSolver(DihuContext context) :
  TimeSteppingImplicit<FiniteElementMethodDiffusion>(context, "MultidomainSolver"),
  finiteElementMethodPotentialFlow_(context),
  cellMLAdapter_(context)
{
  this->data_ = std::make_shared<Data::Multidomain<FiniteElementMethodDiffusion>>(context); // create data object for implicit euler
}

template<typename FiniteElementMethodPotentialFlow,typename CellMLAdapter,typename FiniteElementMethodDiffusion>
void MultidomainSolver<FiniteElementMethodPotentialFlow,CellMLAdapter,FiniteElementMethodDiffusion>::
advanceTimeSpan()
{
  // compute timestep width
  double timeSpan = this->endTime_ - this->startTime_;

  LOG(DEBUG) << "MultidomainSolver::advanceTimeSpan, timeSpan=" << timeSpan<< ", timeStepWidth=" << this->timeStepWidth_
    << " n steps: " << this->numberTimeSteps_;

  Vec solution = this->data_->solution()->valuesGlobal();

  // loop over time steps
  double currentTime = this->startTime_;

  for(int timeStepNo = 0; timeStepNo < this->numberTimeSteps_;)
  {
    if (timeStepNo % this->timeStepOutputInterval_ == 0 && timeStepNo > 0)
    {
      LOG(INFO) << "Multidomain solver, timestep " << timeStepNo << "/" << this->numberTimeSteps_<< ", t=" << currentTime;
    }

    //VLOG(1) << "initial solution: " << *this->data_->solution();

    // advance simulation time
    timeStepNo++;
    currentTime = this->startTime_ + double(timeStepNo) / this->numberTimeSteps_ * timeSpan;

    // adjust rhs vector such that boundary conditions are satisfied
    this->dirichletBoundaryConditions_->applyInRightHandSide(this->data_->solution(), this->dataImplicit_->boundaryConditionsRightHandSideSummand());

    //VLOG(1) << "solution after apply BC: " << *this->data_->solution();

    // advance computed value
    // solve A*u^{t+1} = u^{t} for u^{t+1} where A is the system matrix, solveLinearSystem(b,x)
    this->solveLinearSystem(solution, solution);

    VLOG(1) << "new solution: " << *this->data_->solution();

    // write current output values
    this->outputWriterManager_.writeOutput(*this->data_, timeStepNo, currentTime);

    //this->data_->print();
  }
}

template<typename FiniteElementMethodPotentialFlow,typename CellMLAdapter,typename FiniteElementMethodDiffusion>
void MultidomainSolver<FiniteElementMethodPotentialFlow,CellMLAdapter,FiniteElementMethodDiffusion>::
run()
{
  // initialize everything
  initialize();

  // solve potential flow Laplace problem
  finiteElementMethodPotentialFlow_.run();

  // compute a gradient field from the solution of the potential flow
  data_.flowPotential()->setValues(finiteElementMethodPotentialFlow_.data().solution());
  finiteElementMethodPotentialFlow_.data().solution()->computeGradientField(data_.fibreDirection());

  advanceTimeStepping();
}

template<typename DiscretizableInTime>
void MultidomainSolver<FiniteElementMethodPotentialFlow,CellMLAdapter,FiniteElementMethodDiffusion>::
initialize()
{
  // TimeSteppingImplicit stores the diffusion finite element class
  TimeSteppingImplicit<FiniteElementMethodDiffusion>::initialize();

  // initialize potential flow
  finiteElementMethodPotentialFlow_.initialize();

  // initialize cellml adapter
  cellMLAdapter_.initialize();

  // initialize system matrix
  setSystemMatrix(this->timeStepWidth_);
}

template<typename DiscretizableInTime>
void MultidomainSolver<FiniteElementMethodPotentialFlow,CellMLAdapter,FiniteElementMethodDiffusion>::
setSystemMatrix(double timeStepWidth)
{

}


} // namespace TimeSteppingScheme
