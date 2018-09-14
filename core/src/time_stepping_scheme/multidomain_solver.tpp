#include "time_stepping_scheme/multidomain_solver.h"

#include <Python.h>  // has to be the first included header

#include "utility/python_utility.h"
#include "utility/petsc_utility.h"

namespace TimeSteppingScheme
{

template<typename FiniteElementMethodPotentialFlow,typename CellMLAdapter,typename FiniteElementMethodDiffusion>
MultidomainSolver<FiniteElementMethodPotentialFlow,CellMLAdapter,FiniteElementMethodDiffusion>::
MultidomainSolver(DihuContext context) :
  TimeSteppingImplicit<FiniteElementMethodDiffusion>(context, "MultidomainSolver"),
  finiteElementMethodPotentialFlow_(context),
  cellMLAdapter_(context)
{
  PyObject *topLevelSettings = this->context_.getPythonConfig();
  this->specificSettings_ = PythonUtility::getOptionPyObject(topLevelSettings, "MultidomainSolver");
  this->outputWriterManager_.initialize(this->specificSettings_);
}

template<typename FiniteElementMethodPotentialFlow,typename CellMLAdapter,typename FiniteElementMethodDiffusion>
void MultidomainSolver<FiniteElementMethodPotentialFlow,CellMLAdapter,FiniteElementMethodDiffusion>::
advanceTimeSpan()
{
  // compute timestep width
  double timeSpan = this->endTime_ - this->startTime_;

  LOG(DEBUG) << "MultidomainSolver::advanceTimeSpan, timeSpan=" << timeSpan<< ", timeStepWidth=" << this->timeStepWidth_
    << " n steps: " << this->numberTimeSteps_;

  // loop over time steps
  double currentTime = this->startTime_;
  for(int timeStepNo = 0; timeStepNo < this->numberTimeSteps_;)
  {
    if (timeStepNo % this->timeStepOutputInterval_ == 0)
     LOG(INFO) << "Timestep " << timeStepNo << "/" << this->numberTimeSteps_<< ", t=" << currentTime;

    //LOG(DEBUG) << "solution before integration: " << PetscUtility::getStringVector(this->data_.solution()->valuesGlobal());

    // advance computed value
    // compute next delta_u = f(u)
    this->discretizableInTime_.evaluateTimesteppingRightHandSide(
      this->data_.solution()->valuesGlobal(), this->data_.increment()->valuesGlobal(), timeStepNo, currentTime);

    // integrate, y += dt * delta_u
    VecAXPY(this->data_.solution()->valuesGlobal(), timeStepWidth, this->data_.increment()->valuesGlobal());

    // advance simulation time
    timeStepNo++;
    currentTime = this->startTime_ + double(timeStepNo) / this->numberTimeSteps_ * timeSpan;

    //LOG(DEBUG) << "solution after integration: " << PetscUtility::getStringVector(this->data_.solution()->valuesGlobal());
    // write current output values
    this->outputWriterManager_.writeOutput(this->data_, timeStepNo, currentTime);

    //this->data_.print();
  }
}

template<typename FiniteElementMethodPotentialFlow,typename CellMLAdapter,typename FiniteElementMethodDiffusion>
void MultidomainSolver<FiniteElementMethodPotentialFlow,CellMLAdapter,FiniteElementMethodDiffusion>::
run()
{
  TimeSteppingImplicit<FiniteElementMethodDiffusion>::run();
}

template<typename DiscretizableInTime>
void MultidomainSolver<FiniteElementMethodPotentialFlow,CellMLAdapter,FiniteElementMethodDiffusion>::
initialize()
{
  finiteElementMethodPotentialFlow_.initialize();

  cellMLAdapter_.initialize();

  // TimeSteppingImplicit stores the diffusion finite element class
  TimeSteppingImplicit<FiniteElementMethodDiffusion>::initialize();

}

} // namespace TimeSteppingScheme
