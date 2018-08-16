#include "time_stepping_scheme/implicit_euler.h"

#include <Python.h>  // has to be the first included header

#include "utility/python_utility.h"
#include "utility/petsc_utility.h"

namespace TimeSteppingScheme
{

template<typename DiscretizableInTime>
ImplicitEuler<DiscretizableInTime>::ImplicitEuler(DihuContext context) :
  TimeSteppingSchemeOde<DiscretizableInTime>(context, "ImplicitEuler")
{
  this->data_ = std::make_shared <Data::TimeStepping<typename DiscretizableInTime::BasisOnMesh, DiscretizableInTime::nComponents()>>(context); // create data object for implicit euler
  PyObject *topLevelSettings = this->context_.getPythonConfig();
  this->specificSettings_ = PythonUtility::getOptionPyObject(topLevelSettings, "ImplicitEuler");
  this->outputWriterManager_.initialize(this->specificSettings_);
}


template<typename DiscretizableInTime>
void ImplicitEuler<DiscretizableInTime>::initialize()
{
  if (initialized_)
    return;
  
  LOG(TRACE)<<"ImplicitEuler::initialize";
  
  TimeSteppingSchemeOde<DiscretizableInTime>::initialize();
  //In case required to initialize objects related to the implicit time stepping
  
  initialized_ = true;
  
}


template<typename DiscretizableInTime>
void ImplicitEuler<DiscretizableInTime>::advanceTimeSpan()
{
  // compute timestep width
  double timeSpan = this->endTime_ - this->startTime_;
  double timeStepWidth = timeSpan / this->numberTimeSteps_;

  LOG(DEBUG) << "ImplicitEuler::advanceTimeSpan, timeSpan="<<timeSpan<<", timeStepWidth="<<timeStepWidth
    <<" n steps: "<<this->numberTimeSteps_;
  
  int nEntries;
  VecGetSize(this->data_->solution().values(), &nEntries);

  // loop over time steps
  double currentTime = this->startTime_;
  
  
  for(int timeStepNo = 0; timeStepNo < this->numberTimeSteps_;)
  {
    if (timeStepNo % this->timeStepOutputInterval_ == 0)
     LOG(INFO) << "Timestep "<<timeStepNo<<"/"<<this->numberTimeSteps_<<", t="<<currentTime;
    
    /*
    PetscErrorCode ierr;
    double val_get;
    for (int i=0;i<nEntries;i++)
    {
      ierr=VecGetValues(this->data_->solution().values(),1,&i,&val_get);
      LOG(INFO)<<"val_get solution before solve " << i <<": " << val_get;   
    }
    */

    //prepare rhs for the variant 1 of the implicit Euler
    //this->discretizableInTime_.evaluateTimesteppingRightHandSideImplicit(this->data_->solution().values(),this->data_->rhs().values(),timeStepNo, currentTime);
    
    /*
    ierr=VecCopy(this->data_->solution().values(), this->data_->rhs().values());
    // initialize with 0
    ierr = VecSet(this->data_->solution().values(), 0.0); CHKERRV(ierr);
    */
    
    /*
    for (int i=0;i<nEntries;i++)
    {
      ierr=VecGetValues(this->data_->rhs().values(),1,&i,&val_get);
      LOG(INFO)<<"val_get:rhs "<< val_get;   
    } 
    */
    
    // computed value
    this->discretizableInTime_.solveLinearSystem(this->data_->solution().values(), this->data_->solution().values());
    //Variant 1 of the implicit Euler
    //this->discretizableInTime_.solveLinearSystem(this->data_->rhs().values(), this->data_->solution().values());
    
    /*
    for (int i=0;i<nEntries;i++)
    {
      ierr=VecGetValues(this->data_->solution().values(),1,&i,&val_get);
      LOG(INFO)<<"val_get solution after solve " << i <<": " << val_get;   
    }
    */

    // advance simulation time
    timeStepNo++;
    currentTime = this->startTime_ + double(timeStepNo) / this->numberTimeSteps_ * timeSpan;

    //LOG(DEBUG) << "solution after integration: " << PetscUtility::getStringVector(this->data_->solution().values());
    // write current output values
    this->outputWriterManager_.writeOutput(*this->data_, timeStepNo, currentTime);

    //this->data_->print();
  }
}

template<typename DiscretizableInTime>
void ImplicitEuler<DiscretizableInTime>::run()
{
  TimeSteppingSchemeOde<DiscretizableInTime>::run();
}
} // namespace TimeSteppingScheme