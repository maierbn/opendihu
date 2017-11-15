#include "time_stepping_scheme/explicit_euler.h"

#include <Python.h>

#include "control/python_utility.h"

namespace TimeSteppingScheme
{

template<typename DiscretizableInTime>
ExplicitEuler<DiscretizableInTime>::ExplicitEuler(DihuContext &context) : 
  context_(context), data_(context_), discretizableInTime(context_)
{
  PyObject *topLevelSettings = context_.getPythonConfig();
  specificSettings_ = PythonUtility::extractDict(topLevelSettings, "ExplicitEuler");
}

template<typename DiscretizableInTime>
void ExplicitEuler<DiscretizableInTime>::setInitialValues()
{
  int nDegreesOfFreedom = data_.nDegreesOfFreedom();
  Vec &solution = data_.solution();
  
  if (!discretizableInTime.setInitialValues(solution))
  {
   // initialize with 0
   PetscErrorCode ierr;
   ierr = VecSet(solution, 0.0); CHKERRV(ierr);
   
   // set from settings
   std::vector<double> values;
   PythonUtility::getOptionVector(specificSettings_, "initialValues", nDegreesOfFreedom, values);
   
   // loop over initialValues list
   for(int i=0; i<nDegreesOfFreedom; i++)
   {
     //                 vector    row  value
     ierr = VecSetValue(solution, i, values[i], INSERT_VALUES); CHKERRV(ierr);
   }
  }
}

template<typename DiscretizableInTime>
void ExplicitEuler<DiscretizableInTime>::run()
{
  // compute time stepping width
  endTime_ = PythonUtility::getOptionDouble(specificSettings_, "endTime", 1.0, PythonUtility::Positive);
  numberTimeSteps_ = PythonUtility::getOptionInt(specificSettings_, "numberTimeSteps", 10, PythonUtility::Positive);
  
  timeStepWidth_ = endTime_ / numberTimeSteps_;
 
  // initialize finite element part
  discretizableInTime.initialize();
  data_.setNDegreesOfFreedomPerNode(discretizableInTime.numberDegreesOfFreedomPerNode());
  data_.setMesh(discretizableInTime.mesh());
  
  // set initial values from settings
  setInitialValues();
  data_.print();
  
  // loop over time steps
  double currentTime = 0.0;
  for(int timeStepNo = 0; timeStepNo < numberTimeSteps_; timeStepNo++)
  {
    currentTime = double(timeStepNo) / (numberTimeSteps_-1) * endTime_;
    
    if (timeStepNo % 100 == 0)
     LOG(INFO) << "Timestep "<<timeStepNo<<"/"<<numberTimeSteps_<<", t="<<currentTime;
    
    // compute next delta_u = f(u)
    discretizableInTime.evaluateTimesteppingRightHandSide(data_.solution(), data_.increment(), timeStepNo, currentTime);
    
    // integrate, y += dt * delta_u
    VecAXPY(data_.solution(), timeStepWidth_, data_.increment());
    
    // write current output values
    context_.writeOutput(data_, timeStepNo, currentTime);
    
    data_.print();
  }
}

} // namespace TimeSteppingScheme