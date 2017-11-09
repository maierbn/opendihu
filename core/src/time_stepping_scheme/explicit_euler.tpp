#include "time_stepping_scheme/explicit_euler.h"

#include <Python.h>

#include "control/python_utility.h"

namespace TimeSteppingScheme
{

template<typename DiscretizableInTime>
ExplicitEuler<DiscretizableInTime>::ExplicitEuler(DihuContext &context) : 
  context_(context), data_(context_), discretizableInTime(context_)
{
}
  
template<typename DiscretizableInTime>
void ExplicitEuler<DiscretizableInTime>::run()
{
  PyObject *topLevelSettings = context_.getPythonConfig();
  PyObject *specificSettings = PythonUtility::extractDict(topLevelSettings, "ExplicitEuler");
  
  run(specificSettings);
}

template<typename DiscretizableInTime>
void ExplicitEuler<DiscretizableInTime>::setInitialValues(PyObject *specificSettings)
{
  int nDegreesOfFreedom = data_.mesh()->nDegreesOfFreedom();
  
  Vec &solution = data_.solution();
  
  // initialize with 0
  PetscErrorCode ierr;
  ierr = VecSet(solution, 0.0); CHKERRV(ierr);
  
  // set from settings
  // loop over initialValues list
 
  // get the first dirichlet boundary condition from the list
  double initialValue = PythonUtility::getOptionListBegin<double>(specificSettings, "initialValues");
  node_idx_t nodeNo = 0; 
 
  // loop over initial values and store them in rhs
  for (;
       !PythonUtility::getOptionListEnd(specificSettings, "initialValues")
       && nodeNo < nDegreesOfFreedom; 
       PythonUtility::getOptionListNext<double>(specificSettings, "initialValues", initialValue), nodeNo++)
  {
    //                 vector    row     value
    ierr = VecSetValue(solution, nodeNo, initialValue, INSERT_VALUES); CHKERRV(ierr);
  }
}

template<typename DiscretizableInTime>
void ExplicitEuler<DiscretizableInTime>::run(PyObject *specificSettings)
{
  endTime_ = PythonUtility::getOptionDouble(specificSettings, "endTime", 1.0, PythonUtility::Positive);
  numberTimeSteps_ = PythonUtility::getOptionInt(specificSettings, "numberTimeSteps_", 10, PythonUtility::Positive);
  
  timeStepWidth_ = endTime_ / numberTimeSteps_;
  
  discretizableInTime.initialize();
  data_.setMesh(discretizableInTime.mesh());
  
  setInitialValues(specificSettings);
  
  double currentTime = 0.0;
  
  for(int timeStepNo = 0; timeStepNo < numberTimeSteps_; timeStepNo++)
  {
    currentTime = timeStepNo / (numberTimeSteps_-1) * endTime_;
    
    // compute next delta_u = f(u)
    discretizableInTime.evaluateTimesteppingRightHandSide(data_.solution(), data_.increment());
    
    // integrate, y += dt * delta_u
    VecAXPY(data_.solution(), timeStepWidth_, data_.increment());
    
    // write current output values
    context_.writeOutput(data_, specificSettings, timeStepNo, currentTime);
  }
}

} // namespace TimeSteppingScheme