#include "time_stepping_scheme/time_stepping_scheme_ode.h"

#include <vector>

#include "control/python_utility.h"

namespace TimeSteppingScheme
{

template<typename DiscretizableInTime>
TimeSteppingSchemeOde<DiscretizableInTime>::TimeSteppingSchemeOde(const DihuContext &context) : 
  TimeSteppingScheme(context), data_(context), discretizableInTime(context)
{
}

template<typename DiscretizableInTime>
Data::Data &TimeSteppingSchemeOde<DiscretizableInTime>::data()
{
  return data_;
}

template<typename DiscretizableInTime>
void TimeSteppingSchemeOde<DiscretizableInTime>::setInitialValues()
{
  int nDegreesOfFreedom = data_.nDegreesOfFreedom();
  Vec &solution = data_.solution();
 
  // initialize with 0
  PetscErrorCode ierr;
  ierr = VecSet(solution, 0.0); CHKERRV(ierr);
  
  // set from settings
  std::vector<double> values;
  PythonUtility::getOptionVector(this->specificSettings_, "initialValues", nDegreesOfFreedom, values);
  
  // loop over initialValues list
  for(int i=0; i<nDegreesOfFreedom; i++)
  {
    //                 vector    row  value
    ierr = VecSetValue(solution, i, values[i], INSERT_VALUES); CHKERRV(ierr);
  }
}

template<typename DiscretizableInTime>
SolutionVectorMapping &TimeSteppingSchemeOde<DiscretizableInTime>::solutionVectorMapping()
{
  return discretizableInTime.solutionVectorMapping();
}

template<typename DiscretizableInTime>
Vec &TimeSteppingSchemeOde<DiscretizableInTime>::solution()
{
  return data_.solution();
}

template<typename DiscretizableInTime>
void TimeSteppingSchemeOde<DiscretizableInTime>::
initialize()
{
  TimeSteppingScheme::initialize();
  
  // initialize underlying DiscretizableInTime object
  discretizableInTime.initialize();
  data_.setNDegreesOfFreedomPerNode(discretizableInTime.numberDegreesOfFreedomPerNode());
  data_.setMesh(discretizableInTime.mesh());
  
  timeStepOutputFrequency_ = PythonUtility::getOptionInt(specificSettings_, "timeStepOutputFrequency", 100, PythonUtility::Positive);
  
  // set initial values from settings
  
  Vec &solution = data_.solution();
  
  if (!discretizableInTime.setInitialValues(solution))
  {
    this->setInitialValues();
  }
  data_.print();
}

template<typename DiscretizableInTime>
void TimeSteppingSchemeOde<DiscretizableInTime>::
run()
{
  // initialize
  this->initialize();
  
  // do simulations
  this->advanceTimeSpan();
}

} // namespace