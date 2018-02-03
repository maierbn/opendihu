#include "time_stepping_scheme/time_stepping_scheme_ode.h"

#include <vector>

#include "utility/python_utility.h"

namespace TimeSteppingScheme
{

template<typename DiscretizableInTime>
TimeSteppingSchemeOde<DiscretizableInTime>::TimeSteppingSchemeOde(const DihuContext &context, const std::string name) : 
  TimeSteppingScheme(context), data_(context), discretizableInTime_(context[name])
{
}

template<typename DiscretizableInTime>
Data::TimeStepping<typename DiscretizableInTime::BasisOnMesh> &TimeSteppingSchemeOde<DiscretizableInTime>::data()
{
  return data_;
}

template<typename DiscretizableInTime>
void TimeSteppingSchemeOde<DiscretizableInTime>::setInitialValues()
{
  int nDegreesOfFreedom = data_.nDegreesOfFreedom();
  Vec &solution = data_.solution().values();
 
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
  return discretizableInTime_.solutionVectorMapping();
}

template<typename DiscretizableInTime>
Vec &TimeSteppingSchemeOde<DiscretizableInTime>::solution()
{
  return data_.solution().values();
}

template<typename DiscretizableInTime>
void TimeSteppingSchemeOde<DiscretizableInTime>::
initialize()
{
  TimeSteppingScheme::initialize();
  LOG(TRACE) << "TimeSteppingSchemeOde::initialize";
  
  // initialize underlying DiscretizableInTime object
  discretizableInTime_.initialize();
  data_.setNComponentsPerNode(discretizableInTime_.nComponentsNode());
  std::shared_ptr<Mesh::Mesh> mesh = discretizableInTime_.mesh();
  data_.setMesh(std::static_pointer_cast<typename DiscretizableInTime::BasisOnMesh>(mesh));
  
  timeStepOutputInterval_ = PythonUtility::getOptionInt(specificSettings_, "timeStepOutputInterval", 100, PythonUtility::Positive);
  
  // set initial values from settings
  
  Vec &solution = data_.solution().values();
  
  if (!discretizableInTime_.setInitialValues(solution))
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

template<typename DiscretizableInTime>
bool TimeSteppingSchemeOde<DiscretizableInTime>::
knowsMeshType()
{
  return this->discretizableInTime_.knowsMeshType();
}


} // namespace