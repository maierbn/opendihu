#include "time_stepping_scheme/time_stepping_scheme_ode.h"

#include <vector>

#include "utility/python_utility.h"

namespace TimeSteppingScheme
{

template<typename DiscretizableInTimeType>
TimeSteppingSchemeOde<DiscretizableInTimeType>::
TimeSteppingSchemeOde(const DihuContext &context, const std::string name) : 
  TimeSteppingScheme(context), data_(context), discretizableInTime_(context[name])
{
}

template<typename DiscretizableInTimeType>
Data::TimeStepping<typename DiscretizableInTimeType::BasisOnMesh> &TimeSteppingSchemeOde<DiscretizableInTimeType>::
data()
{
  return data_;
}

template<typename DiscretizableInTimeType>
void TimeSteppingSchemeOde<DiscretizableInTimeType>::
setInitialValues()
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

template<typename DiscretizableInTimeType>
SolutionVectorMapping &TimeSteppingSchemeOde<DiscretizableInTimeType>::
solutionVectorMapping()
{
  return discretizableInTime_.solutionVectorMapping();
}

template<typename DiscretizableInTimeType>
Vec &TimeSteppingSchemeOde<DiscretizableInTimeType>::
solution()
{
  return data_.solution().values();
}

template<typename DiscretizableInTimeType>
void TimeSteppingSchemeOde<DiscretizableInTimeType>::
initialize()
{
  TimeSteppingScheme::initialize();
  LOG(TRACE) << "TimeSteppingSchemeOde::initialize";
  
  // initialize underlying DiscretizableInTime object
  discretizableInTime_.initialize();
  data_.setNComponentsPerNode(discretizableInTime_.nComponentsNode());
  std::shared_ptr<Mesh::Mesh> mesh = discretizableInTime_.mesh();
  data_.setMesh(std::static_pointer_cast<typename DiscretizableInTimeType::BasisOnMesh>(mesh));
  data_.initialize();
  
  timeStepOutputInterval_ = PythonUtility::getOptionInt(specificSettings_, "timeStepOutputInterval", 100, PythonUtility::Positive);
  
  // set initial values from settings
  
  Vec &solution = data_.solution().values();
  
  if (!discretizableInTime_.setInitialValues(solution))
  {
    this->setInitialValues();
  }
  data_.print();
}

template<typename DiscretizableInTimeType>
void TimeSteppingSchemeOde<DiscretizableInTimeType>::
run()
{
  // initialize
  this->initialize();
  
  // do simulations
  this->advanceTimeSpan();
}

template<typename DiscretizableInTimeType>
bool TimeSteppingSchemeOde<DiscretizableInTimeType>::
knowsMeshType()
{
  return this->discretizableInTime_.knowsMeshType();
}


} // namespace