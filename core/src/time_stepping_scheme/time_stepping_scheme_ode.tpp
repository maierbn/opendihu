#include "time_stepping_scheme/time_stepping_scheme_ode.h"

#include <vector>

#include "utility/python_utility.h"

namespace TimeSteppingScheme
{

template<typename DiscretizableInTimeType>
TimeSteppingSchemeOde<DiscretizableInTimeType>::
TimeSteppingSchemeOde(DihuContext context, const std::string name) : 
  TimeSteppingScheme(context), data_(context), discretizableInTime_(context[name])
{
}

template<typename DiscretizableInTimeType>
Data::TimeStepping<typename DiscretizableInTimeType::BasisOnMesh, DiscretizableInTimeType::nComponents()> &TimeSteppingSchemeOde<DiscretizableInTimeType>::
data()
{
  return data_;
}

template<typename DiscretizableInTimeType>
void TimeSteppingSchemeOde<DiscretizableInTimeType>::
setInitialValues()
{
  dof_no_t nUnknowns = data_.nUnknowns();
  Vec &solution = data_.solution().values();
 
  // initialize with 0
  PetscErrorCode ierr;
  ierr = VecSet(solution, 0.0); CHKERRV(ierr);
  
  // set from settings
  std::vector<double> values;
  PythonUtility::getOptionVector(this->specificSettings_, "initialValues", nUnknowns, values);
  
  PetscUtility::setVector(values, solution);
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