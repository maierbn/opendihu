#include "specialized_solver/quasi_static_nonlinear_elasticity_solver_chaste.h"

#include <Python.h>  // has to be the first included header

#include "utility/python_utility.h"
#include "utility/petsc_utility.h"
#include "data_management/specialized_solver/multidomain.h"
#include "control/performance_measurement.h"

namespace TimeSteppingScheme
{

template<int D>
QuasiStaticNonlinearElasticitySolverChaste<D>::
QuasiStaticNonlinearElasticitySolverChaste(DihuContext context) :
  context_(context["QuasiStaticNonlinearElasticitySolverChaste"]), data_(context_),
  initialized_(false)
{
  // get python config
  this->specificSettings_ = this->context_.getPythonConfig();

  // read in the durationLogKey
  if (specificSettings_.hasKey("durationLogKey"))
  {
    this->durationLogKey_ = specificSettings_.getOptionString("durationLogKey", "");
  }

  maximumActiveStress_ = specificSettings_.getOptionDouble("maximumActiveStress", 1.0, PythonUtility::ValidityCriterion::Positive);
  strainScalingCurveWidth_ = specificSettings_.getOptionDouble("strainScalingCurveWidth", 1.0, PythonUtility::ValidityCriterion::Positive);

  LOG(DEBUG) << "QuasiStaticNonlinearElasticitySolverChaste: parsed parameters maximumActiveStress: " << maximumActiveStress_ << ", strainScalingCurveWidth: " << strainScalingCurveWidth_;
  LOG(DEBUG) << "now parse output writers";

  // initialize output writers
  this->outputWriterManager_.initialize(this->context_, this->specificSettings_);
}

template<int D>
void QuasiStaticNonlinearElasticitySolverChaste<D>::
advanceTimeSpan()
{
  // start duration measurement, the name of the output variable can be set by "durationLogKey" in the config
  if (this->durationLogKey_ != "")
    Control::PerformanceMeasurement::start(this->durationLogKey_);

  LOG(TRACE) << "advanceTimeSpan, endTime: " << endTime_;

  // compute active stress from activation

  LOG(DEBUG) << "solve linear elasticity";

  // solve linear elasticity

  // stop duration measurement
  if (this->durationLogKey_ != "")
    Control::PerformanceMeasurement::stop(this->durationLogKey_);

  // write current output values
  this->outputWriterManager_.writeOutput(this->data_, 0, endTime_);
}

template<int D>
void QuasiStaticNonlinearElasticitySolverChaste<D>::
run()
{
  // initialize everything
  initialize();

  this->advanceTimeSpan();
}

template<int D>
void QuasiStaticNonlinearElasticitySolverChaste<D>::
setTimeSpan(double startTime, double endTime)
{
  endTime_ = endTime;
}

template<int D>
void QuasiStaticNonlinearElasticitySolverChaste<D>::
initialize()
{
  if (this->initialized_)
    return;

  LOG(DEBUG) << "initialize QuasiStaticNonlinearElasticitySolverChaste";
  assert(this->specificSettings_.pyObject());

  // create function space / mesh, the geometry is from the settings
  LOG(DEBUG) << "FiniteElementMethodBase constructor, create new function space from settings";
  std::shared_ptr<FunctionSpaceType> functionSpace = context_.meshManager()->functionSpace<FunctionSpaceType>(specificSettings_);

  // initialize the data object
  // store mesh in data
  data_.setFunctionSpace(functionSpace);

  data_.initialize();

  LOG(DEBUG) << "initialization done";
  this->initialized_ = true;
}

template<int D>
void QuasiStaticNonlinearElasticitySolverChaste<D>::reset()
{
  this->initialized_ = false;
}

//! return whether the underlying discretizableInTime object has a specified mesh type and is not independent of the mesh type
template<int D>
bool QuasiStaticNonlinearElasticitySolverChaste<D>::
knowsMeshType()
{
  return true;
}

template<int D>
typename QuasiStaticNonlinearElasticitySolverChaste<D>::Data &QuasiStaticNonlinearElasticitySolverChaste<D>::
data()
{
  return data_;
}

//! get the data that will be transferred in the operator splitting to the other term of the splitting
//! the transfer is done by the solution_vector_mapping class
template<int D>
typename QuasiStaticNonlinearElasticitySolverChaste<D>::TransferableSolutionDataType
QuasiStaticNonlinearElasticitySolverChaste<D>::
getSolutionForTransfer()
{
  return this->data_.activation();
}

//! output the given data for debugging
template<int D>
std::string QuasiStaticNonlinearElasticitySolverChaste<D>::
getString(typename QuasiStaticNonlinearElasticitySolverChaste<D>::TransferableSolutionDataType &data)
{
  std::stringstream s;
  s << "<QuasiStaticNonlinearElasticitySolverChaste:" << *data << ">";
  return s.str();
}

} // namespace TimeSteppingScheme
