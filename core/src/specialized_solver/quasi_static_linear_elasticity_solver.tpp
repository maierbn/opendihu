#include "specialized_solver/static_bidomain_solver.h"

#include <Python.h>  // has to be the first included header

#include "utility/python_utility.h"
#include "utility/petsc_utility.h"
#include "data_management/specialized_solver/multidomain.h"
#include "control/performance_measurement.h"

namespace TimeSteppingScheme
{

template<typename FiniteElementMethod>
QuasiStaticLinearElasticitySolver<FiniteElementMethodPotentialFlow,FiniteElementMethodDiffusion>::
QuasiStaticLinearElasticitySolver(DihuContext context) :
  context_(context["QuasiStaticLinearElasticitySolver"]),
  finiteElementMethodLinearElasticity_(this->context_), initialized_(false)
{
  // get python config
  this->specificSettings_ = this->context_.getPythonConfig();

  // read in the durationLogKey
  if (specificSettings_.hasKey("durationLogKey"))
  {
    this->durationLogKey_ = specificSettings_.getOptionString("durationLogKey", "");
  }

  // initialize output writers
  this->outputWriterManager_.initialize(this->context_, this->specificSettings_);
}

template<typename FiniteElementMethod>
void QuasiStaticLinearElasticitySolver<FiniteElementMethodPotentialFlow,FiniteElementMethodDiffusion>::
advanceTimeSpan()
{
  // start duration measurement, the name of the output variable can be set by "durationLogKey" in the config
  if (this->durationLogKey_ != "")
    Control::PerformanceMeasurement::start(this->durationLogKey_);

  // compute active stress from activation

  //this->data_->activation();
  //this->data_->activeStress();

  // set active stress in linear elasticity class
  finiteElementMethodLinearElasticity_->data().setActiveStress(data_->activeStress());

  // solve linear elasticity
  finiteElementMethodLinearElasticity_->solve();

  // add displacements to geometry
  finiteElementMethodLinearElasticity_->data().updateGeometry();

  // stop duration measurement
  if (this->durationLogKey_ != "")
    Control::PerformanceMeasurement::stop(this->durationLogKey_);

  // write current output values
  this->outputWriterManager_.writeOutput(this->data_, 0, endTime_);
}

template<typename FiniteElementMethod>
void QuasiStaticLinearElasticitySolver<FiniteElementMethodPotentialFlow,FiniteElementMethodDiffusion>::
run()
{
  // initialize everything
  initialize();

  this->advanceTimeSpan();
}

template<typename FiniteElementMethod>
void QuasiStaticLinearElasticitySolver<FiniteElementMethodPotentialFlow,FiniteElementMethodDiffusion>::
setTimeSpan(double startTime, double endTime)
{
  endTime_ = endTime;
}

template<typename FiniteElementMethod>
void QuasiStaticLinearElasticitySolver<FiniteElementMethodPotentialFlow,FiniteElementMethodDiffusion>::
initialize()
{
  if (this->initialized_)
    return;

  LOG(DEBUG) << "initialize static_bidomain_solver";
  assert(this->specificSettings_.pyObject());

  // initialize the finite element method, this also creates the function space
  finiteElementMethodLinearElasticity_.initialize();

  // initialize the data object
  data_.setFunctionSpace(finiteElementMethodLinearElasticity_.functionSpace());
  data_.initialize();
  data_.setData(std::make_shared<FiniteElementMethod>(finiteElementMethodLinearElasticity_->data()));

  LOG(DEBUG) << "initialization done";
  this->initialized_ = true;
}

template<typename FiniteElementMethod>
void QuasiStaticLinearElasticitySolver<FiniteElementMethodPotentialFlow,FiniteElementMethodDiffusion>::reset()
{
  this->initialized_ = false;
}

//! return whether the underlying discretizableInTime object has a specified mesh type and is not independent of the mesh type
template<typename FiniteElementMethod>
bool QuasiStaticLinearElasticitySolver<FiniteElementMethodPotentialFlow,FiniteElementMethodDiffusion>::
knowsMeshType()
{
  return true;
}

template<typename FiniteElementMethod>
typename QuasiStaticLinearElasticitySolver<FiniteElementMethodPotentialFlow,FiniteElementMethodDiffusion>::Data &QuasiStaticLinearElasticitySolver<FiniteElementMethodPotentialFlow,FiniteElementMethodDiffusion>::
data()
{
  return data_;
}

//! get the data that will be transferred in the operator splitting to the other term of the splitting
//! the transfer is done by the solution_vector_mapping class
template<typename FiniteElementMethod>
typename QuasiStaticLinearElasticitySolver<FiniteElementMethodPotentialFlow,FiniteElementMethodDiffusion>::TransferableSolutionDataType
QuasiStaticLinearElasticitySolver<FiniteElementMethodPotentialFlow,FiniteElementMethodDiffusion>::
getSolutionForTransfer()
{
  return this->data_.activation();
}

//! output the given data for debugging
template<typename FiniteElementMethod>
std::string QuasiStaticLinearElasticitySolver<FiniteElementMethodPotentialFlow,FiniteElementMethodDiffusion>::
getString(typename QuasiStaticLinearElasticitySolver<FiniteElementMethodPotentialFlow,FiniteElementMethodDiffusion>::TransferableSolutionDataType &data)
{
  std::stringstream s;
  s << "<QuasiStaticLinearElasticitySolver:" << *data << ">";
  return s.str();
}

} // namespace TimeSteppingScheme
