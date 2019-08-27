#include "specialized_solver/quasi_static_linear_elasticity_solver.h"

#include <Python.h>  // has to be the first included header

#include "utility/python_utility.h"
#include "utility/petsc_utility.h"
#include "data_management/specialized_solver/multidomain.h"
#include "control/performance_measurement.h"

namespace TimeSteppingScheme
{

template<typename FiniteElementMethod>
QuasiStaticLinearElasticitySolver<FiniteElementMethod>::
QuasiStaticLinearElasticitySolver(DihuContext context) :
  context_(context["QuasiStaticLinearElasticitySolver"]), data_(context_),
  finiteElementMethodPotentialFlow_(this->context_["PotentialFlow"]),
  finiteElementMethodLinearElasticity_(this->context_), initialized_(false)
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
  scalingFactor_ = specificSettings_.getOptionDouble("scalingFactor", 1.0);

  LOG(DEBUG) << "QuasiStaticLinearElasticitySolver: parsed parameters maximumActiveStress: " << maximumActiveStress_ << ", strainScalingCurveWidth: " << strainScalingCurveWidth_;
  LOG(DEBUG) << "now parse output writers";

  // initialize output writers
  this->outputWriterManager_.initialize(this->context_, this->specificSettings_);
}

template<typename FiniteElementMethod>
void QuasiStaticLinearElasticitySolver<FiniteElementMethod>::
advanceTimeSpan()
{
  // start duration measurement, the name of the output variable can be set by "durationLogKey" in the config
  if (this->durationLogKey_ != "")
    Control::PerformanceMeasurement::start(this->durationLogKey_);

  LOG(TRACE) << "advanceTimeSpan, endTime: " << endTime_;

  // compute active stress from activation
  computeActiveStress();   //this->data_.activation(); to this->data_.activeStress() = finiteElementMethodLinearElasticity_.data().activeStress()

  LOG(DEBUG) << "solve linear elasticity";

  // solve linear elasticity
  finiteElementMethodLinearElasticity_.setRightHandSide();
  finiteElementMethodLinearElasticity_.applyBoundaryConditions();
  finiteElementMethodLinearElasticity_.run();

  LOG(DEBUG) << "compute strain";

  // add displacements to geometry
  finiteElementMethodLinearElasticity_.data().computeStrain(data_.strain());

  LOG(DEBUG) << "update geometry";

  finiteElementMethodLinearElasticity_.data().updateGeometry(this->scalingFactor_);

  // stop duration measurement
  if (this->durationLogKey_ != "")
    Control::PerformanceMeasurement::stop(this->durationLogKey_);

  // write current output values
  this->outputWriterManager_.writeOutput(this->data_, 0, endTime_);
}

template<typename FiniteElementMethod>
void QuasiStaticLinearElasticitySolver<FiniteElementMethod>::
computeActiveStress()
{
  //this->data_.activation(); to this->data_.activeStress();, use data_.fiberDirection()
  const int nDofsPerElement = FunctionSpace::nDofsPerElement();
  const int D = 3;

  LOG(TRACE) << "";
  LOG(DEBUG) << "activation: " << *this->data_.activation();

  // loop over elements
  for (element_no_t elementNoLocal = 0; elementNoLocal < this->data_.functionSpace()->nElementsLocal(); elementNoLocal++)
  {
    std::array<double,nDofsPerElement> activationValues;
    this->data_.activation()->getElementValues(elementNoLocal, activationValues);

    std::array<Vec3,nDofsPerElement> fiberDirectionValues;
    this->data_.fiberDirection()->getElementValues(elementNoLocal, fiberDirectionValues);

    VLOG(1) << "element " << elementNoLocal << ", activationValues: " << activationValues;

    // loop over dofs of element
    for (int dofIndex = 0; dofIndex < nDofsPerElement; dofIndex++)
    {
      double directionalStress = maximumActiveStress_ * activationValues[dofIndex];

      dof_no_t dofNo = this->data_.functionSpace()->getDofNo(elementNoLocal, dofIndex);

      // get the strain in fiber direction
      int componentNo = D*D - 1;   // get the epsilon_zz entry which is the last
      double strainZ = data_.strain()->getValue(componentNo, dofNo);
      double relativeSarcomereLength = strainZ;

      // scale the stress dependent on the strain which corresponds to the relative sarcomere length
      double scalingFactor = 1 - MathUtility::sqr(relativeSarcomereLength / strainScalingCurveWidth_);   // scaling function 1 - (z/w)^2, w is the width of the parabola
      scalingFactor = 1.0;  // currently disabled for debugging
      directionalStress *= scalingFactor;

      MathUtility::Matrix<D,D> stressTensor = std::array<double,D*D>{
        directionalStress, 0, 0,
        0, 0, 0,
        0, 0, 0
      };

      // rotate stress tensor in fiber direction
      MathUtility::rotateMatrix(stressTensor, fiberDirectionValues[dofIndex]);

      LOG(DEBUG) << " dof " << dofNo << ", strainZ: " << strainZ << ", relativeSarcomereLength: " << relativeSarcomereLength
        << ", activationValues: " << activationValues << ", directionalStress: " << directionalStress
        << ", fiber direction: " << fiberDirectionValues[dofIndex] << ", rotated stress: " << stressTensor;

      std::array<double,D*D> stressTensorValues = std::array<double,D*D>(stressTensor);
      this->data_.activeStress()->setValue(dofNo, stressTensorValues);
    }
  }

}

template<typename FiniteElementMethod>
void QuasiStaticLinearElasticitySolver<FiniteElementMethod>::
run()
{
  // initialize everything
  initialize();

  this->advanceTimeSpan();
}

template<typename FiniteElementMethod>
void QuasiStaticLinearElasticitySolver<FiniteElementMethod>::
setTimeSpan(double startTime, double endTime)
{
  endTime_ = endTime;
}

template<typename FiniteElementMethod>
void QuasiStaticLinearElasticitySolver<FiniteElementMethod>::
initialize()
{
  if (this->initialized_)
    return;

  LOG(DEBUG) << "initialize QuasiStaticLinearElasticitySolver";
  assert(this->specificSettings_.pyObject());

  // initialize the data object
  data_.setFunctionSpace(finiteElementMethodLinearElasticity_.functionSpace());
  data_.initialize();

  // set pointer to active stress in linear elasticity class
  finiteElementMethodLinearElasticity_.data().setActiveStress(data_.activeStress());
  finiteElementMethodLinearElasticity_.data().setRightHandSideActive(data_.rightHandSideActive());

  // initialize the finite element method
  finiteElementMethodLinearElasticity_.initialize();

  data_.setData(std::make_shared<DataLinearElasticityType>(finiteElementMethodLinearElasticity_.data()));

  // initialize the potential flow finite element method
  finiteElementMethodPotentialFlow_.initialize();

  LOG(INFO) << "Run potential flow simulation for fiber directions.";

  // solve potential flow Laplace problem
  finiteElementMethodPotentialFlow_.run();

  // compute a gradient field from the solution of the potential flow
  data_.flowPotential()->setValues(*finiteElementMethodPotentialFlow_.data().solution());
  data_.flowPotential()->computeGradientField(data_.fiberDirection());

  LOG(DEBUG) << "initialization done";
  this->initialized_ = true;
}

template<typename FiniteElementMethod>
void QuasiStaticLinearElasticitySolver<FiniteElementMethod>::reset()
{
  this->initialized_ = false;
}

//! return whether the underlying discretizableInTime object has a specified mesh type and is not independent of the mesh type
template<typename FiniteElementMethod>
bool QuasiStaticLinearElasticitySolver<FiniteElementMethod>::
knowsMeshType()
{
  return true;
}

template<typename FiniteElementMethod>
typename QuasiStaticLinearElasticitySolver<FiniteElementMethod>::Data &QuasiStaticLinearElasticitySolver<FiniteElementMethod>::
data()
{
  return data_;
}

//! get the data that will be transferred in the operator splitting to the other term of the splitting
//! the transfer is done by the solution_vector_mapping class
template<typename FiniteElementMethod>
typename QuasiStaticLinearElasticitySolver<FiniteElementMethod>::OutputConnectorDataType
QuasiStaticLinearElasticitySolver<FiniteElementMethod>::
getOutputConnectorData()
{
  return this->data_.activation();
}

//! output the given data for debugging
template<typename FiniteElementMethod>
std::string QuasiStaticLinearElasticitySolver<FiniteElementMethod>::
getString(typename QuasiStaticLinearElasticitySolver<FiniteElementMethod>::OutputConnectorDataType &data)
{
  std::stringstream s;
  s << "<QuasiStaticLinearElasticitySolver:" << *data << ">";
  return s.str();
}

} // namespace TimeSteppingScheme
