#include "specialized_solver/solid_mechanics/quasi_static/quasi_static_linear_elasticity_solver.h"

#include <Python.h>  // has to be the first included header

#include "utility/python_utility.h"
#include "utility/petsc_utility.h"
#include "data_management/specialized_solver/multidomain.h"
#include "control/diagnostic_tool/performance_measurement.h"

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

  // avoid that solver structure file is created, this should only be done after the whole simulation has finished
  DihuContext::solverStructureVisualizer()->disable();

  finiteElementMethodLinearElasticity_.run();

  // enable again
  DihuContext::solverStructureVisualizer()->enable();


  LOG(DEBUG) << "compute strain";

  finiteElementMethodLinearElasticity_.data().computeStrain(data_.strain());

  //data_.debug();

  LOG(DEBUG) << "update geometry";

  // add displacements to geometry
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

  this->data_.activeStress()->setRepresentationGlobal();
  this->data_.activeStress()->startGhostManipulation();

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
      // get list of local dofNos of the current element
      dof_no_t dofNo = this->data_.functionSpace()->getDofNo(elementNoLocal, dofIndex);

      // compute the active stress
      double directionalStress = maximumActiveStress_ * activationValues[dofIndex];

      // scale stress by strain to obtain a characteristic force-length relationship
#if 0
      // get the strain in fiber direction
      int componentNo = D*D - 1;   // get the epsilon_zz entry which is the last
      double strainZ = data_.strain()->getValue(componentNo, dofNo);
      double relativeSarcomereLength = strainZ;

      // scale the stress dependent on the strain which corresponds to the relative sarcomere length
      double scalingFactor = 1 - MathUtility::sqr(relativeSarcomereLength / strainScalingCurveWidth_);   // scaling function 1 - (z/w)^2, w is the width of the parabola
      scalingFactor = 1.0;  // currently disabled for debugging
      directionalStress *= scalingFactor;
#endif

      MathUtility::Matrix<D,D> stressTensor = directionalStress * anisotropyTensor_.value(elementNoLocal);

      /*
      std::array<double,D*D>{
        directionalStress, 0, 0,
        0, 0, 0,
        0, 0, 0
      };*/

      VLOG(1) << "initial stressTensor: " << stressTensor;

      // rotate stress tensor in fiber direction by change of basis
      MathUtility::rotateMatrix(stressTensor, fiberDirectionValues[dofIndex]);

#if 0
      VLOG(1) << " dof " << dofNo << ", strainZ: " << strainZ << ", relativeSarcomereLength: " << relativeSarcomereLength
        << ", activationValues: " << activationValues << ", directionalStress: " << directionalStress
        << ", fiber direction: " << fiberDirectionValues[dofIndex] << ", rotated stress: " << stressTensor;
#endif

      // assign the computed active stress to the field variable
      std::array<double,D*D> stressTensorValues = std::array<double,D*D>(stressTensor);
      this->data_.activeStress()->setValue(dofNo, stressTensorValues, INSERT_VALUES);
    }
  }

  this->data_.activeStress()->finishGhostManipulation();
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

  // add this solver to the solvers diagram
  DihuContext::solverStructureVisualizer()->addSolver("QuasiStaticLinearElasticitySolver");

  // indicate in solverStructureVisualizer that now a child solver will be initialized
  DihuContext::solverStructureVisualizer()->beginChild();

  // set pointer to active stress in linear elasticity class
  finiteElementMethodLinearElasticity_.data().setActiveStress(data_.activeStress());
  finiteElementMethodLinearElasticity_.data().setRightHandSideActive(data_.rightHandSideActive());

  // initialize the finite element method
  finiteElementMethodLinearElasticity_.initialize();

  data_.setData(std::make_shared<DataLinearElasticityType>(finiteElementMethodLinearElasticity_.data()));

  // indicate in solverStructureVisualizer that the child solver initialization is done
  DihuContext::solverStructureVisualizer()->endChild();

  // if fiberDirection was given in settings, use this vector directly
  if (specificSettings_.hasKey("fiberDirection"))
  {
    fiberDirection_ = specificSettings_.getOptionArray<double,3>("fiberDirection", Vec3{0,0,1});

    // set fiber direction in field variable data_.fiberDirection()
    std::vector<Vec3> fiberDirectionValues(data_.functionSpace()->nDofsLocalWithoutGhosts(), fiberDirection_);
    data_.fiberDirection()->setValuesWithoutGhosts(fiberDirectionValues);

    usePotentialFlowForFiberDirection_ = false;
  }
  else
  {
    // fiberDirection is not given in settings, use a potential flow simulation (Laplace equation) and get the fiber direction from the streamlines
    usePotentialFlowForFiberDirection_ = true;
    DihuContext::solverStructureVisualizer()->beginChild("PotentialFlow");

    // initialize the potential flow finite element method
    finiteElementMethodPotentialFlow_.initialize();

    // indicate in solverStructureVisualizer that the child solver initialization is done
    DihuContext::solverStructureVisualizer()->endChild();

    LOG(INFO) << "Run potential flow simulation for fiber directions.";

    // avoid that solver structure file is created, this should only be done after the whole simulation has finished
    DihuContext::solverStructureVisualizer()->disable();

    // solve potential flow Laplace problem
    finiteElementMethodPotentialFlow_.run();

    // enable again
    DihuContext::solverStructureVisualizer()->enable();

    // compute a gradient field from the solution of the potential flow
    data_.flowPotential()->setValues(*finiteElementMethodPotentialFlow_.data().solution());
    data_.flowPotential()->computeGradientField(data_.fiberDirection());
  }

  // set the outputConnectorData for the solverStructureVisualizer to appear in the solver diagram
  DihuContext::solverStructureVisualizer()->setOutputConnectorData(getOutputConnectorData());

  // parse anisotropy tensor
  MathUtility::Matrix<3,3,double> defaultValue(std::array<double,9>{1, 0, 0,    0, 0, 0,    0, 0, 0});
  this->anisotropyTensor_.initialize(this->specificSettings_, "anisotropyTensor", defaultValue, this->data_.functionSpace());

  LOG(DEBUG) << "initialization done";
  this->initialized_ = true;
}

template<typename FiniteElementMethod>
void QuasiStaticLinearElasticitySolver<FiniteElementMethod>::reset()
{
  this->initialized_ = false;
}

template<typename FiniteElementMethod>
typename QuasiStaticLinearElasticitySolver<FiniteElementMethod>::Data &QuasiStaticLinearElasticitySolver<FiniteElementMethod>::
data()
{
  return data_;
}

//! get the data that will be transferred in the operator splitting to the other term of the splitting
//! the transfer is done by the output_connector_data_transfer class
template<typename FiniteElementMethod>
std::shared_ptr<typename QuasiStaticLinearElasticitySolver<FiniteElementMethod>::OutputConnectorDataType>
QuasiStaticLinearElasticitySolver<FiniteElementMethod>::
getOutputConnectorData()
{
  return this->data_.getOutputConnectorData();
}

//! output the given data for debugging
template<typename FiniteElementMethod>
std::string QuasiStaticLinearElasticitySolver<FiniteElementMethod>::
getString(std::shared_ptr<typename QuasiStaticLinearElasticitySolver<FiniteElementMethod>::OutputConnectorDataType> data)
{
  std::stringstream s;
  s << "<QuasiStaticLinearElasticitySolver:" << *data << ">";
  return s.str();
}

} // namespace TimeSteppingScheme
