#include "specialized_solver/hyperelasticity/quasi_static_hyperelasticity_solver.h"

#include <Python.h>  // has to be the first included header

#include "utility/python_utility.h"
#include "utility/petsc_utility.h"
#include "data_management/specialized_solver/multidomain.h"
#include "control/performance_measurement.h"

namespace TimeSteppingScheme
{

QuasiStaticHyperelasticitySolver::
QuasiStaticHyperelasticitySolver(DihuContext context) :
  context_(context["QuasiStaticHyperelasticitySolver"]), data_(context_), initialized_(false)
{
  // get python config
  this->specificSettings_ = this->context_.getPythonConfig();

  // read in the durationLogKey
  if (specificSettings_.hasKey("durationLogKey"))
  {
    this->durationLogKey_ = specificSettings_.getOptionString("durationLogKey", "");
  }
  c0_ = specificSettings_.getOptionDouble("c0", 1.0, PythonUtility::ValidityCriterion::Positive);
  c1_ = specificSettings_.getOptionDouble("c1", 1.0, PythonUtility::ValidityCriterion::Positive);

  LOG(DEBUG) << "QuasiStaticHyperelasticitySolver: parsed parameters c0: " << c0_ << ", c1: " << c1_;
  LOG(DEBUG) << "now parse output writers";

  // initialize output writers
  this->outputWriterManager_.initialize(this->context_, this->specificSettings_);
}


void QuasiStaticHyperelasticitySolver::
advanceTimeSpan()
{
  // start duration measurement, the name of the output variable can be set by "durationLogKey" in the config
  if (this->durationLogKey_ != "")
    Control::PerformanceMeasurement::start(this->durationLogKey_);

  LOG(TRACE) << "advanceTimeSpan, endTime: " << endTime_;

  // stop duration measurement
  if (this->durationLogKey_ != "")
    Control::PerformanceMeasurement::stop(this->durationLogKey_);

  // write current output values
  this->outputWriterManager_.writeOutput(this->data_, 0, endTime_);
}

void QuasiStaticHyperelasticitySolver::
run()
{
  // initialize everything
  initialize();

  this->advanceTimeSpan();
}


void QuasiStaticHyperelasticitySolver::
setTimeSpan(double startTime, double endTime)
{
  endTime_ = endTime;
}


void QuasiStaticHyperelasticitySolver::
initialize()
{
  if (this->initialized_)
    return;

  LOG(DEBUG) << "initialize QuasiStaticHyperelasticitySolver";
  assert(this->specificSettings_.pyObject());

  // create function space / mesh, the geometry is from the settings
  displacementsFunctionSpace_ = context_.meshManager()->functionSpace<DisplacementsFunctionSpace>(specificSettings_);

  // create 3D function space with linear basis functions
  //FunctionSpaceDofsNodes(std::shared_ptr<Partition::Manager> partitionManager, const std::vector<Vec3> &nodePositions,
  //                       const std::array<element_no_t,D> nElementsPerCoordinateDirection, const std::array<int,D> nRanksPerCoordinateDirection);

  std::vector<Vec3> nodePositionsLinearMesh;

  // loop over nodes of quadratic function space and extract nodes for linear function space
  for (node_no_t nodeIndexZ = 0; nodeIndexZ < displacementsFunctionSpace_->meshPartition()->nNodesLocalWithoutGhosts(2); nodeIndexZ++)
  {
    for (node_no_t nodeIndexY = 0; nodeIndexY < displacementsFunctionSpace_->meshPartition()->nNodesLocalWithoutGhosts(1); nodeIndexY++)
    {
      for (node_no_t nodeIndexX = 0; nodeIndexX < displacementsFunctionSpace_->meshPartition()->nNodesLocalWithoutGhosts(0); nodeIndexX++)
      {
        node_no_t nodeNoLocal = nodeIndexZ*displacementsFunctionSpace_->meshPartition()->nNodesLocalWithoutGhosts(1)*displacementsFunctionSpace_->meshPartition()->nNodesLocalWithoutGhosts(0)
          + nodeIndexY*displacementsFunctionSpace_->meshPartition()->nNodesLocalWithoutGhosts(0) + nodeIndexX;

        dof_no_t dofNoLocal = nodeNoLocal;   // no Hermite, i.e. no multiple dofs per node

        if (nodeIndexX % 2 == 0 && nodeIndexY % 2 == 0 && nodeIndexZ % 2 == 0)
        {
          nodePositionsLinearMesh.push_back(displacementsFunctionSpace_->geometryField().getValue(dofNoLocal));
        }
      }
    }
  }

  std::array<element_no_t,3> nElementsPerCoordinateDirection;
  std::array<int,3> nRanksPerCoordinateDirection;

  for (int i = 0; i < 3; i++)
  {
    nElementsPerCoordinateDirection[i] = displacementsFunctionSpace_->meshPartition()->nElementsLocal(i);
    nRanksPerCoordinateDirection[i] = displacementsFunctionSpace_->meshPartition()->nRanks(i);
  }

  pressureFunctionSpace_ = context_.meshManager()->createFunctionSpace<PressureFunctionSpace>("pressureFunctionSpace",nodePositionsLinearMesh, nElementsPerCoordinateDirection, nRanksPerCoordinateDirection);

  // initialize the data object
  // store mesh in data
  data_.setDisplacementsFunctionSpace(displacementsFunctionSpace_);
  data_.setPressureFunctionSpace(pressureFunctionSpace_);

  data_.initialize();


  LOG(DEBUG) << "initialization done";
  this->initialized_ = true;
}


void QuasiStaticHyperelasticitySolver::reset()
{
  this->initialized_ = false;
}

//! return whether the underlying discretizableInTime object has a specified mesh type and is not independent of the mesh type

bool QuasiStaticHyperelasticitySolver::
knowsMeshType()
{
  return true;
}

typename QuasiStaticHyperelasticitySolver::Data &
QuasiStaticHyperelasticitySolver::
data()
{
  return data_;
}

} // namespace TimeSteppingScheme
