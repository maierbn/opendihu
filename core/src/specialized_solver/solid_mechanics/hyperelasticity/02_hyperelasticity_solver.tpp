#include "specialized_solver/solid_mechanics/hyperelasticity/02_hyperelasticity_solver.h"

#include <Python.h>  // has to be the first included header

#include "utility/python_utility.h"
#include "utility/petsc_utility.h"
#include "solver/solver_manager.h"
#include "data_management/specialized_solver/multidomain.h"
#include "control/diagnostic_tool/performance_measurement.h"
#include "control/diagnostic_tool/solver_structure_visualizer.h"
#include "partition/mesh_partition/01_mesh_partition_structured.h"

namespace SpatialDiscretization
{

template<typename Term,bool withLargeOutput,typename MeshType,int nDisplacementComponents>
HyperelasticitySolver<Term,withLargeOutput,MeshType,nDisplacementComponents>::
HyperelasticitySolver(DihuContext context, std::string settingsKey) :
  HyperelasticityMaterialComputations<Term,withLargeOutput,MeshType,nDisplacementComponents>::HyperelasticityMaterialComputations(context, settingsKey)
{
}

template<typename Term,bool withLargeOutput,typename MeshType,int nDisplacementComponents>
void HyperelasticitySolver<Term,withLargeOutput,MeshType,nDisplacementComponents>::
advanceTimeSpan(bool withOutputWritersEnabled)
{
  LOG_SCOPE_FUNCTION;

  // start duration measurement, the name of the output variable can be set by "durationLogKey" in the config
  if (this->durationLogKey_ != "")
    Control::PerformanceMeasurement::start(this->durationLogKey_);

  LOG(TRACE) << "advanceTimeSpan, endTime: " << endTime_;

  // write reference output values but don't increment counter
  // if (withOutputWritersEnabled)
  // {
  //   this->outputWriterManager_.writeOutput(this->data_, 0, 0.0, 0);
  //   this->outputWriterManagerPressure_.writeOutput(this->pressureDataCopy_, 0, 0.0, 0);
  // }

  this->nonlinearSolve();
  postprocessSolution();

  // stop duration measurement
  if (this->durationLogKey_ != "")
    Control::PerformanceMeasurement::stop(this->durationLogKey_);

  // write current output values
  if (withOutputWritersEnabled)
  {
    this->outputWriterManager_.writeOutput(this->data_, 1, endTime_);
    this->outputWriterManagerPressure_.writeOutput(this->pressureDataCopy_, 1, endTime_);
  }
}

template<typename Term,bool withLargeOutput,typename MeshType,int nDisplacementComponents>
void HyperelasticitySolver<Term,withLargeOutput,MeshType,nDisplacementComponents>::
run()
{
  // initialize everything
  LOG(DEBUG) << "call initialize in run()";
  this->initialize();

  this->advanceTimeSpan();
}

template<typename Term,bool withLargeOutput,typename MeshType,int nDisplacementComponents>
void HyperelasticitySolver<Term,withLargeOutput,MeshType,nDisplacementComponents>::
setTimeSpan(double startTime, double endTime)
{
  endTime_ = endTime;
}

template<typename Term,bool withLargeOutput,typename MeshType,int nDisplacementComponents>
void HyperelasticitySolver<Term,withLargeOutput,MeshType,nDisplacementComponents>::
callOutputWriter(int timeStepNo, double currentTime, int callCountIncrement)
{
  // call the own output writer
  this->outputWriterManager_.writeOutput(this->data_, 1, endTime_);
  this->outputWriterManagerPressure_.writeOutput(this->pressureDataCopy_, 1, endTime_);
}

} // namespace SpatialDiscretization
