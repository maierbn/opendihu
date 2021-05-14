#pragma once

#include <Python.h>  // has to be the first included header

#include "interfaces/runnable.h"
#include "control/precice/volume_coupling/01_read_write.h"

namespace Control
{

/** Generic Precice adapter for volume coupling.
 */
template<typename NestedSolver>
class PreciceAdapterVolumeCoupling :
  public PreciceAdapterVolumeCouplingReadWrite<NestedSolver>
{
public:

  //! define the type of the data object
  typedef typename NestedSolver::Data Data;
  typedef typename NestedSolver::SlotConnectorDataType SlotConnectorDataType;

  //! constructor
  using PreciceAdapterVolumeCouplingReadWrite<NestedSolver>::PreciceAdapterVolumeCouplingReadWrite;

  //! run solution process, this first calls initialize() and then run()
  void run();

  //! reset state of this object, such that a new initialize() is necessary ("uninitialize")
  void reset();

  //! return the data object, with the call to this method the output writers get the data to create their output files
  Data &data();

  //! get the data that will be transferred in the operator splitting to the other term of the splitting
  //! the transfer is done by the slot_connector_data_transfer class
  std::shared_ptr<SlotConnectorDataType> getSlotConnectorData();

  //! set a new time interval that will be simulated by next call to advanceTimeSpan. This also potentially changes the time step width (it preserves the number of timesteps in the new time span)
  void setTimeSpan(double startTime, double endTime);

  //! advance simulation by the given time span [startTime_, endTime_] with given numberTimeSteps
  void advanceTimeSpan(bool withOutputWritersEnabled=true);

  //! call the output writer on the data object, output files will contain currentTime, with callCountIncrement !=1 output timesteps can be skipped
  void callOutputWriter(int timeStepNo, double currentTime, int callCountIncrement = 1);

};

}  // namespace

#include "control/precice/volume_coupling/precice_adapter_volume_coupling.tpp"
