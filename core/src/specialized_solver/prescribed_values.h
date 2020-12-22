#pragma once

#include <Python.h>  // has to be the first included header

#include "interfaces/runnable.h"
#include "time_stepping_scheme/00_time_stepping_scheme.h"
#include "data_management/specialized_solver/prescribed_values.h"   // adjust this include

/** A dummy solver that sets its field variables to prescribed values that are given by a callback function in the python settings.
  */
template<typename FunctionSpaceType, int nComponents1=1, int nComponents2=1>
class PrescribedValues :
  public Runnable,
  public TimeSteppingScheme::TimeSteppingScheme
{
public:
  //! make the FunctionSpace available
  typedef FunctionSpaceType FunctionSpace;

  //! define the type of the data object
  typedef ::Data::PrescribedValues<FunctionSpaceType,nComponents1,nComponents2> Data;

  //! Define the type of data that will be transferred between solvers when there is a coupling scheme.
  typedef typename Data::SlotConnectorDataType SlotConnectorDataType;

  //! constructor, gets the DihuContext object which contains all python settings
  PrescribedValues(DihuContext context);

  //! advance simulation by the given time span [startTime_, endTime_]
  void advanceTimeSpan(bool withOutputWritersEnabled = true);

  //! initialize time span from specificSettings_
  void initialize();

  //! run solution process, this first calls initialize() and then run()
  void run();

  //! reset state of this object, such that a new initialize() is necessary ("uninitialize")
  void reset();

  //! call the output writer on the data object, output files will contain currentTime, with callCountIncrement !=1 output timesteps can be skipped
  void callOutputWriter(int timeStepNo, double currentTime, int callCountIncrement = 1);

  //! return the data object of the timestepping scheme, with the call to this method the output writers get the data to create their output files
  Data &data();

  //! Get the data that will be transferred in the operator splitting or coupling to the other term of the splitting/coupling.
  //! the transfer is done by the slot_connector_data_transfer class
  std::shared_ptr<SlotConnectorDataType> getSlotConnectorData();

protected:

  //! call the callback function
  void callCallbacks(int timeStepNo, double currentTime);

  Data data_;                                   //< the data object that stores at least all field variables that should be output by output writers.
  OutputWriter::Manager outputWriterManager_;   //< manager object holding all output writers

  std::vector<PyObject *> callbackFunctions1_;  //< callback functions that set the values of field variable 1
  std::vector<PyObject *> callbackFunctions2_;  //< callback functions that set the values of field variable 2

  PyObject *pyCallbackAdditionalParameter_;     //< the last argument to the callback functions which is given by the "additionalArgument" setting
  PyObject *pyGlobalNaturalDofsList_;           //< python list of global dof nos
};

#include "specialized_solver/prescribed_values.tpp"
