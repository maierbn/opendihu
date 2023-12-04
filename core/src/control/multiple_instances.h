#pragma once

#include <Python.h>  // has to be the first included header
#include <petscvec.h>
#include <functional>
#include <vector>

#include "interfaces/runnable.h"
#include "interfaces/multipliable.h"
#include "control/dihu_context.h"
#include "data_management/control/multiple_instances.h"
#include "output_writer/manager.h"
#include "partition/mesh_partition/02_mesh_partition.h"

namespace Control
{

/** This class holds multiple instances of the template type, e.g. for having multiple fibers, which are each as in example electrophysiology
  */
template<typename TimeSteppingScheme>
class MultipleInstances: public Runnable, public Multipliable
{
public:

  typedef std::vector<std::shared_ptr<typename TimeSteppingScheme::SlotConnectorDataType>> SlotConnectorDataType;
  typedef typename TimeSteppingScheme::FunctionSpace FunctionSpace;
  typedef typename ::Data::MultipleInstances<typename TimeSteppingScheme::FunctionSpace, TimeSteppingScheme> Data;
  typedef TimeSteppingScheme TimeSteppingSchemeType;

  //! constructor
  MultipleInstances(DihuContext context);

  //! advance simulation by the given time span [startTime_, endTime_]
  void advanceTimeSpan(bool withOutputWritersEnabled = true);

  //! set a new time interval that will be simulated by next call to advanceTimeSpan. This also potentially changes the time step width (it preserves the number of timesteps in the new time span)
  void setTimeSpan(double startTime, double endTime);

  //! initialize time span from specificSettings_
  void initialize();

  //! return the data object
  Data &data();

  //! get the data that will be transferred in the operator splitting to the other term of the splitting
  //! the transfer is done by the slot_connector_data_transfer class
  std::shared_ptr<SlotConnectorDataType> getSlotConnectorData();

  //! run solution process
  void run();

  //! reset the objects state
  void reset();

  //! output the given data for debugging
  std::string getString(std::shared_ptr<SlotConnectorDataType> data);

  //! the FastMonodomainSolver accesses the internals of MultipleInstances
  std::vector<TimeSteppingScheme> &instancesLocal();

  //! time of simulation
  double startTime();

  //! time of simulation
  double endTime();

  //! number of time steps in simulation time
  int numberTimeSteps();

  //! call the own output writer on the data object, output files will contain currentTime, with callCountIncrement !=1 output timesteps can be skipped
  void writeOwnOutput(int timeStepNo, double currentTime, int callCountIncrement = 1);

  //! call the output writer on the data object and all nested solvers, output files will contain currentTime, with callCountIncrement !=1 output timesteps can be skipped
  void callOutputWriter(int timeStepNo, double currentTime, int callCountIncrement = 1);

  void saveFiberDataCheckpoint();

  void restoreFiberDataCheckpoint();


protected:

  DihuContext context_;                         //< the context object that holds the config for this class
  PythonConfig specificSettings_;               //< config for this object
  OutputWriter::Manager outputWriterManager_;   //< manager object holding all output write
  Data data_;                                   //< the data object

  std::shared_ptr<SlotConnectorDataType> slotConnectorData_;       //< the slotConnectorData as vector of references o f teh slotConnectorData's of the instances

  int nInstances_;                              //< number of instances that are given by config
  int nInstancesComputedGlobally_;              //< number of instances that any process will compute
  int nInstancesLocal_;                         //< the number of local instances, i.e. the size of the instancesLocal_ vector

  std::vector<TimeSteppingScheme> instancesLocal_;   //< the instances of the problem that are computed on the local rank
  std::vector<std::shared_ptr<Partition::RankSubset>> rankSubsetsLocal_;  //< the rankSubset corresponding to the instances in instancesLocal_

  std::shared_ptr<Partition::RankSubset> rankSubsetAllComputedInstances_;   //< the rank nos of all computed instances of this MultipleInstances object
  std::string logKey_;                          //< the key under which the duration of all instances together is saved in the log

  bool outputInitializeThisInstance_;           //< if this instance displays progress of initialization
};

extern bool outputInitialize_;                  //< if the message about initialization was already printed

}  // namespace

#include "control/multiple_instances.tpp"
