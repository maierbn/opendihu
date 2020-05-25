#pragma once

#include <Python.h>  // has to be the first included header
#include <petscvec.h>
#include <functional>
#include <vector>

#include "interfaces/runnable.h"
#include "control/dihu_context.h"
#include "data_management/multiple_instances.h"
#include "output_writer/manager.h"
#include "partition/mesh_partition/02_mesh_partition.h"

namespace Control
{

/** This class holds multiple instances of the template type, e.g. for having multiple fibers, which are each as in example electrophysiology
  */
template<class TimeSteppingScheme>
class MultipleInstances: public Runnable
{
public:

  typedef std::vector<std::shared_ptr<typename TimeSteppingScheme::OutputConnectorDataType>> OutputConnectorDataType;
  typedef typename TimeSteppingScheme::FunctionSpace FunctionSpace;
  typedef typename ::Data::MultipleInstances<typename TimeSteppingScheme::FunctionSpace, TimeSteppingScheme> Data;
  typedef TimeSteppingScheme TimeSteppingSchemeType;

  //! constructor
  MultipleInstances(DihuContext context);

  //! advance simulation by the given time span [startTime_, endTime_]
  void advanceTimeSpan();

  //! set a new time interval that will be simulated by next call to advanceTimeSpan. This also potentially changes the time step width (it preserves the number of timesteps in the new time span)
  void setTimeSpan(double startTime, double endTime);

  //! initialize time span from specificSettings_
  void initialize();

  //! return the data object
  Data &data();

  //! get the data that will be transferred in the operator splitting to the other term of the splitting
  //! the transfer is done by the output_connector_data_transfer class
  std::shared_ptr<OutputConnectorDataType> getOutputConnectorData();

  //! run solution process
  void run();

  //! reset the objects state
  void reset();

  //! output the given data for debugging
  std::string getString(std::shared_ptr<OutputConnectorDataType> data);

  //! the FastMonodomainSolver accesses the internals of MultipleInstances
  std::vector<TimeSteppingScheme> &instancesLocal();
 
  /** write data to file using the output writer manager
   * @param callCountIncrement increment to the counter to reduce the output frequency in settings, this is only needed for FastMonodomainSolver
   */
  void writeOutput(int timeStepNo, double currentTime, int callCountIncrement = 1);

protected:

  DihuContext context_; //< the context object that holds the config for this class
  PythonConfig specificSettings_;    //< config for this object
  OutputWriter::Manager outputWriterManager_; //< manager object holding all output writer

  int nInstances_; //< number of instances that are given by config
  int nInstancesComputedGlobally_; //< number of instances that any process will compute
  std::vector<TimeSteppingScheme> instancesLocal_;   //< the instances of the problem that are computed on the local rank
  int nInstancesLocal_;   //< the number of local instances, i.e. the size of the instancesLocal_ vector

  std::shared_ptr<Partition::RankSubset> rankSubsetAllComputedInstances_;   //< the rank nos of all computed instances of this MultipleInstances object
  std::string logKey_;   //< the key under which the duration of all instances together is saved in the log

  Data data_;  //< the data object
  std::shared_ptr<OutputConnectorDataType> outputConnectorData_;       //< the outputConnectorData as vector of references o f teh outputConnectorData's of the instances

  bool outputInitializeThisInstance_;   //< if this instance displays progress of initialization
};

extern bool outputInitialize_;   //< if the message about initialization was already printed

}  // namespace

#include "control/multiple_instances.tpp"
