#pragma once

#include <Python.h>  // has to be the first included header
#include <petscvec.h>

#include "interfaces/runnable.h"
#include "control/dihu_context.h"
#include "data_management/multiple_instances.h"
#include "output_writer/manager.h"
#include "partition/mesh_partition/01_mesh_partition.h"

namespace Control
{

/** This class holds multiple instances of the template type, e.g. for having multiple fibers, which are each as in example electrophysiology
  */
template<class TimeSteppingScheme>
class MultipleInstances: public Runnable
{
public:

  typedef std::vector<typename TimeSteppingScheme::TransferableSolutionDataType> TransferableSolutionDataType;
  typedef typename TimeSteppingScheme::FunctionSpace FunctionSpace;
  typedef typename TimeSteppingScheme::Data Data;

  //! constructor
  MultipleInstances(DihuContext context);

  //! advance simulation by the given time span [startTime_, endTime_]
  void advanceTimeSpan();

  //! set a new time interval that will be simulated by next call to advanceTimeSpan. This also potentially changes the time step width (it preserves the number of timesteps in the new time span)
  void setTimeSpan(double startTime, double endTime);

  //! initialize time span from specificSettings_
  void initialize();

  //! return whether the scheme has a specified mesh type and is not independent of the mesh type
  bool knowsMeshType();
/*
  //! returns the Petsc solution vector
  Vec &solution();
*/

  //! get the data that will be transferred in the operator splitting to the other term of the splitting
  //! the transfer is done by the solution_vector_mapping class
  TransferableSolutionDataType getSolutionForTransferInOperatorSplitting();

  //! run solution process
  void run();

  //! reset the objects state
  void reset();

protected:

  DihuContext context_; ///< the context object that holds the config for this class
  PythonConfig specificSettings_;    ///< config for this object
  OutputWriter::Manager outputWriterManager_; ///< manager object holding all output writer

  int nInstances_; ///< number of instances that are given by config
  int nInstancesComputedGlobally_; ///< number of instances that any process will compute
  std::vector<TimeSteppingScheme> instancesLocal_;   ///< the instances of the problem that are computed on the local rank
  int nInstancesLocal_;   ///< the number of local instances, i.e. the size of the instancesLocal_ vector
  
  std::shared_ptr<Partition::RankSubset> rankSubsetAllComputedInstances_;   ///< the rank nos of all computed instances of this MultipleInstances object

  ::Data::MultipleInstances<typename TimeSteppingScheme::FunctionSpace, TimeSteppingScheme> data_;  ///< the data object
};

};

#include "control/multiple_instances.tpp"
