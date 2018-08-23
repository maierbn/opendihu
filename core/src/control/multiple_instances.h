#pragma once

#include <Python.h>  // has to be the first included header
#include <petscvec.h>

#include "control/runnable.h"
#include "control/dihu_context.h"
#include "data_management/multiple_instances.h"
#include "output_writer/manager.h"
#include "partition/01_mesh_partition.h"

namespace Control
{

/** This class holds multiple instances of the template type, e.g. for having multiple fibres, which are each as in example electrophysiology
  */
template<class TimeSteppingScheme>
class MultipleInstances: public Runnable
{
public:

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

  //! returns the Petsc solution vector
  Vec &solution();

  //! run solution process
  void run();

protected:

  DihuContext context_; ///< the context object that holds the config for this class
  PyObject *specificSettings_;    ///< config for this object
  OutputWriter::Manager outputWriterManager_; ///< manager object holding all output writer

  int nInstances_; ///< number of instances
  std::vector<TimeSteppingScheme> instancesLocal_;   ///< the instances of the problem that are computed on the local rank
  int nInstancesLocal_;   ///< the number of local instances, i.e. the size of the instancesLocal_ vector
  
  Data::MultipleInstances<typename TimeSteppingScheme::BasisOnMesh, TimeSteppingScheme> data_;  ///< the data object
};

};

#include "control/multiple_instances.tpp"