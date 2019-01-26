#pragma once

#include <Python.h>  // has to be the first included header
#include <petscvec.h>

#include "interfaces/runnable.h"
#include "interfaces/splittable.h"
#include "interfaces/discretizable_in_time.h"
#include "control/dihu_context.h"
#include "control/python_config.h"

namespace Control
{

/** This class implements adaptive time stepping and load balancing between multiple fibers.
  */
template<class TimeSteppingScheme>
class LoadBalancing:
  public Runnable,
  public Splittable
{
public:

  typedef typename TimeSteppingScheme::FunctionSpace FunctionSpace;
  typedef typename TimeSteppingScheme::Data Data;
  typedef typename TimeSteppingScheme::TransferableSolutionDataType TransferableSolutionDataType;

  //! constructor
  LoadBalancing(DihuContext context);

  //! advance simulation by the given time span [startTime_, endTime_]
  void advanceTimeSpan();

  //! set a new time interval that will be simulated by next call to advanceTimeSpan. This also potentially changes the time step width (it preserves the number of timesteps in the new time span)
  void setTimeSpan(double startTime, double endTime);

  //! initialize time span from specificSettings_
  void initialize();

  //! return whether the scheme has a specified mesh type and is not independent of the mesh type
  bool knowsMeshType();

  //! run solution process
  void run();

  //! reset state
  void reset();

  //! return the data object of the timestepping scheme
  Data &data();

  /*
  //! returns the Petsc solution vector
  std::shared_ptr<typename Data::FieldVariableType> solution();
  */

  //! get the data that will be transferred in the operator splitting to the other term of the splitting
  //! the transfer is done by the solution_vector_mapping class
  TransferableSolutionDataType getSolutionForTransfer();

protected:

  DihuContext context_;           ///< the context object that holds the config for this class
  PythonConfig specificSettings_;    ///< the python dictionary under "LoadBalancing"

  TimeSteppingScheme timeSteppingScheme_;   ///< the underlying timestepping method that is controlled by this class, e.g. Heun
};

}  // namespace

#include "control/load_balancing/load_balancing.tpp"
