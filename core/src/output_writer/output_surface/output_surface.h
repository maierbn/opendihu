#pragma once

#include <Python.h>  // has to be the first included header

#include "function_space/function_space.h"
#include "data_management/output_surface/output_surface.h"

namespace OutputWriter
{

/** A class that has own OutputWriters and writes the 2D surface of a 3D mesh.
 */
template<typename Solver>
class OutputSurface :
  public Runnable
{
public:
  typedef typename Solver::FunctionSpace FunctionSpace;
  typedef typename Solver::Data Data3D;
  typedef typename Solver::OutputConnectorDataType OutputConnectorDataType;

  //! constructor
  OutputSurface(DihuContext context);

  //! advance simulation by the given time span [startTime_, endTime_]
  void advanceTimeSpan();

  //! initialize time span from specificSettings_
  void initialize();

  //! return whether the scheme has a specified mesh type and is not independent of the mesh type
  bool knowsMeshType();

  //! run solution process
  void run();

  //! reset state
  void reset();

  //! set a new time interval
  void setTimeSpan(double startTime, double endTime);

  //! return the data object of the timestepping scheme
  Data3D &data();

  //! get the data that will be transferred in the operator splitting to the other term of the splitting
  //! the transfer is done by the solution_vector_mapping class
  OutputConnectorDataType getOutputConnectorData();

protected:

  DihuContext context_;    ///< object that contains the python config for the current context and the global singletons meshManager and solverManager
  Solver solver_;     ///< the contained solver object

  bool initialized_ = false;   ///< if this object is initialized
  Data::OutputSurface<Data3D> data_;   ///< data object
  bool ownRankInvolvedInOutput_;   ///< if the own rank should call the output writer, because surface meshes are output, it can be that the surface is only contained on a subset of ranks
  int timeStepNo_;     ///< time step no for output writer
  double currentTime_;   ///< current simulation time for output writer

  Manager outputWriterManager_; ///< manager object holding all output writers

};

}  // namespace OutputWriter

#include "output_writer/output_surface/output_surface.tpp"
