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
  typedef typename Solver::TransferableSolutionDataType TransferableSolutionDataType;

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

  //! return the data object of the timestepping scheme
  Data3D &data();

  //! get the data that will be transferred in the operator splitting to the other term of the splitting
  //! the transfer is done by the solution_vector_mapping class
  TransferableSolutionDataType getSolutionForTransfer();

protected:

  DihuContext context_;    ///< object that contains the python config for the current context and the global singletons meshManager and solverManager
  Solver solver_;     ///< the contained solver object

  Data::OutputSurface<Data3D> data_;   ///< data object

  Manager outputWriterManager_; ///< manager object holding all output writers

};

}  // namespace OutputWriter

#include "output_writer/output_surface/output_surface.tpp"
