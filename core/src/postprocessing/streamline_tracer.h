#pragma once

#include <Python.h>  // has to be the first included header
#include <vector>

#include "time_stepping_scheme/discretizable_in_time.h"
#include "control/runnable.h"
#include "data_management/streamline_tracer.h"

namespace Postprocessing
{

/** A class that traces streamlines through a given solution field.
 *  DiscretizableInTimeType can be e.g. FiniteElements module.
 */
template<typename DiscretizableInTimeType>
class StreamlineTracer : public Runnable
{
public:
  //! constructor
  StreamlineTracer(DihuContext context);

  //! initialize
  void initialize();

  //! run tracing of stream lines
  void run();

protected:

  //! trace the streamlines starting from seed points
  void traceStreamlines();

  const DihuContext context_;    ///< object that contains the python config for the current context and the global singletons meshManager and solverManager
  DiscretizableInTimeType problem_;   ///< the DiscretizableInTime object that is managed by this class

  Data::StreamlineTracer<typename DiscretizableInTimeType::BasisOnMesh, typename DiscretizableInTimeType::Data> data_;    ///< the data object that holds the gradient field variable

  PyObject *specificSettings_;   ///< the specific python config for this module
  double lineStepWidth_;     ///< the line step width used for integrating the streamlines
  std::vector<Vec3> seedPositions_;  ///< the seed points from where the streamlines start

  OutputWriter::Manager outputWriterManager_; ///< manager object holding all output writer
};

};  // namespace

#include "postprocessing/streamline_tracer.tpp"