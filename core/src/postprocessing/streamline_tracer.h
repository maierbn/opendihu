#pragma once

#include <Python.h>  // has to be the first included header
#include <vector>

#include "interfaces/discretizable_in_time.h"
#include "interfaces/runnable.h"
#include "postprocessing/streamline_tracer_base.h"
#include "data_management/streamline_tracer.h"

namespace Postprocessing
{

/** A class that traces streamlines through a given solution field.
 *  DiscretizableInTimeType can be e.g. FiniteElements module.
 *  This class also does some postprocessing like resampling of the streamlines to create equidistant elements.
 *  The actual tracing is done in the base class.
 */
template<typename DiscretizableInTimeType>
class StreamlineTracer :
  public StreamlineTracerBase<typename DiscretizableInTimeType::FunctionSpace>,
  public Runnable
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
  
  //! drop streamlines that are smaller than  discardRelativeLength_*median streamline length and resample to match targetElementLength_
  void postprocessStreamlines(std::vector<std::vector<Vec3>> &nodePositions);

  const DihuContext context_;    ///< object that contains the python config for the current context and the global singletons meshManager and solverManager
  DiscretizableInTimeType problem_;   ///< the DiscretizableInTime object that is managed by this class

  Data::StreamlineTracer<typename DiscretizableInTimeType::FunctionSpace, typename DiscretizableInTimeType::Data> data_;    ///< the data object that holds the gradient field variable

  PythonConfig specificSettings_;   ///< the specific python config for this module
  std::vector<Vec3> seedPositions_;  ///< the seed points from where the streamlines start

  OutputWriter::Manager outputWriterManager_; ///< manager object holding all output writer

  double targetElementLength_;   ///< the final length of each element of the traced streamlines. After the streamlines were traced using the fine lineStepWidth_, it gets resampled with this width.
  double targetLength_;           ///< the final length of the longest streamline, 0 means disabled
  double discardRelativeLength_;   ///< a relative length (in [0,1]), at the end streamlines are dropped that are smaller than this relative length times the median fiber length
  std::string csvFilename_;      ///< a csv output filename to write the node positions of the streamlines to (after postprocessing)
  std::string csvFilenameBeforePostprocessing_;      ///< a csv output filename to write the node positions of the streamlines to (before postprocessing)
};

} // namespace

#include "postprocessing/streamline_tracer.tpp"
