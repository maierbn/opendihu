#pragma once

#include <Python.h>  // has to be the first included header
#include <vector>

#include "discretizable_in_time/discretizable_in_time.h"
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

  //! trace the streamline starting from startingPoint in the element initialElementNo, direction is either 1. or -1. depending on the direction
  void traceStreamline(element_no_t initialElementNo, std::array<double,(unsigned long int)3> xi, Vec3 startingPoint, double direction, std::vector<Vec3> &points);
  
  //! trace the streamlines starting from seed points
  void traceStreamlines();
  
  //! drop streamlines that are smaller than  discardRelativeLength_*median streamline length and resample to match targetElementLength_
  void postprocessStreamlines(std::vector<std::vector<Vec3>> &nodePositions);

  const DihuContext context_;    ///< object that contains the python config for the current context and the global singletons meshManager and solverManager
  DiscretizableInTimeType problem_;   ///< the DiscretizableInTime object that is managed by this class

  Data::StreamlineTracer<typename DiscretizableInTimeType::FunctionSpace, typename DiscretizableInTimeType::Data> data_;    ///< the data object that holds the gradient field variable

  PyObject *specificSettings_;   ///< the specific python config for this module
  double lineStepWidth_;     ///< the line step width used for integrating the streamlines
  std::vector<Vec3> seedPositions_;  ///< the seed points from where the streamlines start

  OutputWriter::Manager outputWriterManager_; ///< manager object holding all output writer
  int maxNIterations_;   ///< the maximum number of iterations to trace for a streamline
  bool useGradientField_;  ///< There are 2 implementations of streamline tracing. The first one (useGradientField_) uses a precomputed gradient field that is interpolated linearly and the second uses the gradient directly from the Laplace solution field. // The first one seems more stable, because the gradient is zero and the position of the boundary conditions.

  double targetElementLength_;   ///< the final length of each element of the traced streamlines. After the streamlines were traced using the fine lineStepWidth_, it gets resampled with this width.
  double targetLength_;           ///< the final length of the longest streamline, 0 means disabled
  double discardRelativeLength_;   ///< a relative length (in [0,1]), at the end streamlines are dropped that are smaller than this relative length times the median fibre length
  std::string csvFilename_;      ///< a csv output filename to write the node positions of the streamlines to (after postprocessing)
  std::string csvFilenameBeforePostprocessing_;      ///< a csv output filename to write the node positions of the streamlines to (before postprocessing)
};

};  // namespace

#include "postprocessing/streamline_tracer.tpp"