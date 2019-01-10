#pragma once

#include <Python.h>  // has to be the first included header
#include <vector>

#include "interfaces/discretizable_in_time.h"
#include "interfaces/runnable.h"
#include "data_management/streamline_tracer.h"

namespace Postprocessing
{

/** A class that traces streamlines through a given solution field. This base class only performs the tracing of the streamlines.
 */
template<typename FunctionSpace>
class StreamlineTracerBase
{
public:

  //! trace the streamline starting from startingPoint in the element initialElementNo, direction is either 1. or -1. depending on the direction
  void traceStreamline(Vec3 startingPoint, double direction, std::vector<Vec3> &points);

protected:

  std::shared_ptr<FunctionSpace> functionSpace_;   ///< function space of the solution field in which the tracing is performed
  std::shared_ptr<FieldVariable::FieldVariable<FunctionSpace,1>> solution_;   ///< solution field in which the tracing is performed
  std::shared_ptr<FieldVariable::FieldVariable<FunctionSpace,3>> gradient_;   ///< gradient field which can be used to trace the streamlines (if useGradient_ is set to true)

  std::array<std::shared_ptr<FunctionSpace>,4> ghostMesh_;   ///< the ghost  meshes around the subdomain elements, for faces Mesh::face_t::face0Minus, face0Plus, face1Minus, face1Plus
  std::array<std::shared_ptr<FieldVariable::FieldVariable<FunctionSpace,3>>,4> ghostMeshGradient_;    ///< gradient field in ghost meshes, ghost meshes are surrounding the regular subdomain by one layer of elements
  std::array<std::shared_ptr<FieldVariable::FieldVariable<FunctionSpace,1>>,4> ghostMeshSolution_;    ///< solution field in ghost meshes, ghost meshes are surrounding the regular subdomain by one layer of elements

  double lineStepWidth_;     ///< the line step width used for integrating the streamlines

  int maxNIterations_;   ///< the maximum number of iterations to trace for a streamline
  bool useGradientField_;  ///< There are 2 implementations of streamline tracing. The first one (useGradientField_) uses a precomputed gradient field that is interpolated linearly and the second uses the gradient directly from the Laplace solution field. // The first one seems more stable, because the gradient is zero and the position of the boundary conditions.
};

};  // namespace

#include "postprocessing/streamline_tracer_base.tpp"
