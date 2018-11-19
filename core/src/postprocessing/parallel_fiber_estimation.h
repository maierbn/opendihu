#pragma once

#include <Python.h>  // has to be the first included header
#include <vector>

#include "function_space/function_space.h"
#include "postprocessing/streamline_tracer_base.h"
#include "interfaces/discretizable_in_time.h"
#include "interfaces/runnable.h"
#include "data_management/parallel_fiber_estimation.h"
#include "quadrature/gauss.h"

namespace Postprocessing
{

/** A class that creates a parallel mesh and at the same time traces streamlines through a Δu=0 field.
 */
template<typename BasisFunctionType>
class ParallelFiberEstimation :
  public StreamlineTracerBase<FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<3>,BasisFunctionType>>,
  public Runnable
{
public:
  //! constructor
  ParallelFiberEstimation(DihuContext context);

  //! initialize
  void initialize();

  //! run tracing of stream lines
  void run();

  //! function space to use, i.e. 3D structured deformable grid
  typedef FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<3>, BasisFunctionType> FunctionSpaceType;
  typedef SpatialDiscretization::FiniteElementMethod<
    Mesh::StructuredDeformableOfDimension<3>,
    BasisFunctionType,
    Quadrature::Gauss<3>,
    Equation::Static::Laplace
  > FiniteElementMethodType;

protected:

  //! perform the algorithm to recursively, collectively refine the mesh and trace fibers for the boundaries of the subdomains
  void generateParallelMesh();

  //! recursive part of the algorithm
  void generateParallelMeshRecursion(std::array<std::vector<std::vector<Vec3>>,4> &borderPoints, int level, std::array<bool,4> subdomainIsAtBorder);

  //! take the streamlines at equidistant z points in streamlineZPoints and copy them to the array borderPointsSubdomain
  void reorganizeStreamlinePoints(std::vector<std::vector<Vec3>> &streamlineZPoints, std::array<std::array<std::vector<std::vector<Vec3>>,4>,8> &borderPointsSubdomain, std::array<bool,4> &subdomainIsAtBorder);

  const DihuContext context_;    ///< object that contains the python config for the current context and the global singletons meshManager and solverManager
  std::shared_ptr<FiniteElementMethodType> problem_;   ///< the DiscretizableInTime object that is managed by this class

  Data::ParallelFiberEstimation<FunctionSpaceType> data_;    ///< the data object that holds the gradient field variable

  //std::shared_ptr<FunctionSpaceType> functionSpace_;   ///< current function space / mesh

  PythonConfig specificSettings_;   ///< the specific python config for this module
  std::vector<Vec3> seedPositions_;  ///< the seed points from where the streamlines start

  std::string stlFilename_;   ///< the filename of the STL file
  int bottomZClip_;   ///< bottom z-value of the volume to consider
  int topZClip_;   ///< top z-value of the volume to consider
  int nElementsZPerSubdomain_;   ///< number of elements per subdomain in z direction
  int nBorderPointsX_;    ///< number of subdivisions of the line

  PyObject* moduleStlCreateMesh_;   ///< python module, file "stl_create_mesh.py"
  std::shared_ptr<Partition::RankSubset> currentRankSubset_;  ///< the rank subset of the ranks that are used at the current stage of the algorithm

  OutputWriter::Manager outputWriterManager_; ///< manager object holding all output writer
};

};  // namespace

#include "postprocessing/parallel_fiber_estimation.tpp"
