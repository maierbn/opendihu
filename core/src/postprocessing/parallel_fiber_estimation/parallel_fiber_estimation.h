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
  void generateParallelMeshRecursion(std::array<std::vector<std::vector<Vec3>>,4> &borderPoints, std::array<bool,4> subdomainIsAtBorder);

  //! take the streamlines at equidistant z points in streamlineZPoints and copy them to the array borderPointsSubdomain
  void rearrangeStreamlinePoints(std::vector<std::vector<Vec3>> &streamlineZPoints, std::array<std::array<std::vector<std::vector<Vec3>>,4>,8> &borderPointsSubdomain,
                                 std::array<std::array<std::vector<bool>,4>,8> &borderPointsSubdomainAreValid, std::array<bool,4> &subdomainIsAtBorder);

  //! interpolate points for invalid streamlines
  void fixIncompleteStreamlines(std::array<std::array<std::vector<std::vector<Vec3>>,4>,8> &borderPointsSubdomain, std::array<std::array<std::vector<bool>,4>,8> &borderPointsSubdomainAreValid, bool streamlineDirectionUpwards);

  //! refine the given border points (borderPointsOld) in x and y direction
  void refineBorderPoints(std::array<std::vector<std::vector<Vec3>>,4> &borderPointsOld, std::array<std::vector<std::vector<Vec3>>,4> &borderPoints);

  //! create the mesh with given borderPoints, using harmonic maps by calling the python script
  void createMesh(std::array<std::vector<std::vector<Vec3>>,4> &borderPoints, std::vector<Vec3> &nodePositions, std::array<int,3> &nElementsPerCoordinateDirectionLocal);

  //! check if the algorithm is at the stage where no more subdomains are created and the final fibers are traced
  bool checkTraceFinalFibers(int &level);

  //! create Dirichlet BC object
  void createDirichletBoundaryConditions(const std::array<int,3> &nElementsPerCoordinateDirectionLocal, std::shared_ptr<SpatialDiscretization::DirichletBoundaryConditions<FunctionSpaceType,1>> &dirichletBoundaryConditions);

  //! communicate ghost values for gradient and solution value to neighbouring processes, the ghost elements are obtained from the mesh partition
  void exchangeGhostValues(const std::array<bool,4> &subdomainIsAtBorder);

  //! create the seed points in form of a cross at the center of the current domain
  void createSeedPoints(const std::array<bool,4> &subdomainIsAtBorder, int seedPointsZIndex, const std::vector<Vec3> &nodePositions, std::vector<Vec3> &seedPoints);

  //! trace the fibers that are evenly distributed in the subdomain, this is the final step of the algorithm
  void traceResultFibers(double streamlineDirection, int seedPointsZIndex, const std::vector<Vec3> &nodePositions);

  //! trace the streamlines starting from the seed points, this uses functionality from the parent class
  void traceStreamlines(int nRanksZ, int rankZNo, double streamlineDirection, bool streamlineDirectionUpwards, std::vector<Vec3> &seedPoints, std::vector<std::vector<Vec3>> &streamlinePoints);

  //! sample the streamlines at equidistant z points, if the streamline does not run from bottom to top, only add seedPoint
  void sampleAtEquidistantZPoints(std::vector<std::vector<Vec3>> &streamlinePoints, const std::vector<Vec3> &seedPoints, std::vector<std::vector<Vec3>> &streamlineZPoints);

  //! fill in missing points at the borders, where no streamlines were traced
  void fillBorderPoints(std::array<std::vector<std::vector<Vec3>>,4> &borderPoints, std::array<std::array<std::vector<std::vector<Vec3>>,4>,8> &borderPointsSubdomain,
                        std::array<std::array<std::vector<bool>,4>,8> &borderPointsSubdomainAreValid, std::array<bool,4> &subdomainIsAtBorder);

  //! send border points to those ranks that will handle them in the next subdomain
  void sendBorderPoints(std::array<std::array<std::vector<std::vector<Vec3>>,4>,8> &borderPointsSubdomain, std::vector<MPI_Request> &sendRequests);

  //! receive the border points
  void receiveBorderPoints(int nRanksPerCoordinateDirectionPreviously, std::array<std::vector<std::vector<Vec3>>,4> &borderPointsNew, std::array<bool,4> &subdomainIsAtBorderNew);

  //! write all border points to a common file
  void outputBorderPoints(std::array<std::array<std::vector<std::vector<Vec3>>,4>,8> &borderPointsSubdomain, std::string name);

  const DihuContext context_;    ///< object that contains the python config for the current context and the global singletons meshManager and solverManager
  std::shared_ptr<FiniteElementMethodType> problem_;   ///< the DiscretizableInTime object that is managed by this class

  Data::ParallelFiberEstimation<FunctionSpaceType> data_;    ///< the data object that holds the gradient field variable

  //std::shared_ptr<FunctionSpaceType> functionSpace_;   ///< current function space / mesh

  PythonConfig specificSettings_;   ///< the specific python config for this module
  std::vector<Vec3> seedPositions_;  ///< the seed points from where the streamlines start

  std::string stlFilename_;   ///< the filename of the STL file
  double bottomZClip_;   ///< bottom z-value of the volume to consider
  double topZClip_;   ///< top z-value of the volume to consider
  int nBorderPointsX_;    ///< number of subdivisions of the line
  int nBorderPointsZ_;    ///< number of subdivisions in z direction
  int maxLevel_;   ///< the maximum level up to which the domain will be subdivided, number of final domains is 8^maxLevel_ (octree structure)
  int nBorderPointsXNew_;  ///< the value of nBorderPointsX_ in the next subdomain
  int nBorderPointsZNew_;  ///< the value of nBorderPointsZ_ in the next subdomain

  PyObject *moduleStlCreateMesh_;   ///< python module, file "stl_create_mesh.py"
  PyObject *moduleStlCreateRings_;   ///< python module, file "stl_create_rings.py"
  PyObject *moduleStlDebugOutput_;   ///< python module, file "stl_debug_output.py"
  PyObject *functionCreateRingSection_;  ///< python function create_ring_section
  PyObject *functionCreateRingSectionMesh_;  ///< python function create_ring_section_mesh
  PyObject *functionGetStlMesh_;  ///< python function get_stl_mesh
  PyObject *functionCreateRings_;               ///< create_rings
  PyObject *functionRingsToBorderPoints_;       ///< rings_to_border_points
  PyObject *functionBorderPointLoopsToList_;    ///< border_point_loops_to_list
  PyObject *functionOutputPoints_;              ///< output_points
  PyObject *functionOutputBorderPoints_;        ///< output_border_points
  PyObject *functionOutputGhostElements_;        ///< output_ghots_elements
  PyObject *functionCreate3dMeshFromBorderPointsFaces_;       ///< create_3d_mesh_from_border_points_faces

  std::shared_ptr<Partition::RankSubset> currentRankSubset_;  ///< the rank subset of the ranks that are used at the current stage of the algorithm
  std::array<int,3> nRanksPerCoordinateDirection_;   ///< the numbers of ranks in each coordinate direction at the current stage of the algorithm
  std::shared_ptr<Partition::MeshPartition<FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<3>, BasisFunctionType>>> meshPartition_; //< the mesh partition of this subdomain which contains information about the neighbouring ranks and the own index in the ranks grid

  OutputWriter::Manager outputWriterManager_; ///< manager object holding all output writer
};

};  // namespace

#include "postprocessing/parallel_fiber_estimation/parallel_fiber_estimation.tpp"
