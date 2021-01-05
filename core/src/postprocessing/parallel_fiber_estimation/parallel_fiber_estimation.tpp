#include "postprocessing/parallel_fiber_estimation/parallel_fiber_estimation.h"

#include <algorithm>
#include <petscvec.h>

#include "utility/python_utility.h"
#include "mesh/face_t.h"
#include "partition/mesh_partition/01_mesh_partition.h"
#include "spatial_discretization/dirichlet_boundary_conditions/01_dirichlet_boundary_conditions.h"

// write or load various checkpoints, this is for debugging to only run part of the algorithm on prescribed data
//#define USE_CHECKPOINT_BOUNDARY_POINTS
//#define USE_CHECKPOINT_MESH
//#define WRITE_CHECKPOINT_MESH
//#define WRITE_CHECKPOINT_BOUNDARY_POINTS
//#define WRITE_CHECKPOINT_GHOST_MESH
//#define USE_CHECKPOINT_GHOST_MESH

// output STl files for debugging
#define STL_OUTPUT                // output some stl files
//#define STL_OUTPUT_VERBOSE     // output more stl files

//#define FILE_COMMUNICATION        // when sending boundary points between ranks, do not use MPI but file I/O instead

// include files that implement various methods of this class, these make use the previous defines
#include "postprocessing/parallel_fiber_estimation/01_refine_boundary_points.tpp"
#include "postprocessing/parallel_fiber_estimation/02_create_mesh.tpp"
#include "postprocessing/parallel_fiber_estimation/03_create_dirichlet_boundary_conditions.tpp"
#include "postprocessing/parallel_fiber_estimation/03_create_neumann_boundary_conditions.tpp"
#include "postprocessing/parallel_fiber_estimation/04_exchange_ghost_values.tpp"
#include "postprocessing/parallel_fiber_estimation/05_create_seed_points.tpp"
#include "postprocessing/parallel_fiber_estimation/06_trace_streamlines.tpp"
#include "postprocessing/parallel_fiber_estimation/07_sample_at_equidistant_z_points.tpp"
#include "postprocessing/parallel_fiber_estimation/08_rearrange_streamline_points.tpp"
#include "postprocessing/parallel_fiber_estimation/09_fill_boundary_points.tpp"
#include "postprocessing/parallel_fiber_estimation/10_fix_incomplete_streamlines.tpp"
#include "postprocessing/parallel_fiber_estimation/11_fix_invalid.tpp"
#include "postprocessing/parallel_fiber_estimation/11_interpolate_fine_fibers_from_file.tpp"
#include "postprocessing/parallel_fiber_estimation/11_resample_fibers_in_file.tpp"
#include "postprocessing/parallel_fiber_estimation/11_trace_result_fibers.tpp"
#include "postprocessing/parallel_fiber_estimation/12_send_receive_boundary_points.tpp"
#include "postprocessing/parallel_fiber_estimation/adjust_filename.tpp"
#include "postprocessing/parallel_fiber_estimation/exchange_seed_points.tpp"
#include "postprocessing/parallel_fiber_estimation/output_boundary_points.tpp"

namespace Postprocessing
{

template<typename BasisFunctionType>
ParallelFiberEstimation<BasisFunctionType>::
ParallelFiberEstimation(DihuContext context) :
  context_(context["ParallelFiberEstimation"]), problem_(nullptr), data_(context_), specificSettings_(context_.getPythonConfig()), level_(0)
{
  LOG(TRACE) << "ParallelFiberEstimation::ParallelFiberEstimation()";

  outputWriterManager_.initialize(context_, specificSettings_);

  inputMeshFilename_ = specificSettings_.getOptionString("inputMeshFilename", "");
  resultFilename_ = specificSettings_.getOptionString("resultFilename", "fibers.bin");
  bottomZClip_ = specificSettings_.getOptionDouble("bottomZClip", 0);
  topZClip_ = specificSettings_.getOptionDouble("topZClip", 100);
  nBoundaryPointsZ_ = specificSettings_.getOptionInt("nElementsZPerSubdomain", 13)+1;
  nBoundaryPointsX_ = specificSettings_.getOptionInt("nElementsXPerSubdomain", 13)+1;
  maxLevel_ = specificSettings_.getOptionInt("maxLevel", 2);
  nFineGridFibers_ = specificSettings_.getOptionInt("nFineGridFibers", 2);
  improveMesh_ = specificSettings_.getOptionBool("improveMesh", true);
  useNeumannBoundaryConditions_ = specificSettings_.getOptionBool("useNeumannBoundaryConditions", false);
  nNodesPerFiber_ = specificSettings_.getOptionInt("nNodesPerFiber", 1000);
  finalBottomZClip_ = specificSettings_.getOptionDouble("finalBottomZClip", bottomZClip_);
  finalTopZClip_ = specificSettings_.getOptionDouble("finalTopZClip", topZClip_);
  refinementFactors_ = specificSettings_.getOptionArray<int,3>("refinementFactors", std::array<int,3>({1,1,1}));
  laplacianSmoothingNIterations_ = specificSettings_.getOptionInt("laplacianSmoothingNIterations", 10);
  ghostLayerWidth_ = specificSettings_.getOptionInt("ghostLayerWidth", 1);
  maxAreaFactor_ = specificSettings_.getOptionDouble("maxAreaFactor", 100.0);

  this->lineStepWidth_ = specificSettings_.getOptionDouble("lineStepWidth", 1e-2, PythonUtility::Positive);
  this->maxNIterations_ = specificSettings_.getOptionInt("maxIterations", 100000, PythonUtility::Positive);
  this->useGradientField_ = specificSettings_.getOptionBool("useGradientField", false);

  if ((nBoundaryPointsX_-1) % 2 != 0)
  {
    LOG(WARNING) << "Value for nElementsXPerSubdomain (" << (nBoundaryPointsX_-1) << ") is odd, should be even. "
      << "Now changing to " << 2*int(nBoundaryPointsX_/2);
  }

  // ensure nBoundaryPointsX is odd
  nBoundaryPointsX_ = 2*int(nBoundaryPointsX_/2)+1;
  // ensure nBoundaryPointsZ is odd
  nBoundaryPointsZ_ = 2*int(nBoundaryPointsZ_/2)+1;


  // run stl_create_mesh.rings_to_boundary_points
  // "Standardize every ring to be in counter-clockwise direction and starting with the point with lowest x coordinate, then sample boundary points"
}

template<typename BasisFunctionType>
void ParallelFiberEstimation<BasisFunctionType>::
initialize()
{
  LOG(TRACE) << "ParallelFiberEstimation::initialize";

  currentRankSubset_ = std::make_shared<Partition::RankSubset>(0);  // rankSubset with only the master rank (rank no 0)

  // load python modules
  moduleStlCreateRings_ = PyImport_ImportModule("stl_create_rings");
  PythonUtility::checkForError();
  if (moduleStlCreateRings_ == NULL)
  {
    LOG(FATAL) << "Could not load python module stl_create_rings. Ensure that the PYTHONPATH environment variable contains the path to opendihu/scripts/geometry_manipulation" << std::endl
      << "Execute the following: " << std::endl << "export PYTHONPATH=$PYTHONPATH:" << OPENDIHU_HOME << "/scripts/geometry_manipulation";
  }

  moduleStlCreateMesh_ = PyImport_ImportModule("stl_create_mesh");
  PythonUtility::checkForError();
  assert(moduleStlCreateMesh_);

  moduleStlDebugOutput_ = PyImport_ImportModule("stl_debug_output");
  PythonUtility::checkForError();
  assert(moduleStlDebugOutput_);

  // load functions from those modules
  // actual functions
  functionCreateRingSection_ = PyObject_GetAttrString(moduleStlCreateRings_, "create_ring_section");
  assert(functionCreateRingSection_);

  functionCreateBoundaryPoints_ = PyObject_GetAttrString(moduleStlCreateRings_, "create_boundary_points");
  assert(functionCreateBoundaryPoints_);

  // debugging functions
  functionOutputPoints_ = PyObject_GetAttrString(moduleStlDebugOutput_, "output_points");
  assert(functionOutputPoints_);

  functionOutputBoundaryPoints_ = PyObject_GetAttrString(moduleStlDebugOutput_, "output_boundary_points");
  assert(functionOutputBoundaryPoints_);

  functionOutputStreamline_ = PyObject_GetAttrString(moduleStlDebugOutput_, "output_streamline");
  assert(functionOutputStreamline_);

  functionOutputStreamlines_ = PyObject_GetAttrString(moduleStlDebugOutput_, "output_streamlines");
  assert(functionOutputStreamlines_);

  functionOutputRings_ = PyObject_GetAttrString(moduleStlDebugOutput_, "output_rings");
  assert(functionOutputRings_);

  functionOutputGhostElements_ = PyObject_GetAttrString(moduleStlDebugOutput_, "output_ghost_elements");
  assert(functionOutputGhostElements_);

  functionCreate3dMeshFromBoundaryPointsFaces_ = PyObject_GetAttrString(moduleStlCreateMesh_, "create_3d_mesh_from_boundary_points_faces");
  assert(functionCreate3dMeshFromBoundaryPointsFaces_);
}

template<typename BasisFunctionType>
void ParallelFiberEstimation<BasisFunctionType>::
run()
{
  initialize();

  generateParallelMesh();

  // output
  outputWriterManager_.writeOutput(data_);
}

template<typename BasisFunctionType>
void ParallelFiberEstimation<BasisFunctionType>::
generateParallelMesh()
{
  LOG(DEBUG) << "generateParallelMesh, nBoundaryPointsX: " << nBoundaryPointsX_;

  //MPI_Barrier(MPI_COMM_WORLD);

  // prepare rankSubset
  //currentRankSubset_ = std::make_shared<Partition::RankSubset>(0);  // rankSubset with only the master rank (rank no 0)
  nRanksPerCoordinateDirection_.fill(1);

  // get loops of whole domain

  std::array<std::vector<std::vector<Vec3>>,4> boundaryPoints;  // boundaryPoints[face_t][z-level][pointIndex]
  std::array<bool,4> subdomainIsAtBoundary;

  loadInitialCheckpoints(boundaryPoints, subdomainIsAtBoundary);

  generateParallelMeshRecursion(boundaryPoints, subdomainIsAtBoundary);
}

template<typename BasisFunctionType>
void ParallelFiberEstimation<BasisFunctionType>::
generateParallelMeshRecursion(std::array<std::vector<std::vector<Vec3>>,4> &boundaryPointsOld, std::array<bool,4> subdomainIsAtBoundary)
{
  // ----------------------------
  // algorithm:
  // create mesh in own domain, using python, harmonic maps

  // solve laplace problem globally

  // communicate ghost elements to neighbouring subdomains

  // trace fibers to determine new subdomain boundaries, also on the outside (shared between domains)

  // create subdomains
    // create new communicator
    // communicate all old elements to the processes of the new communicator
    // on the new processes create new meshes using the coarse data

  // call method recursively


  LOG(DEBUG) << "generateParallelMeshRecursion, n boundary points: "
    << boundaryPointsOld[0].size() << "," << boundaryPointsOld[1].size() << "," << boundaryPointsOld[2].size() << "," << boundaryPointsOld[3].size();

#ifndef NDEBUG
#ifdef STL_OUTPUT
//#ifdef STL_OUTPUT_VERBOSE
  if (currentRankSubset_->ownRankIsContained())
  {
    PyObject_CallFunction(functionOutputBoundaryPoints_, "s i i O f", "01_boundary_points_old", currentRankSubset_->ownRankNo(), level_,
                          PythonUtility::convertToPython<std::array<std::vector<std::vector<Vec3>>,4>>::get(boundaryPointsOld), 0.2);
    PythonUtility::checkForError();
  }

#endif
//#endif
#endif

  // check if all vectors are filled with points
  for (int face = (int)Mesh::face_t::face0Minus; face <= (int)Mesh::face_t::face1Plus; face++)
  {
    for (int zLevelIndex = 0; zLevelIndex < boundaryPointsOld[face].size(); zLevelIndex++)
    {
      if (boundaryPointsOld[face][zLevelIndex].empty())
      {
        LOG(ERROR) << "boundaryPoints face: " << Mesh::getString((Mesh::face_t)face) << " z: " << zLevelIndex << " contains no points.";
      }
    }
  }

  //! new boundary points for the next subdomain, these are computed on the local subdomain and then send to the sub-subdomains
  std::array<std::array<std::vector<std::vector<Vec3>>,4>,8> boundaryPointsSubdomain;  // [subdomain index][face_t][z-level][point index]

  // check if the algorithm is at the stage where no more subdomains are created and the final fibers are traced
  // level_ = log2(nRanksPerCoordinateDirection_)
  bool traceFinalFibers = checkTraceFinalFibers();
  LOG(DEBUG) << "level " << level_;

  bool refineSubdomainsOnThisRank = false;

  if (!currentRankSubset_->ownRankIsContained())
  {
    LOG(DEBUG) << "own rank is not contained in currentRankSubset_";
  }

  // here we have the boundary points of our subdomain, now the boundarys that split this subdomain in 8 more are computed
  if (currentRankSubset_->ownRankIsContained())
  {
    assert(!boundaryPointsOld[0].empty());

    refineSubdomainsOnThisRank = true;
    nBoundaryPointsXNew_ = nBoundaryPointsX_*2-1;   // number of boundary points per coordinate direction of the logical rectangle, for the subdomain

    // refine the given boundary points in x and y direction, i.e. horizontally, to yield twice the number of points on each boundary
    std::array<std::vector<std::vector<Vec3>>,4> boundaryPoints;    // [face_t][z-level][pointIndex]
    refineBoundaryPoints(boundaryPointsOld, boundaryPoints);

    // dimensions of boundaryPoints[face_t][z-level][pointIndex]: boundaryPoints[4][nBoundaryPointsZ_][nBoundaryPointsXNew_]

    // create 3D mesh in own subdomain, using python, algorithm with harmonic maps
    std::vector<Vec3> nodePositions;
    std::array<int,3> nElementsPerCoordinateDirectionLocal;
    createMesh(boundaryPoints, nodePositions, nElementsPerCoordinateDirectionLocal);

    std::array<global_no_t,3> nElementsPerCoordinateDirectionGlobal;

    // create meshPartition for function space with the currently used ranks
    context_.partitionManager()->setRankSubsetForNextCreatedPartitioning(currentRankSubset_);

    std::vector<int> rankNos;
    //std::shared_ptr<Partition::MeshPartition<FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<3>, BasisFunctionType>>> meshPartition
    meshPartition_ = context_.partitionManager()->template createPartitioningStructuredLocal<FunctionSpaceType>(
        this->specificSettings_, nElementsPerCoordinateDirectionGlobal, nElementsPerCoordinateDirectionLocal, nRanksPerCoordinateDirection_, rankNos);

    LOG(DEBUG) << "nElementsPerCoordinateDirectionGlobal: " << nElementsPerCoordinateDirectionGlobal;

    // output ghost indices
    std::stringstream stream2;
    std::vector<PetscInt> ghostDofNosGlobalPetsc = meshPartition_->ghostDofNosGlobalPetsc();
    for (int i = 0; i < ghostDofNosGlobalPetsc.size(); i++)
    {
      stream2 << ghostDofNosGlobalPetsc[i] << " ";
    }
    LOG(DEBUG) << "ghostDofNosGlobalPetsc: " << stream2.str();

    // create function space
    std::stringstream meshName;
    meshName << "meshLevel" << level_;

    int nNodePositionsWithoutGhostsX = nBoundaryPointsXNew_ - 1;
    int nNodePositionsWithoutGhostsY = nBoundaryPointsXNew_ - 1;
    int nNodePositionsWithoutGhostsZ = nBoundaryPointsZ_ - 1;

    if (meshPartition_->hasFullNumberOfNodes(0))
      nNodePositionsWithoutGhostsX += 1;
    if (meshPartition_->hasFullNumberOfNodes(1))
      nNodePositionsWithoutGhostsY += 1;
    if (meshPartition_->hasFullNumberOfNodes(2))
      nNodePositionsWithoutGhostsZ += 1;

    int nNodePositionsWithoutGhosts = nNodePositionsWithoutGhostsX * nNodePositionsWithoutGhostsY * nNodePositionsWithoutGhostsZ;
    std::vector<Vec3> nodePositionsWithoutGhosts(nNodePositionsWithoutGhosts);

    // extract the node positions without ghosts from all node positions
    LOG(DEBUG) << "nNodePositionsWithoutGhosts: " << nNodePositionsWithoutGhosts << ", nNodePositionsWithGhosts: " << nodePositions.size();
    LOG(DEBUG) << nNodePositionsWithoutGhostsX << "x" << nNodePositionsWithoutGhostsY << "x" << nNodePositionsWithoutGhostsZ;
    LOG(DEBUG) << nBoundaryPointsXNew_ << "x" << nBoundaryPointsXNew_ << "x" << nBoundaryPointsZ_ << "=" << nBoundaryPointsZ_*nBoundaryPointsXNew_*nBoundaryPointsXNew_;

    for (dof_no_t z = 0; z < nBoundaryPointsZ_; z++)
    {
      for (dof_no_t y = 0; y < nBoundaryPointsXNew_; y++)
      {
        for (dof_no_t x = 0; x < nBoundaryPointsXNew_; x++)
        {
          dof_no_t indexWithoutGhosts = z*nNodePositionsWithoutGhostsX*nNodePositionsWithoutGhostsY + y*nNodePositionsWithoutGhostsX + x;
          dof_no_t indexWithGhosts = z*nBoundaryPointsXNew_*nBoundaryPointsXNew_ + y*nBoundaryPointsXNew_ + x;

          if (z < nNodePositionsWithoutGhostsZ && y < nNodePositionsWithoutGhostsY && x < nNodePositionsWithoutGhostsX)
          {
            assert(indexWithGhosts < nodePositions.size());
            assert(indexWithoutGhosts < nodePositionsWithoutGhosts.size());

            nodePositionsWithoutGhosts[indexWithoutGhosts] = nodePositions[indexWithGhosts];
          }
        }
      }
    }

    // create the new function space
    context_.partitionManager()->setRankSubsetForNextCreatedPartitioning(currentRankSubset_);
    this->functionSpace_ = context_.meshManager()->template createFunctionSpaceWithGivenMeshPartition<FunctionSpaceType>(
      meshName.str(), meshPartition_, nodePositionsWithoutGhosts, nElementsPerCoordinateDirectionLocal, nRanksPerCoordinateDirection_, false);

    LOG(DEBUG) << "after refinement: ";
    LOG(DEBUG) << "n nodePositions with ghosts: " << nodePositions.size();
    LOG(DEBUG) << "n nodePositions without ghosts: " << nodePositionsWithoutGhosts.size();
    LOG(DEBUG) << "n dofs local without ghosts: " << this->functionSpace_->meshPartition()->nDofsLocalWithoutGhosts();
    LOG(DEBUG) << "nElementsPerCoordinateDirectionLocal: " << nElementsPerCoordinateDirectionLocal;

    if (refinementFactors_[0] > 1 || refinementFactors_[1] > 1 || refinementFactors_[2] > 1)
    {
      this->functionSpace_->refineMesh(refinementFactors_);
    }

    LOG(DEBUG) << "after refinement: ";
    LOG(DEBUG) << "n dofs local without ghosts: " << this->functionSpace_->meshPartition()->nDofsLocalWithoutGhosts();

    /*
    std::array<global_no_t,FunctionSpace::dim()> &nElementsGlobal,
      const std::array<element_no_t,FunctionSpace::dim()> nElementsLocal,
      const std::array<int,FunctionSpace::dim()> nRanks
    * */

    // solve laplace problem globally
    // create problem
    problem_ = std::make_shared<FiniteElementMethodType>(context_, this->functionSpace_);
    data_.setProblem(problem_);

#ifndef NDEBUG
#ifdef STL_OUTPUT
    std::vector<Vec3> geometryValues;
    problem_->data().functionSpace()->geometryField().getValuesWithGhosts(geometryValues);

    PyObject_CallFunction(functionOutputPoints_, "s i i O f", "03_geometry_values", currentRankSubset_->ownRankNo(), level_,
                          PythonUtility::convertToPython<std::vector<Vec3>>::get(geometryValues), 0.1);
    PythonUtility::checkForError();
#endif
#endif

    // create boundary condition objects, either dirichlet or neumann
    std::shared_ptr<SpatialDiscretization::DirichletBoundaryConditions<FunctionSpaceType,1>> dirichletBoundaryConditions;
    std::shared_ptr<SpatialDiscretization::NeumannBoundaryConditions<FunctionSpaceType,Quadrature::Gauss<3>,1>> neumannBoundaryConditions;
    if (useNeumannBoundaryConditions_)
    {
      // create neumann boundary condition object, using this->functionSpace_
      createNeumannBoundaryConditions(neumannBoundaryConditions);
    }
    else
    {
      // create dirichlet boundary condition object, using this->functionSpace_
      createDirichletBoundaryConditions(dirichletBoundaryConditions);
    }

    // initialize data object to allocate gradient field
    data_.setFunctionSpace(this->functionSpace_);
    data_.reset();
    data_.initialize();

    // set dirichlet values in data for debugging output of vector of dirichlet BC values
    if (!useNeumannBoundaryConditions_)
    {
      const std::vector<dof_no_t> &dirichletDofs = dirichletBoundaryConditions->boundaryConditionNonGhostDofLocalNos();
      const std::vector<VecD<1>> &dirichletValues = dirichletBoundaryConditions->boundaryConditionValues();
      data_.dirichletValues()->setValues(-1);
      data_.dirichletValues()->setValues(dirichletDofs, dirichletValues);
    }

    problem_->reset();

    // set boundary conditions to the problem
    if (useNeumannBoundaryConditions_)
    {
      problem_->setNeumannBoundaryConditions(neumannBoundaryConditions);
    }
    else
    {
      problem_->setDirichletBoundaryConditions(dirichletBoundaryConditions);
    }

    problem_->initialize();

    // solve the laplace problem, globally on all subdomains on all ranks of the current rank subset
    problem_->run();

    // compute a gradient field from the solution
    problem_->data().solution()->computeGradientField(data_.gradient(), data_.jacobianConditionNumber());

    data_.gradient()->setRepresentationGlobal();

    // output the results
    std::stringstream filename;
    filename << "out/level_" << level_ << "/biceps";
    this->outputWriterManager_.setFilename(filename.str());
    this->outputWriterManager_.writeOutput(data_);  // output gradient

    //if (level > 1)
    //  LOG(FATAL) << "end";

    // exchange layer of ghost elements (not nodes) between neighbouring ranks,
    // then every subdomain knows one layer of elements around it (only in x/y direction)
    LOG(DEBUG) << "\nConstruct ghost elements";

    exchangeGhostValues(subdomainIsAtBoundary);

    // initialize values in base class
    this->solution_ = problem_->data().solution();
    this->gradient_ = data_.gradient();

    // determine seed points for streamlines
    // nodePositions contains all node positions in the current 3D mesh
    // naming of boundarys: (0-,1-,0+,1+)
    //     ___1+__
    //    |   |   |
    // 0- |___|___| 0+
    // ^  |   |   |
    // |  |___|___|
    // +-->   1-
    std::vector<Vec3> seedPoints;
    // seedPoints contains in this order:
    // face0Minus, face0Plus, face1Minus (with corner points), face1Plus (with corner points),
    // horizontal center line (with corner points), vertical center line (with corner points, with center point)

    int rankZNo = meshPartition_->ownRankPartitioningIndex(2);
    int nRanksZ = meshPartition_->nRanks(2);
    bool streamlineDirectionUpwards = rankZNo >= nRanksZ/2;
    // 2^level = nRanksz
    // level 0: 1 rank, nRanksZ/2=0, upwards
    // level 1: 2x2x2 ranks, nRanksZ/2=1,  rankZNo=0 | rankZNo=1
    // level 2: 4x4x4 ranks, nRanksZ/2=2,  rankZNo=0,1 | rankZNo=2,3

    LOG(DEBUG) << " nRanksZ: " << nRanksZ;

    int seedPointsZIndex = 0;
    double streamlineDirection = 1.0;
    if (!streamlineDirectionUpwards)
    {
      seedPointsZIndex = nBoundaryPointsZ_-1;
      streamlineDirection = -1.0;
    }

    if (nRanksZ == 1)
    {
      // if there is only one rank, start the streamlines at the center and trace towards top and bottom
      seedPointsZIndex = nBoundaryPointsZ_/2;
    }

    LOG(DEBUG) << "createSeedPoints";

    // determine the seed points of the streamlines
    createSeedPoints(subdomainIsAtBoundary, seedPointsZIndex, nodePositions, seedPoints);

    LOG(DEBUG) << "traceStreamlines";
    LOG(INFO) << "tracing streamlines ...";

    // trace streamlines from seed points, this also exchanges the seed points
    std::vector<std::vector<Vec3>> streamlinePoints;
    traceStreamlines(nRanksZ, rankZNo, streamlineDirection, streamlineDirectionUpwards, seedPoints, streamlinePoints);

    // sample streamlines at equidistant z points
    nBoundaryPointsZNew_ = nBoundaryPointsZ_*2 - 1;

    std::vector<std::vector<Vec3>> streamlineZPoints;
    sampleAtEquidistantZPoints(streamlinePoints, seedPoints, streamlineZPoints);

    // save streamline points to the portions of the faces
    //std::array<std::array<std::vector<std::vector<Vec3>>,4>,8> boundaryPointsSubdomain;  // [subdomain index][face_t][z-level][point index]

    LOG(DEBUG) << "rearrangeStreamlinePoints";

    // assign sampled points to the data structure boundaryPointsSubdomain, which contains the points for each subdomain and face, as list of points for each z level
    std::array<std::array<std::vector<bool>,4>,8> boundaryPointsSubdomainAreValid;
    std::array<std::vector<Vec3>,4> cornerStreamlines;
    rearrangeStreamlinePoints(streamlineZPoints, boundaryPointsSubdomain, cornerStreamlines, boundaryPointsSubdomainAreValid, subdomainIsAtBoundary);

    LOG(DEBUG) << "nBoundaryPointsZ_: " << nBoundaryPointsZ_ << ", nBoundaryPointsZNew_: " << nBoundaryPointsZNew_;

    // write boundary points to file
    outputStreamlines(boundaryPointsSubdomain, "09_traced");

    // fill the streamline points that are at the boundary
    fillBoundaryPoints(boundaryPoints, boundaryPointsSubdomain, cornerStreamlines, boundaryPointsSubdomainAreValid, subdomainIsAtBoundary);

    // write boundary points to file
    outputStreamlines(boundaryPointsSubdomain, "10_filled");

    // fill in missing points
    fixIncompleteStreamlines(boundaryPointsSubdomain, boundaryPointsSubdomainAreValid, streamlineDirectionUpwards, subdomainIsAtBoundary, boundaryPoints);

    // write boundary points to file
    outputStreamlines(boundaryPointsSubdomain, "11_fixed");

    LOG(DEBUG) << "call exchangeSeedPointsAfterTracing with " << seedPoints.size()
      << " seed points, " << streamlinePoints.size() << " streamlines from traceStreamlines";

    // send end points of streamlines to next rank that continues the streamline
    exchangeBoundarySeedPointsAfterTracing(nRanksZ, rankZNo, streamlineDirectionUpwards, subdomainIsAtBoundary, boundaryPointsSubdomain, cornerStreamlines);

    // if this is the end of the recursion, trace final fibers, reuse the fibers at the boundarys that were already traced
    if (traceFinalFibers)
    {
      traceResultFibers(streamlineDirection, seedPointsZIndex, nodePositions, boundaryPointsSubdomain, true);

      //MPI_Barrier(this->currentRankSubset_->mpiCommunicator());
      LOG(DEBUG) << "algorithm is finished";

      // algorithm is finished
      return;
    }
    else
    {
      // trace resulting fibers, for debugging output
      traceResultFibers(streamlineDirection, seedPointsZIndex, nodePositions, boundaryPointsSubdomain, false);
    }
  }  // if own rank is part of this stage of the algorithm

  LOG(DEBUG) << "create subdomains";

  //MPI_Barrier(MPI_COMM_WORLD);

  // create subdomains
  // create new rank subset i.e. a new MPI communicator
  int nRanksPerCoordinateDirectionPreviously = nRanksPerCoordinateDirection_[0];
  nRanksPerCoordinateDirection_.fill(nRanksPerCoordinateDirectionPreviously*2);

  int nRanks = nRanksPerCoordinateDirection_[0] * nRanksPerCoordinateDirection_[1] * nRanksPerCoordinateDirection_[2];
  std::vector<int> ranks(nRanks);
  std::iota(ranks.begin(), ranks.end(), 0);
  currentRankSubset_ = std::make_shared<Partition::RankSubset>(ranks.begin(), ranks.end(), context_.rankSubset());

  // sendBoundaryPoints and receiveBoundaryPoints need the same level
  determineLevel();
  LOG(DEBUG) << "new currentRankSubset_ created, refineSubdomainsOnThisRank: " << refineSubdomainsOnThisRank << ", rankSubset: " << *currentRankSubset_ << ", new level: " << level_;

  // send boundary points to ranks that will handle the new subdomains
  std::vector<MPI_Request> sendRequests;
  std::vector<std::vector<double>> sendBuffers;
  if (refineSubdomainsOnThisRank)     // if this rank has created new subdomains (was part of the previous rank subset)
  {
    sendBoundaryPoints(boundaryPointsSubdomain, sendBuffers, sendRequests);
  }

  // receive boundary points
  std::array<std::vector<std::vector<Vec3>>,4> boundaryPointsNew;
  std::array<bool,4> subdomainIsAtBoundaryNew;

  if (currentRankSubset_->ownRankIsContained())
  {
    receiveBoundaryPoints(nRanksPerCoordinateDirectionPreviously, boundaryPointsNew, subdomainIsAtBoundaryNew);
  }

#ifndef FILE_COMMUNICATION
  MPIUtility::handleReturnValue(MPI_Waitall(sendRequests.size(), sendRequests.data(), MPI_STATUSES_IGNORE), "MPI_Waitall");
#endif

  sendBuffers.clear();
  //LOG(FATAL) << "done";

  // call method recursively
  generateParallelMeshRecursion(boundaryPointsNew, subdomainIsAtBoundaryNew);

/*
  LOG(FATAL) << "SUCCESS";*/
}

template<typename BasisFunctionType>
void ParallelFiberEstimation<BasisFunctionType>::
loadInitialCheckpoints(std::array<std::vector<std::vector<Vec3>>,4> &boundaryPoints, std::array<bool,4> &subdomainIsAtBoundary)
{
  bool useCheckpointsBoundaryPoints = false;

#ifdef USE_CHECKPOINT_BOUNDARY_POINTS     // normal execution

  // start at higher level
  level_ = 0;
  std::vector<int> ranks((level_+1)*(level_+1)*(level_+1));
  std::iota(ranks.begin(), ranks.end(), 0);
  currentRankSubset_ = std::make_shared<Partition::RankSubset>(ranks.begin(), ranks.end());
  nRanksPerCoordinateDirection_.fill(level_+1);
  determineLevel();

  // parse file
  int subdomainIndex = currentRankSubset_->ownRankNo();
  std::stringstream filename;
  filename << "checkpoints/checkpoint_boundaryPoints_l" << level_ << "_subdomain_" << subdomainIndex << ".csv";
  std::ifstream file(filename.str().c_str(), std::ios::in);
  if (!file.is_open())
  {
    LOG(ERROR) << "Could not open file \"" << filename.str() << "\" for reading'";
  }
  else
  {
    useCheckpointsBoundaryPoints = true;
    assert(file.is_open());

    char c;
    int size1, size2;
    file >> size1 >> c >> size2 >> c;
    for (int i = 0; i < 4; i++)
    {
      int value;
      file >> value >> c;
      subdomainIsAtBoundary[i] = (value == 1);
    }
    std::string line;
    std::getline(file, line);


    LOG(DEBUG) << "parse values from file \"" << filename.str() << "\".";
    LOG(DEBUG) << "size1: " << size1 << ", size2: " << size2;


    for (int faceNo = (int)Mesh::face_t::face0Minus; faceNo <= (int)Mesh::face_t::face1Plus; faceNo++)
    {
      boundaryPoints[faceNo].resize(size1);
      for (int zLevelIndex = 0; zLevelIndex < boundaryPoints[faceNo].size(); zLevelIndex++)
      {
        boundaryPoints[faceNo][zLevelIndex].resize(size2);
        std::string line;
        std::getline(file, line);

        for (int pointIndex = 0; pointIndex < boundaryPoints[faceNo][zLevelIndex].size(); pointIndex++)
        {
          for (int i = 0; i < 3; i++)
          {
            boundaryPoints[faceNo][zLevelIndex][pointIndex][i] = atof(line.substr(0, line.find(";")).c_str());
            line = line.substr(line.find(";")+1);
          }
        }
      }
    }

    file.close();

    LOG(DEBUG) << "subdomainIndex: " << subdomainIndex;
    LOG(DEBUG) << "subdomainIsAtBoundary: " << subdomainIsAtBoundary;
    // LOG(DEBUG) << "boundaryPoints: " << boundaryPoints;



    //MPI_Barrier(MPI_COMM_WORLD);
  }
#endif

  // if checkpoint file could not be opened, continue with normal execution
  if (!useCheckpointsBoundaryPoints)
  {
    // only rank 0 creates the first boundary points
    if (DihuContext::ownRankNoCommWorld() == 0)
    {
      // run python script to generate loops for the whole volume

      PyObject *boundaryPointsPy = PyObject_CallFunction(functionCreateBoundaryPoints_, "s f f i i",
        inputMeshFilename_.c_str(), bottomZClip_, topZClip_, nBoundaryPointsZ_, 4*(nBoundaryPointsX_-1));
      PythonUtility::checkForError();
      assert(boundaryPointsPy);

      std::vector<std::vector<Vec3>> loops = PythonUtility::convertFromPython<std::vector<std::vector<Vec3>>>::get(boundaryPointsPy);
      //std::vector<double> lengths = PythonUtility::convertFromPython<std::vector<double>>::get(lengthsPy);
      //LOG(DEBUG) << "loops: " << loops << " (size: " << loops.size() << ")" << ", lengths: " << lengths;

  #ifndef NDEBUG
  #ifdef STL_OUTPUT
      // output the loops
      PyObject_CallFunction(functionOutputRings_, "s i i O f", "00_loops", currentRankSubset_->ownRankNo(), level_,
                            PythonUtility::convertToPython<std::vector<std::vector<Vec3>>>::get(loops), 0.1);
      PythonUtility::checkForError();
  #endif
  #endif

      // rearrange the boundary points from the loops to the portions of the faces
      //std::array<std::vector<std::vector<Vec3>>,4> boundaryPoints;  // boundaryPoints[face_t][z-level][pointIndex]

      // loop over faces: face0Minus = 0, face0Plus, face1Minus, face1Plus

      //   ^ --(1+)-> ^
      // ^ 0-         0+
      // | | --(1-)-> |
      // +-->

      for (int face = Mesh::face_t::face0Minus; face <= Mesh::face_t::face1Plus; face++)
      {
        boundaryPoints[face].resize(nBoundaryPointsZ_);
        for (int zIndex = 0; zIndex < nBoundaryPointsZ_; zIndex++)
        {
          boundaryPoints[face][zIndex].resize(nBoundaryPointsX_);

          VLOG(1) << "face " << face << ", zIndex " << zIndex << " nloops: " << loops[zIndex].size();

          if (face == Mesh::face_t::face1Minus)
          {
            std::copy(loops[zIndex].begin(), loops[zIndex].begin() + (nBoundaryPointsX_-1)+1, boundaryPoints[face][zIndex].begin());
          }
          else if (face == Mesh::face_t::face0Plus)
          {
            std::copy(loops[zIndex].begin() + (nBoundaryPointsX_-1), loops[zIndex].begin() + 2*(nBoundaryPointsX_-1)+1, boundaryPoints[face][zIndex].begin());
          }
          else if (face == Mesh::face_t::face1Plus)
          {
            std::reverse_copy(loops[zIndex].begin() + 2*(nBoundaryPointsX_-1), loops[zIndex].begin() + 3*(nBoundaryPointsX_-1)+1, boundaryPoints[face][zIndex].begin());
          }
          else if (face == Mesh::face_t::face0Minus)
          {
            std::reverse_copy(loops[zIndex].begin() + 3*(nBoundaryPointsX_-1), loops[zIndex].begin() + 4*(nBoundaryPointsX_-1), boundaryPoints[face][zIndex].begin()+1);
            boundaryPoints[face][zIndex][0] = loops[zIndex].front();
          }
        }
      }

      subdomainIsAtBoundary[Mesh::face_t::face0Minus] = true;
      subdomainIsAtBoundary[Mesh::face_t::face0Plus] = true;
      subdomainIsAtBoundary[Mesh::face_t::face1Minus] = true;
      subdomainIsAtBoundary[Mesh::face_t::face1Plus] = true;
    }
  } // if not using checkpoint data

}

} // namespace
