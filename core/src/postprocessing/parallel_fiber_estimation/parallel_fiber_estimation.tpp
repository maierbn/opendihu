#include "postprocessing/parallel_fiber_estimation/parallel_fiber_estimation.h"

#include <algorithm>
#include <petscvec.h>

#include "utility/python_utility.h"
#include "mesh/face_t.h"
#include "partition/mesh_partition/01_mesh_partition.h"
#include "spatial_discretization/dirichlet_boundary_conditions.h"

// write or load various checkpoints, this is for debugging to only run part of the algorithm on prescribed data
//#define USE_CHECKPOINT_BORDER_POINTS
//#define USE_CHECKPOINT_MESH
//#define WRITE_CHECKPOINT_MESH
#define WRITE_CHECKPOINT_BORDER_POINTS
//#define WRITE_CHECKPOINT_GHOST_MESH
//#define USE_CHECKPOINT_GHOST_MESH

// output STl files for debugging
#define STL_OUTPUT                // output some stl files
//#define STL_OUTPUT_VERBOSE     // output more stl files

//#define FILE_COMMUNICATION        // when sending border points between ranks, do not use MPI but file I/O instead

// include files that implement various methods of this class, these make use the previous defines
#include "postprocessing/parallel_fiber_estimation/create_dirichlet_boundary_conditions.tpp"
#include "postprocessing/parallel_fiber_estimation/create_mesh.tpp"
#include "postprocessing/parallel_fiber_estimation/create_seed_points.tpp"
#include "postprocessing/parallel_fiber_estimation/exchange_ghost_values.tpp"
#include "postprocessing/parallel_fiber_estimation/exchange_seed_points.tpp"
#include "postprocessing/parallel_fiber_estimation/fill_border_points.tpp"
#include "postprocessing/parallel_fiber_estimation/fix_incomplete_streamlines.tpp"
#include "postprocessing/parallel_fiber_estimation/output_border_points.tpp"
#include "postprocessing/parallel_fiber_estimation/rearrange_streamline_points.tpp"
#include "postprocessing/parallel_fiber_estimation/refine_border_points.tpp"
#include "postprocessing/parallel_fiber_estimation/sample_at_equidistant_z_points.tpp"
#include "postprocessing/parallel_fiber_estimation/send_receive_border_points.tpp"
#include "postprocessing/parallel_fiber_estimation/trace_result_fibers.tpp"
#include "postprocessing/parallel_fiber_estimation/trace_streamlines.tpp"


namespace Postprocessing
{

template<typename BasisFunctionType>
ParallelFiberEstimation<BasisFunctionType>::
ParallelFiberEstimation(DihuContext context) :
  context_(context["ParallelFiberEstimation"]), problem_(nullptr), data_(context_), specificSettings_(context_.getPythonConfig())
{
  LOG(TRACE) << "ParallelFiberEstimation::ParallelFiberEstimation()";

  outputWriterManager_.initialize(context_, specificSettings_);

  stlFilename_ = specificSettings_.getOptionString("stlFilename", "");
  resultFilename_ = specificSettings_.getOptionString("resultFilename", "fibers.bin");
  bottomZClip_ = specificSettings_.getOptionDouble("bottomZClip", 0);
  topZClip_ = specificSettings_.getOptionDouble("topZClip", 100);
  nBorderPointsZ_ = specificSettings_.getOptionInt("nElementsZPerSubdomain", 13)+1;
  nBorderPointsX_ = specificSettings_.getOptionInt("nElementsXPerSubdomain", 13)+1;
  maxLevel_ = specificSettings_.getOptionInt("maxLevel", 2);
  nFineGridFibers_ = specificSettings_.getOptionInt("nFineGridFibers", 2);
  improveMesh_ = specificSettings_.getOptionBool("improveMesh", true);
  nNodesPerFiber_ = specificSettings_.getOptionInt("nNodesPerFiber", 1000);

  this->lineStepWidth_ = specificSettings_.getOptionDouble("lineStepWidth", 1e-2, PythonUtility::Positive);
  this->maxNIterations_ = specificSettings_.getOptionInt("maxIterations", 100000, PythonUtility::Positive);
  this->useGradientField_ = specificSettings_.getOptionBool("useGradientField", false);

  // ensure nBorderPointsX is odd
  nBorderPointsX_ = 2*int(nBorderPointsX_/2)+1;
  // ensure nBorderPointsZ is odd
  nBorderPointsZ_ = 2*int(nBorderPointsZ_/2)+1;


  // run stl_create_mesh.rings_to_border_points
  // "Standardize every ring to be in counter-clockwise direction and starting with the point with lowest x coordinate, then sample border points"
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

  functionCreateRingSection_ = PyObject_GetAttrString(moduleStlCreateRings_, "create_ring_section");
  assert(functionCreateRingSection_);

  functionCreateRingSectionMesh_ = PyObject_GetAttrString(moduleStlCreateRings_, "create_ring_section_mesh");
  assert(functionCreateRingSectionMesh_);

  functionGetStlMesh_ = PyObject_GetAttrString(moduleStlCreateRings_, "get_stl_mesh");
  assert(functionGetStlMesh_);

  functionCreateRings_ = PyObject_GetAttrString(moduleStlCreateRings_, "create_rings");
  assert(functionCreateRings_);

  functionRingsToBorderPoints_ = PyObject_GetAttrString(moduleStlCreateMesh_, "rings_to_border_points");
  assert(functionRingsToBorderPoints_);

  functionBorderPointLoopsToList_ = PyObject_GetAttrString(moduleStlCreateMesh_, "border_point_loops_to_list");
  assert(functionRingsToBorderPoints_);

  functionOutputPoints_ = PyObject_GetAttrString(moduleStlDebugOutput_, "output_points");
  assert(functionOutputPoints_);

  functionOutputBorderPoints_ = PyObject_GetAttrString(moduleStlDebugOutput_, "output_border_points");
  assert(functionOutputBorderPoints_);

  functionOutputStreamline_ = PyObject_GetAttrString(moduleStlDebugOutput_, "output_streamline");
  assert(functionOutputStreamline_);

  functionOutputStreamlines_ = PyObject_GetAttrString(moduleStlDebugOutput_, "output_streamlines");
  assert(functionOutputStreamlines_);

  functionOutputGhostElements_ = PyObject_GetAttrString(moduleStlDebugOutput_, "output_ghost_elements");
  assert(functionOutputGhostElements_);

  functionCreate3dMeshFromBorderPointsFaces_ = PyObject_GetAttrString(moduleStlCreateMesh_, "create_3d_mesh_from_border_points_faces");
  assert(functionCreate3dMeshFromBorderPointsFaces_);
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
  LOG(DEBUG) << "generateParallelMesh, nBorderPointsX: " << nBorderPointsX_;

  MPI_Barrier(MPI_COMM_WORLD);

  // prepare rankSubset
  //currentRankSubset_ = std::make_shared<Partition::RankSubset>(0);  // rankSubset with only the master rank (rank no 0)
  nRanksPerCoordinateDirection_.fill(1);

  // get loops of whole domain

#ifndef USE_CHECKPOINT_BORDER_POINTS     // normal execution

  std::array<std::vector<std::vector<Vec3>>,4> borderPoints;  // borderPoints[face_t][z-level][pointIndex]
  std::array<bool,4> subdomainIsAtBorder;

  // only rank 0 creates the first border points
  if (this->context_.partitionManager()->rankNoCommWorld() == 0)
  {
    // run python script to generate loops for the whole volume
    // run stl_create_rings.create_rings
    // "Create n_loops rings/loops (slices) on a closed surface, in equidistant z-values between bottom_clip and top_clip"

    PyObject* loopsPy = PyObject_CallFunction(functionCreateRings_, "s f f i O", stlFilename_.c_str(), bottomZClip_, topZClip_, nBorderPointsZ_, Py_False);
    PythonUtility::checkForError();
    assert(loopsPy);

    // run stl_create_mesh.rings_to_border_points
    // "Standardize every ring to be in counter-clockwise direction and starting with the point with lowest x coordinate, then sample border points"

    PyObject* returnValue = PyObject_CallFunction(functionRingsToBorderPoints_, "O i", loopsPy, 4*(nBorderPointsX_-1));
    PythonUtility::checkForError();
    if (returnValue == NULL)
    {
      LOG(FATAL) << "rings_to_border_points did not return a valid value";
    }

    // unpack the return value which is a tuple (borderPoints, lengths)
    PyObject *borderPointsPy = NULL;
    PyObject *lengthsPy = NULL;
    if (PyTuple_Check(returnValue))
    {
      borderPointsPy = PyTuple_GetItem(returnValue, (Py_ssize_t)0);
      lengthsPy = PyTuple_GetItem(returnValue, (Py_ssize_t)1);
    }
    else
    {
      assert(false);
    }

    // run stl_create_mesh.border_point_loops_to_list
    // "transform the points from numpy array to list, such that they can be extracted from the opendihu C++ code"
    PyObject* borderPointLoops = PyObject_CallFunction(functionBorderPointLoopsToList_, "O", borderPointsPy);
    PythonUtility::checkForError();

    std::vector<std::vector<Vec3>> loops = PythonUtility::convertFromPython<std::vector<std::vector<Vec3>>>::get(borderPointLoops);
    std::vector<double> lengths = PythonUtility::convertFromPython<std::vector<double>>::get(lengthsPy);
    //LOG(DEBUG) << "loops: " << loops << " (size: " << loops.size() << ")" << ", lengths: " << lengths;

#ifndef NDEBUG
#ifdef STL_OUTPUT
    // output the loops
    PyObject_CallFunction(functionOutputStreamlines_, "s i O f", "00_loops", currentRankSubset_->ownRankNo(),
                          PythonUtility::convertToPython<std::vector<std::vector<Vec3>>>::get(loops), 0.1);
    PythonUtility::checkForError();
#endif
#endif

    // rearrange the border points from the loops to the portions of the faces
    //std::array<std::vector<std::vector<Vec3>>,4> borderPoints;  // borderPoints[face_t][z-level][pointIndex]

    // loop over faces: face0Minus = 0, face0Plus, face1Minus, face1Plus

    //   ^ --(1+)-> ^
    // ^ 0-         0+
    // | | --(1-)-> |
    // +-->

    for (int face = Mesh::face_t::face0Minus; face <= Mesh::face_t::face1Plus; face++)
    {
      borderPoints[face].resize(nBorderPointsZ_);
      for (int zIndex = 0; zIndex < nBorderPointsZ_; zIndex++)
      {
        borderPoints[face][zIndex].resize(nBorderPointsX_);

        VLOG(1) << "face " << face << ", zIndex " << zIndex << " nloops: " << loops[zIndex].size();

        if (face == Mesh::face_t::face1Minus)
        {
          std::copy(loops[zIndex].begin(), loops[zIndex].begin() + (nBorderPointsX_-1)+1, borderPoints[face][zIndex].begin());
        }
        else if (face == Mesh::face_t::face0Plus)
        {
          std::copy(loops[zIndex].begin() + (nBorderPointsX_-1), loops[zIndex].begin() + 2*(nBorderPointsX_-1)+1, borderPoints[face][zIndex].begin());
        }
        else if (face == Mesh::face_t::face1Plus)
        {
          std::reverse_copy(loops[zIndex].begin() + 2*(nBorderPointsX_-1), loops[zIndex].begin() + 3*(nBorderPointsX_-1)+1, borderPoints[face][zIndex].begin());
        }
        else if (face == Mesh::face_t::face0Minus)
        {
          std::reverse_copy(loops[zIndex].begin() + 3*(nBorderPointsX_-1), loops[zIndex].begin() + 4*(nBorderPointsX_-1), borderPoints[face][zIndex].begin()+1);
          borderPoints[face][zIndex][0] = loops[zIndex].front();
        }
      }
    }

    subdomainIsAtBorder[Mesh::face_t::face0Minus] = true;
    subdomainIsAtBorder[Mesh::face_t::face0Plus] = true;
    subdomainIsAtBorder[Mesh::face_t::face1Minus] = true;
    subdomainIsAtBorder[Mesh::face_t::face1Plus] = true;
  }

#else       // using previously stored data

  std::vector<int> ranks(64);
  std::iota(ranks.begin(), ranks.end(), 0);
  currentRankSubset_ = std::make_shared<Partition::RankSubset>(ranks.begin(), ranks.end());
  nRanksPerCoordinateDirection_.fill(2);


  std::array<bool,4> subdomainIsAtBorder;

  // parse file
  std::array<std::vector<std::vector<Vec3>>,4> borderPointsSubdomain;

  int subdomainIndex = currentRankSubset_->ownRankNo();
  std::stringstream filename;
  filename << "checkpoint_borderPoints_subdomain_" << subdomainIndex << ".csv";
  std::ifstream file(filename.str().c_str(), std::ios::in);
  if (!file.is_open())
  {
    LOG(ERROR) << "Could not open file \"" << filename.str() << "\" for reading'";
  }
  assert(file.is_open());

  char c;
  int size1, size2;
  file >> size1 >> c >> size2 >> c;
  for (int i = 0; i < 4; i++)
  {
    int value;
    file >> value >> c;
    subdomainIsAtBorder[i] = (value == 1);
  }
  std::string line;
  std::getline(file, line);


  LOG(DEBUG) << "parse values from file \"" << filename.str() << "\".";
  LOG(DEBUG) << "size1: " << size1 << ", size2: " << size2;


  for (int faceNo = (int)Mesh::face_t::face0Minus; faceNo <= (int)Mesh::face_t::face1Plus; faceNo++)
  {
    borderPointsSubdomain[faceNo].resize(size1);
    for (int zLevelIndex = 0; zLevelIndex < borderPointsSubdomain[faceNo].size(); zLevelIndex++)
    {
      borderPointsSubdomain[faceNo][zLevelIndex].resize(size2);
      std::string line;
      std::getline(file, line);

      for (int pointIndex = 0; pointIndex < borderPointsSubdomain[faceNo][zLevelIndex].size(); pointIndex++)
      {
        for (int i = 0; i < 3; i++)
        {
          borderPointsSubdomain[faceNo][zLevelIndex][pointIndex][i] = atof(line.substr(0, line.find(";")).c_str());
          line = line.substr(line.find(";")+1);
        }
      }
    }
  }

  file.close();

  LOG(DEBUG) << "subdomainIndex: " << subdomainIndex;
  LOG(DEBUG) << "subdomainIsAtBorder: " << subdomainIsAtBorder;
  // LOG(DEBUG) << "borderPointsSubdomain: " << borderPointsSubdomain;

  std::array<std::vector<std::vector<Vec3>>,4> &borderPoints = borderPointsSubdomain; // borderPoints[face_t][z-level][pointIndex]

  MPI_Barrier(MPI_COMM_WORLD);
#endif

  generateParallelMeshRecursion(borderPoints, subdomainIsAtBorder);
}

template<typename BasisFunctionType>
void ParallelFiberEstimation<BasisFunctionType>::
generateParallelMeshRecursion(std::array<std::vector<std::vector<Vec3>>,4> &borderPointsOld, std::array<bool,4> subdomainIsAtBorder)
{
  LOG(DEBUG) << "generateParallelMeshRecursion, n border points: "
    << borderPointsOld[0].size() << "," << borderPointsOld[1].size() << "," << borderPointsOld[2].size() << "," << borderPointsOld[3].size();

#ifndef NDEBUG
#ifdef STL_OUTPUT
//#ifdef STL_OUTPUT_VERBOSE
  if (currentRankSubset_->ownRankIsContained())
  {
    PyObject_CallFunction(functionOutputBorderPoints_, "s i O f", "01_border_points_old", currentRankSubset_->ownRankNo(),
                          PythonUtility::convertToPython<std::array<std::vector<std::vector<Vec3>>,4>>::get(borderPointsOld), 0.2);
    PythonUtility::checkForError();
  }

#endif
//#endif
#endif

  // check if all vectors are filled with points
  for (int face = (int)Mesh::face_t::face0Minus; face <= (int)Mesh::face_t::face1Plus; face++)
  {
    for (int zLevelIndex = 0; zLevelIndex < borderPointsOld[face].size(); zLevelIndex++)
    {
      if (borderPointsOld[face][zLevelIndex].empty())
      {
        LOG(ERROR) << "borderPoints face: " << Mesh::getString((Mesh::face_t)face) << " z: " << zLevelIndex << " contains no points.";
      }
    }
  }

  //! new border points for the next subdomain, these are computed on the local subdomain and then send to the sub-subdomains
  std::array<std::array<std::vector<std::vector<Vec3>>,4>,8> borderPointsSubdomain;  // [subdomain index][face_t][z-level][point index]

  // check if the algorithm is at the stage where no more subdomains are created and the final fibers are traced
  // level_ = log2(nRanksPerCoordinateDirection_)
  bool traceFinalFibers = checkTraceFinalFibers();
  LOG(DEBUG) << "level " << level_;

  bool refineSubdomainsOnThisRank = false;

  if (!currentRankSubset_->ownRankIsContained())
  {
    LOG(DEBUG) << "own rank is not contained in currentRankSubset_";
  }

  // here we have the border points of our subdomain, now the borders that split this subdomain in 8 more are computed
  if (currentRankSubset_->ownRankIsContained())
  {
    assert(!borderPointsOld[0].empty());

    refineSubdomainsOnThisRank = true;
    nBorderPointsXNew_ = nBorderPointsX_*2-1;   // number of border points per coordinate direction of the logical rectangle, for the subdomain

    // refine the given border points in x and y direction, i.e. horizontally, to yield twice the number of points on each border
    std::array<std::vector<std::vector<Vec3>>,4> borderPoints;    // [face_t][z-level][pointIndex]
    refineBorderPoints(borderPointsOld, borderPoints);

    // dimensions of borderPoints[face_t][z-level][pointIndex]: borderPoints[4][nBorderPointsZ_][nBorderPointsXNew_]

    // create 3D mesh in own subdomain, using python, algorithm with harmonic maps
    std::vector<Vec3> nodePositions;
    std::array<int,3> nElementsPerCoordinateDirectionLocal;
    createMesh(borderPoints, nodePositions, nElementsPerCoordinateDirectionLocal);

    std::array<global_no_t,3> nElementsPerCoordinateDirectionGlobal;

    // create meshPartition for function space with the currently used ranks
    context_.partitionManager()->setRankSubsetForNextCreatedPartitioning(currentRankSubset_);

    //std::shared_ptr<Partition::MeshPartition<FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<3>, BasisFunctionType>>> meshPartition
    meshPartition_ = context_.partitionManager()->template createPartitioningStructuredLocal<FunctionSpaceType>(
        nElementsPerCoordinateDirectionGlobal, nElementsPerCoordinateDirectionLocal, nRanksPerCoordinateDirection_);

    LOG(DEBUG) << "nElementsPerCoordinateDirectionGlobal: " << nElementsPerCoordinateDirectionGlobal;

    // create function space
    std::stringstream meshName;
    meshName << "meshLevel" << level_;

    int nNodePositionsWithoutGhostsX = nBorderPointsXNew_ - 1;
    int nNodePositionsWithoutGhostsY = nBorderPointsXNew_ - 1;
    int nNodePositionsWithoutGhostsZ = nBorderPointsZ_ - 1;

    if (meshPartition_->hasFullNumberOfNodes(0))
      nNodePositionsWithoutGhostsX += 1;
    if (meshPartition_->hasFullNumberOfNodes(1))
      nNodePositionsWithoutGhostsY += 1;
    if (meshPartition_->hasFullNumberOfNodes(2))
      nNodePositionsWithoutGhostsZ += 1;

    int nNodePositionsWithoutGhosts = nNodePositionsWithoutGhostsX * nNodePositionsWithoutGhostsY * nNodePositionsWithoutGhostsZ;
    std::vector<Vec3> nodePositionsWithoutGhosts(nNodePositionsWithoutGhosts);

    LOG(DEBUG) << "nNodePositionsWithoutGhosts: " << nNodePositionsWithoutGhosts << ", nNodePositionsWithGhosts: " << nodePositions.size();
    LOG(DEBUG) << nNodePositionsWithoutGhostsX << "x" << nNodePositionsWithoutGhostsY << "x" << nNodePositionsWithoutGhostsZ;
    LOG(DEBUG) << nBorderPointsXNew_ << "x" << nBorderPointsXNew_ << "x" << nBorderPointsZ_ << "=" << nBorderPointsZ_*nBorderPointsXNew_*nBorderPointsXNew_;

    for (dof_no_t z = 0; z < nBorderPointsZ_; z++)
    {
      for (dof_no_t y = 0; y < nBorderPointsXNew_; y++)
      {
        for (dof_no_t x = 0; x < nBorderPointsXNew_; x++)
        {
          dof_no_t indexWithoutGhosts = z*nNodePositionsWithoutGhostsX*nNodePositionsWithoutGhostsY + y*nNodePositionsWithoutGhostsX + x;
          dof_no_t indexWithGhosts = z*nBorderPointsXNew_*nBorderPointsXNew_ + y*nBorderPointsXNew_ + x;

          if (z < nNodePositionsWithoutGhostsZ && y < nNodePositionsWithoutGhostsY && x < nNodePositionsWithoutGhostsX)
          {
            assert(indexWithGhosts < nodePositions.size());
            assert(indexWithoutGhosts < nodePositionsWithoutGhosts.size());

            nodePositionsWithoutGhosts[indexWithoutGhosts] = nodePositions[indexWithGhosts];
          }
        }
      }
    }

    context_.partitionManager()->setRankSubsetForNextCreatedPartitioning(currentRankSubset_);
    this->functionSpace_ = context_.meshManager()->template createFunctionSpaceWithGivenMeshPartition<FunctionSpaceType>(
      meshName.str(), meshPartition_, nodePositionsWithoutGhosts, nElementsPerCoordinateDirectionLocal, nRanksPerCoordinateDirection_);

    LOG(DEBUG) << "n nodePositions with ghosts: " << nodePositions.size();
    LOG(DEBUG) << "n nodePositions without ghosts: " << nodePositionsWithoutGhosts.size();
    LOG(DEBUG) << "n dofs local without ghosts: " << this->functionSpace_->meshPartition()->nDofsLocalWithoutGhosts();
    LOG(DEBUG) << "nElementsPerCoordinateDirectionLocal: " << nElementsPerCoordinateDirectionLocal;

    //LOG(DEBUG) << "nodePositions: " << nodePositions;
    std::vector<Vec3> geometryFieldValues;
    this->functionSpace_->geometryField().getValuesWithoutGhosts(geometryFieldValues);
    LOG(DEBUG) << "geometryFieldValues: " << geometryFieldValues;
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

    PyObject_CallFunction(functionOutputPoints_, "s i O f", "03_geometry_values", currentRankSubset_->ownRankNo(),
                          PythonUtility::convertToPython<std::vector<Vec3>>::get(geometryValues), 0.1);
    PythonUtility::checkForError();
#endif
#endif

    // create dirichlet boundary condition object
    std::shared_ptr<SpatialDiscretization::DirichletBoundaryConditions<FunctionSpaceType,1>> dirichletBoundaryConditions;
    createDirichletBoundaryConditions(nElementsPerCoordinateDirectionLocal, dirichletBoundaryConditions);

    // initialize data object to allocate gradient field
    data_.setFunctionSpace(this->functionSpace_);
    data_.reset();
    data_.initialize();

    // set dirichlet values in data for debugging output of vector of dirichlet BC values
    const std::vector<dof_no_t> &dirichletDofs = dirichletBoundaryConditions->boundaryConditionNonGhostDofLocalNos();
    const std::vector<VecD<1>> &dirichletValues = dirichletBoundaryConditions->boundaryConditionValues();
    data_.dirichletValues()->setValues(-1);
    data_.dirichletValues()->setValues(dirichletDofs, dirichletValues);

    // set boundary conditions to the problem
    problem_->setDirichletBoundaryConditions(dirichletBoundaryConditions);
    problem_->reset();
    problem_->initialize();

    // solve the laplace problem, globally on all subdomains on all ranks of the current rank subset
    problem_->run();

    // compute a gradient field from the solution
    problem_->data().solution()->computeGradientField(data_.gradient());

    data_.gradient()->setRepresentationGlobal();

    // output the results
    this->outputWriterManager_.writeOutput(data_);  // output gradient

    //if (level > 1)
    //  LOG(FATAL) << "end";

    // exchange layer of ghost elements (not nodes) between neighbouring ranks,
    // then every subdomain knows one layer of elements around it (only in x/y direction)
    LOG(DEBUG) << "\nConstruct ghost elements";

    exchangeGhostValues(subdomainIsAtBorder);

    // initialize values in base class
    this->solution_ = problem_->data().solution();
    this->gradient_ = data_.gradient();

    // determine seed points for streamlines
    // nodePositions contains all node positions in the current 3D mesh
    // naming of borders: (0-,1-,0+,1+)
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
      seedPointsZIndex = nBorderPointsZ_-1;
      streamlineDirection = -1.0;
    }

    if (nRanksZ == 1)
    {
      // if there is only one rank, start the streamlines at the center and trace towards top and bottom
      seedPointsZIndex = nBorderPointsZ_/2;
    }

    // determine the seed points of the streamlines
    createSeedPoints(subdomainIsAtBorder, seedPointsZIndex, nodePositions, seedPoints);

    // trace streamlines from seed points
    std::vector<std::vector<Vec3>> streamlinePoints;
    traceStreamlines(nRanksZ, rankZNo, streamlineDirection, streamlineDirectionUpwards, seedPoints, streamlinePoints);

    // sample streamlines at equidistant z points
    nBorderPointsZNew_ = nBorderPointsZ_*2 - 1;
    std::vector<std::vector<Vec3>> streamlineZPoints;
    sampleAtEquidistantZPoints(streamlinePoints, seedPoints, streamlineZPoints);

    // save streamline points to the portions of the faces
    //std::array<std::array<std::vector<std::vector<Vec3>>,4>,8> borderPointsSubdomain;  // [subdomain index][face_t][z-level][point index]

    // assign sampled points to the data structure borderPointsSubdomain, which contains the points for each subdomain and face, as list of points for each z level
    std::array<std::array<std::vector<bool>,4>,8> borderPointsSubdomainAreValid;
    std::vector<std::vector<Vec3>> cornerStreamlines;
    rearrangeStreamlinePoints(streamlineZPoints, borderPointsSubdomain, cornerStreamlines, borderPointsSubdomainAreValid, subdomainIsAtBorder);

    LOG(DEBUG) << "nBorderPointsZ_: " << nBorderPointsZ_ << ", nBorderPointsZNew_: " << nBorderPointsZNew_;

    // write border points to file
    outputStreamlines(borderPointsSubdomain, "09_traced");

#if 1
    // fill the streamline points that are at the boundary
    fillBorderPoints(borderPoints, borderPointsSubdomain, cornerStreamlines, borderPointsSubdomainAreValid, subdomainIsAtBorder);

    MPI_Barrier(this->currentRankSubset_->mpiCommunicator());
    // write border points to file
    outputStreamlines(borderPointsSubdomain, "10_filled");
#endif
    MPI_Barrier(this->currentRankSubset_->mpiCommunicator());

    // fill in missing points
    fixIncompleteStreamlines(borderPointsSubdomain, borderPointsSubdomainAreValid, streamlineDirectionUpwards, subdomainIsAtBorder, borderPoints);

    // write border points to file
    outputStreamlines(borderPointsSubdomain, "11_fixed");

    // if this is the end of the recursion, trace final fibers
    if (traceFinalFibers)
    {
      traceResultFibers(streamlineDirection, seedPointsZIndex, nodePositions, borderPointsSubdomain);

      MPI_Barrier(this->currentRankSubset_->mpiCommunicator());
      LOG(DEBUG) << "algorithm is finished";

      // algorithm is finished
      return;
    }
  }  // if own rank is part of this stage of the algorithm

  LOG(DEBUG) << "create subdomains";

  // create subdomains
  // create new rank subset i.e. a new MPI communicator
  int nRanksPerCoordinateDirectionPreviously = nRanksPerCoordinateDirection_[0];
  nRanksPerCoordinateDirection_.fill(nRanksPerCoordinateDirectionPreviously*2);

  int nRanks = nRanksPerCoordinateDirection_[0] * nRanksPerCoordinateDirection_[1] * nRanksPerCoordinateDirection_[2];
  std::vector<int> ranks(nRanks);
  std::iota(ranks.begin(), ranks.end(), 0);
  currentRankSubset_ = std::make_shared<Partition::RankSubset>(ranks.begin(), ranks.end());

  LOG(DEBUG) << "refineSubdomainsOnThisRank: " << refineSubdomainsOnThisRank << ", rankSubset: " << *currentRankSubset_;

  // send border points to ranks that will handle the new subdomains
  std::vector<MPI_Request> sendRequests;
  std::vector<std::vector<double>> sendBuffers;
  if (refineSubdomainsOnThisRank)     // if this rank has created new subdomains (was part of the previous rank subset)
  {
    sendBorderPoints(borderPointsSubdomain, sendBuffers, sendRequests);
  }

  // receive border points
  std::array<std::vector<std::vector<Vec3>>,4> borderPointsNew;
  std::array<bool,4> subdomainIsAtBorderNew;

  if (currentRankSubset_->ownRankIsContained())
  {
    receiveBorderPoints(nRanksPerCoordinateDirectionPreviously, borderPointsNew, subdomainIsAtBorderNew);
  }

#ifndef FILE_COMMUNICATION
  MPIUtility::handleReturnValue(MPI_Waitall(sendRequests.size(), sendRequests.data(), MPI_STATUSES_IGNORE), "MPI_Waitall");
#endif

  sendBuffers.clear();
  //LOG(FATAL) << "done";

  // call method recursively
  generateParallelMeshRecursion(borderPointsNew, subdomainIsAtBorderNew);


  // ----------------------------
  // algorithm:
  // create mesh in own domain, using python, harmonic maps

  // solve laplace problem globally

  // communicate ghost elements to neighbouring subdomains

  // trace fibers to determine new subdomain boundaries, also on the outside (shared between domains)

  // create subdomains
    // create new communicator
    // communicate all old elements to the processes of the new communcator
    // on the new processes create new meshes using the coarse data

  // call method recursively

/*
  LOG(FATAL) << "SUCCESS";*/
}

} // namespace
