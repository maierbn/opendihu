#include "postprocessing/parallel_fiber_estimation.h"

#include <algorithm>
#include <petscvec.h>

#include "utility/python_utility.h"
#include "mesh/face_t.h"
#include "partition/mesh_partition/01_mesh_partition.h"

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
  bottomZClip_ = specificSettings_.getOptionDouble("bottomZClip", 0);
  topZClip_ = specificSettings_.getOptionDouble("topZClip", 100);
  nBorderPointsZ_ = specificSettings_.getOptionInt("nElementsZPerSubdomain", 13)+1;
  nBorderPointsX_ = specificSettings_.getOptionInt("nElementsXPerSubdomain", 13)+1;
  maxLevel_ = specificSettings_.getOptionInt("maxLevel", 2);

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

  // prepare rankSubset
  currentRankSubset_ = std::make_shared<Partition::RankSubset>(0);  // rankSubset with only the master rank (rank no 0)
  nRanksPerCoordinateDirection_.fill(1);

  // get loops of whole domain

  // run python script to generate loops for the whole volume
  // run stl_create_rings.create_rings
  // "Create n_loops rings/loops (slices) on a closed surface, in equidistant z-values between bottom_clip and top_clip"


  PyObject* functionCreateRings = PyObject_GetAttrString(moduleStlCreateRings_, "create_rings");
  assert(functionCreateRings);

  PyObject* loopsPy = PyObject_CallFunction(functionCreateRings, "s f f i O", stlFilename_.c_str(), bottomZClip_, topZClip_, nBorderPointsZ_, Py_False);
  PythonUtility::checkForError();
  assert(loopsPy);

  // run stl_create_mesh.rings_to_border_points
  // "Standardize every ring to be in counter-clockwise direction and starting with the point with lowest x coordinate, then sample border points"

  PyObject* functionRingsToBorderPoints = PyObject_GetAttrString(moduleStlCreateMesh_, "rings_to_border_points");
  assert(functionRingsToBorderPoints);

  PyObject* returnValue = PyObject_CallFunction(functionRingsToBorderPoints, "O i", loopsPy, 4*(nBorderPointsX_-1));
  PythonUtility::checkForError();
  if (returnValue == NULL)
  {
    LOG(FATAL) << "rings_to_border_points did not return a valid value";
  }

  // unpack the return value which is a tuple (borderPoints, lengths)
  PyObject *borderPointsPy;
  PyObject *lengthsPy;
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
  PyObject* functionBorderPointLoopsToList = PyObject_GetAttrString(moduleStlCreateMesh_, "border_point_loops_to_list");
  assert(functionRingsToBorderPoints);

  PyObject* borderPointLoops = PyObject_CallFunction(functionBorderPointLoopsToList, "O", borderPointsPy);
  PythonUtility::checkForError();

  std::vector<std::vector<Vec3>> loops = PythonUtility::convertFromPython<std::vector<std::vector<Vec3>>>::get(borderPointLoops);
  std::vector<double> lengths = PythonUtility::convertFromPython<std::vector<double>>::get(lengthsPy);
  //LOG(DEBUG) << "loops: " << loops << " (size: " << loops.size() << ")" << ", lengths: " << lengths;

#ifndef NDEBUG
  // output the loops
  PyObject* functionOutputPoints = PyObject_GetAttrString(moduleStlDebugOutput_, "output_points");
  assert(functionOutputPoints);

  std::vector<Vec3> outputPoints;
  for (int loopIndex = 0; loopIndex < loops.size(); loopIndex++)
  {
    for (int pointIndex = 0; pointIndex < loops[loopIndex].size(); pointIndex++)
    {
      outputPoints.push_back(loops[loopIndex][pointIndex]);
    }
  }
  PyObject_CallFunction(functionOutputPoints, "s i O f", "00_loops", currentRankSubset_->ownRankNo(),
                        PythonUtility::convertToPython<std::vector<Vec3>>::get(outputPoints), 3.0);
  PythonUtility::checkForError();
#endif

  // rearrange the border points from the loops to the portions of the faces
  std::array<std::vector<std::vector<Vec3>>,4> borderPoints;  // borderPoints[face_t][z-level][pointIndex]

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

  std::array<bool,4> subdomainIsAtBorder;
  subdomainIsAtBorder[Mesh::face_t::face0Minus] = true;
  subdomainIsAtBorder[Mesh::face_t::face0Plus] = true;
  subdomainIsAtBorder[Mesh::face_t::face1Minus] = true;
  subdomainIsAtBorder[Mesh::face_t::face1Plus] = true;

  generateParallelMeshRecursion(borderPoints, subdomainIsAtBorder);
}

template<typename BasisFunctionType>
void ParallelFiberEstimation<BasisFunctionType>::
generateParallelMeshRecursion(std::array<std::vector<std::vector<Vec3>>,4> &borderPointsOld, std::array<bool,4> subdomainIsAtBorder)
{
  LOG(DEBUG) << "generateParallelMeshRecursion";

#ifndef NDEBUG
  PyObject* functionOutputBorderPoints = PyObject_GetAttrString(moduleStlDebugOutput_, "output_border_points");
  assert(functionOutputBorderPoints);

  PyObject_CallFunction(functionOutputBorderPoints, "s i O f", "01_border_points_old", currentRankSubset_->ownRankNo(),
                        PythonUtility::convertToPython<std::array<std::vector<std::vector<Vec3>>,4>>::get(borderPointsOld), 2.0);
  PythonUtility::checkForError();
#endif

  //! new border points for the next subdomain, these are computed on the local subdomain and then send to the sub-subdomains
  std::array<std::array<std::vector<std::vector<Vec3>>,4>,8> borderPointsSubdomain;  // [subdomain index][face_t][z-level][point index]
  std::shared_ptr<Partition::MeshPartition<FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<3>, BasisFunctionType>>> meshPartition;  // the mesh partition of this subdomain which contains information about the neighbouring ranks and the own index in the ranks grid
  bool refineSubdomainsOnThisRank = false;

  // determine current level = log2(nRanksPerCoordinateDirection_)
  int level = 0;
  int nRanksPerCoordinateDirection = nRanksPerCoordinateDirection_[0];
  while (nRanksPerCoordinateDirection >>= 1)
  {
    level++;
  }

  // here we have the border points of our subdomain, now the borders that split this subdomain in 8 more are computed
  if (currentRankSubset_->ownRankIsContained())
  {
    refineSubdomainsOnThisRank = true;
    int nBorderPointsXNew = nBorderPointsX_*2-1;

    // refine borderPoints to twice the precision, only in x and y direction, stays the same in z direction
    // then we have nBorderPointsXNew points per x,y-direction and nBorderPointsZ_ in z direction, each time also including first and list border point
    // borderPoints[face_t][z-level][pointIndex]
    std::array<std::vector<std::vector<Vec3>>,4> borderPoints;
    for (int face = Mesh::face_t::face0Minus; face <= Mesh::face_t::face1Plus; face++)
    {

      borderPoints[face].resize(nBorderPointsZ_);
      for (int zLevelIndex = 0; zLevelIndex < nBorderPointsZ_; zLevelIndex++)
      {
        borderPoints[face][zLevelIndex].resize(nBorderPointsXNew);
        Vec3 previousPoint = borderPointsOld[face][zLevelIndex][0];

        borderPoints[face][zLevelIndex][0] = previousPoint;

        for (int pointIndex = 1; pointIndex < nBorderPointsX_; pointIndex++)
        {
          const Vec3 &currentPoint = borderPointsOld[face][zLevelIndex][pointIndex];
          borderPoints[face][zLevelIndex][2*pointIndex-1] = 0.5*(currentPoint + previousPoint);

          // interpolate point in between
          borderPoints[face][zLevelIndex][2*pointIndex+0] = currentPoint;

          previousPoint = currentPoint;
        }
      }
    }

#ifndef NDEBUG
    PyObject_CallFunction(functionOutputBorderPoints, "s i O f", "02_border_points", currentRankSubset_->ownRankNo(),
                          PythonUtility::convertToPython<std::array<std::vector<std::vector<Vec3>>,4>>::get(borderPoints), 1.0);
    PythonUtility::checkForError();
#endif

    // dimensions of borderPoints[face_t][z-level][pointIndex]: borderPoints[4][nBorderPointsZ_][nBorderPointsXNew]

    // create mesh in own domain, using python, harmonic maps

    // call stl_create_mesh.create_3d_mesh_from_border_points_faces
    PyObject *functionCreate3dMeshFromBorderPointsFaces = PyObject_GetAttrString(moduleStlCreateMesh_, "create_3d_mesh_from_border_points_faces");
    assert(functionCreate3dMeshFromBorderPointsFaces);
    PythonUtility::checkForError();

    PyObject *borderPointsFacesPy = PythonUtility::convertToPython<std::array<std::vector<std::vector<Vec3>>,4>>::get(borderPoints);
    PythonUtility::checkForError();

    //LOG(DEBUG) << PythonUtility::getString(borderPointsFacesPy);
    LOG(DEBUG) << "call function create_3d_mesh_from_border_points_faces";

    PyObject *meshData = PyObject_CallFunction(functionCreate3dMeshFromBorderPointsFaces, "(O)", borderPointsFacesPy);
    PythonUtility::checkForError();

    //LOG(DEBUG) << PythonUtility::getString(meshData);
    // return value:
    //data = {
    //  "node_positions": node_positions,
    //  "linear_elements": linear_elements,
    //  "quadratic_elements": quadratic_elements,
    //  "seed_points": seed_points,
    //  "bottom_nodes": bottom_node_indices,
    //  "top_nodes": top_node_indices,
    //  "n_linear_elements_per_coordinate_direction": n_linear_elements_per_coordinate_direction,
    //  "n_quadratic_elements_per_coordinate_direction": n_quadratic_elements_per_coordinate_direction,
    //}

    std::vector<Vec3> nodePositions;
    PyObject *object = PythonUtility::getOptionPyObject(meshData, "node_positions", "");
    nodePositions = PythonUtility::convertFromPython<std::vector<Vec3>>::get(object);
    std::array<int,3> nElementsPerCoordinateDirectionLocal = PythonUtility::getOptionArray<int,3>(meshData, "n_linear_elements_per_coordinate_direction", "", std::array<int,3>({0,0,0}));

    LOG(DEBUG) << "nodePositions: " << nodePositions;
    LOG(DEBUG) << "nElementsPerCoordinateDirectionLocal: " << nElementsPerCoordinateDirectionLocal;

    std::array<global_no_t,3> nElementsPerCoordinateDirectionGlobal;

    // create meshPartition for function space
    //std::shared_ptr<Partition::MeshPartition<FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<3>, BasisFunctionType>>> meshPartition
    meshPartition = context_.partitionManager()->template createPartitioningStructuredLocal<FunctionSpaceType>(
        nElementsPerCoordinateDirectionGlobal, nElementsPerCoordinateDirectionLocal, nRanksPerCoordinateDirection_);

    LOG(DEBUG) << "nElementsPerCoordinateDirectionGlobal: " << nElementsPerCoordinateDirectionGlobal;

    // create function space
    std::stringstream meshName;
    meshName << "meshLevel" << level;

    context_.partitionManager()->setRankSubsetForNextCreatedMesh(currentRankSubset_);
    this->functionSpace_ = context_.meshManager()->template createFunctionSpaceWithGivenMeshPartition<FunctionSpaceType>(
      meshName.str(), meshPartition, nodePositions, nElementsPerCoordinateDirectionLocal);

    /*
    std::array<global_no_t,FunctionSpace::dim()> &nElementsGlobal,
      const std::array<element_no_t,FunctionSpace::dim()> nElementsLocal,
      const std::array<int,FunctionSpace::dim()> nRanks
    * */
    // solve laplace problem globally
    // create problem
    problem_ = std::make_shared<FiniteElementMethodType>(context_, this->functionSpace_);

    // create dirichlet boundary condition object
    std::shared_ptr<SpatialDiscretization::DirichletBoundaryConditions<FunctionSpaceType,1>> dirichletBoundaryConditions
      = std::make_shared<SpatialDiscretization::DirichletBoundaryConditions<FunctionSpaceType,1>>();

    typedef typename SpatialDiscretization::DirichletBoundaryConditions<FunctionSpaceType,1>::ElementWithNodes ElementWithNodes;

    std::vector<ElementWithNodes> boundaryConditionElements;
    std::vector<dof_no_t> boundaryConditionNonGhostDofLocalNos;
    std::vector<std::array<double,1>> boundaryConditionValues;

    // fill dirichlet boundary condition object
    // set bottom nodes to 0
    std::set<dof_no_t> boundaryConditionNonGhostDofLocalNosSet;
    const int nDofsPerElement1D = FunctionSpace::FunctionSpaceBaseDim<1,BasisFunctionType>::nDofsPerElement();

    // loop over bottom elements
    for (element_no_t elementIndexX = 0; elementIndexX < nElementsPerCoordinateDirectionLocal[0]; elementIndexX++)
    {
      for (element_no_t elementIndexY = 0; elementIndexY < nElementsPerCoordinateDirectionLocal[1]; elementIndexY++)
      {
        element_no_t elementNoLocal = elementIndexY*nElementsPerCoordinateDirectionLocal[1] + elementIndexX;

        ElementWithNodes elementWithNodes;
        elementWithNodes.elementNoLocal = elementNoLocal;

        // loop over dofs of element that are at bottom
        for (dof_no_t elementalDofIndexX = 0; elementalDofIndexX < nDofsPerElement1D; elementalDofIndexX++)
        {
          for (dof_no_t elementalDofIndexY = 0; elementalDofIndexY < nDofsPerElement1D; elementalDofIndexY++)
          {
            dof_no_t elementalDofIndex = elementalDofIndexY*nDofsPerElement1D + elementalDofIndexX;
            elementWithNodes.elementalDofIndex.push_back(std::pair<int,std::array<double,1>>(elementalDofIndex, std::array<double,1>({0.0})));

            dof_no_t dofLocalNo = this->functionSpace_->getDofNo(elementNoLocal, elementalDofIndex);
            boundaryConditionNonGhostDofLocalNosSet.insert(dofLocalNo);
          }
        }
        boundaryConditionElements.push_back(elementWithNodes);
      }
    }

    // create the vectors for dofs and values

    // transfer the data in the set to the vector for the BC indices
    boundaryConditionNonGhostDofLocalNos.resize(boundaryConditionNonGhostDofLocalNosSet.size());
    std::copy(boundaryConditionNonGhostDofLocalNosSet.begin(), boundaryConditionNonGhostDofLocalNosSet.end(), boundaryConditionNonGhostDofLocalNos.begin());

    // set the same amount of values 0.0 for the BC values
    boundaryConditionValues.resize(boundaryConditionNonGhostDofLocalNosSet.size());
    std::fill(boundaryConditionValues.begin(), boundaryConditionValues.end(), std::array<double,1>({0.0}));
    boundaryConditionNonGhostDofLocalNosSet.clear();

    // set top nodes to 1
    // loop over top elements
    element_no_t elementIndexZ = nElementsPerCoordinateDirectionLocal[2]-1;
    for (element_no_t elementIndexX = 0; elementIndexX < nElementsPerCoordinateDirectionLocal[0]; elementIndexX++)
    {
      for (element_no_t elementIndexY = 0; elementIndexY < nElementsPerCoordinateDirectionLocal[1]; elementIndexY++)
      {
        element_no_t elementNoLocal = elementIndexZ*nElementsPerCoordinateDirectionLocal[0]*nElementsPerCoordinateDirectionLocal[1]
          + elementIndexY*nElementsPerCoordinateDirectionLocal[1] + elementIndexX;

        ElementWithNodes elementWithNodes;
        elementWithNodes.elementNoLocal = elementNoLocal;

        // loop over dofs of element that are at top
        dof_no_t elementalDofIndexZ = nDofsPerElement1D-1;
        for (dof_no_t elementalDofIndexX = 0; elementalDofIndexX < nDofsPerElement1D; elementalDofIndexX++)
        {
          for (dof_no_t elementalDofIndexY = 0; elementalDofIndexY < nDofsPerElement1D; elementalDofIndexY++)
          {
            dof_no_t elementalDofIndex = elementalDofIndexZ*nDofsPerElement1D*nDofsPerElement1D
              + elementalDofIndexY*nDofsPerElement1D + elementalDofIndexX;
            elementWithNodes.elementalDofIndex.push_back(std::pair<int,std::array<double,1>>(elementalDofIndex, std::array<double,1>({1.0})));

            dof_no_t dofLocalNo = this->functionSpace_->getDofNo(elementNoLocal, elementalDofIndex);
            boundaryConditionNonGhostDofLocalNosSet.insert(dofLocalNo);
          }
        }
        boundaryConditionElements.push_back(elementWithNodes);
      }
    }

    // fill the vectors for dofs and values
    // transfer the data in the set to the vector for the BC indices
    int nBottomDofs = boundaryConditionNonGhostDofLocalNos.size();
    boundaryConditionNonGhostDofLocalNos.resize(nBottomDofs + boundaryConditionNonGhostDofLocalNosSet.size());
    std::copy(boundaryConditionNonGhostDofLocalNosSet.begin(), boundaryConditionNonGhostDofLocalNosSet.end(), boundaryConditionNonGhostDofLocalNos.begin()+nBottomDofs);

    // add the same amount of 1.0 values for the BC values
    boundaryConditionValues.resize(nBottomDofs + boundaryConditionNonGhostDofLocalNosSet.size());
    std::fill(boundaryConditionValues.begin() + nBottomDofs, boundaryConditionValues.end(), std::array<double,1>({1.0}));

    dirichletBoundaryConditions->initialize(this->functionSpace_, boundaryConditionElements, boundaryConditionNonGhostDofLocalNos, boundaryConditionValues);

    // set boundary conditions to the problem
    problem_->setDirichletBoundaryConditions(dirichletBoundaryConditions);
    problem_->initialize();

    // solve the laplace problem, globally
    problem_->run();

    // initialize data object to allocate gradient field
    data_.setFunctionSpace(this->functionSpace_);
    data_.initialize();

    // compute a gradient field from the solution
    problem_->data().solution()->computeGradientField(data_.gradient());

    // output the results
    this->outputWriterManager_.writeOutput(problem_->data());

    LOG(DEBUG) << "\nConstruct ghost elements";

    struct GhostValues
    {
      std::vector<double> nodePositionValues;   ///< the values of the node positions of the ghost elements, in sequential order like [x,x,x,x, ... y,y,y,y ... z,z,z,...]
      std::vector<double> solutionValues;      ///< values of solution field variable inside the ghost elements
      std::vector<double> gradientValues;      ///< values of the gradient field variable, consecutive like nodePositionValues
      std::array<element_no_t,3> nElementsPerCoordinateDirection;   ///< size of the ghost mesh
    };
    std::array<GhostValues,4> ghostValuesBuffer;  ///< [face], data for meshes containing ghost elements for the sides, face0Minus, face0Plus, face1Minus, face1Plus

    std::vector<MPI_Request> sendRequests, receiveRequests;

    // communicate ghost elements to neighbouring subdomains
    // determine elements on own domain that are to be send to the neighbouring domain
    // loop over faces: face0Minus = 0, face0Plus, face1Minus, face1Plus
    for (int face = Mesh::face_t::face0Minus; face <= Mesh::face_t::face1Plus; face++)
    {
      if (!subdomainIsAtBorder[face])
      {
        LOG(DEBUG) << "make ghost elements for face " << face;

        // get information about neighbouring rank and boundary elements for face
        int neighbourRankNo;
        std::vector<dof_no_t> dofNos;
        meshPartition->getBoundaryElements((Mesh::face_t)face, neighbourRankNo, ghostValuesBuffer[face].nElementsPerCoordinateDirection, dofNos);

        // get relevant values in own domain that will be send to the neighbouring domain
        GhostValues boundaryValues;
        problem_->data().functionSpace()->geometryField().getValues(dofNos, boundaryValues.nodePositionValues);
        problem_->data().solution()->getValues(dofNos, boundaryValues.solutionValues);
        data_.gradient()->getValues(dofNos, boundaryValues.gradientValues);

        int nNodePositionValues = boundaryValues.nodePositionValues.size();
        int nSolutionValues = boundaryValues.solutionValues.size();
        int nGradientValues = boundaryValues.gradientValues.size();
        assert(nSolutionValues*3 == nGradientValues);
        assert(nNodePositionValues == nGradientValues);

        // post receive requests from neighbouring process
        // post non-blocking receive call to receive node position values
        ghostValuesBuffer[face].nodePositionValues.resize(nNodePositionValues);
        receiveRequests.push_back(0);
        MPIUtility::handleReturnValue(MPI_Irecv(ghostValuesBuffer[face].nodePositionValues.data(), nNodePositionValues, MPI_DOUBLE,
                                                neighbourRankNo, 0, currentRankSubset_->mpiCommunicator(), &receiveRequests.back()), "MPI_Irecv");

        // post non-blocking receive call to receive solution values
        ghostValuesBuffer[face].solutionValues.resize(nSolutionValues);
        receiveRequests.push_back(0);
        MPIUtility::handleReturnValue(MPI_Irecv(ghostValuesBuffer[face].solutionValues.data(), nSolutionValues, MPI_DOUBLE,
                                                neighbourRankNo, 0, currentRankSubset_->mpiCommunicator(), &receiveRequests.back()), "MPI_Irecv");

        // post non-blocking receive call to receive gradient values
        ghostValuesBuffer[face].gradientValues.resize(nGradientValues);
        receiveRequests.push_back(0);
        MPIUtility::handleReturnValue(MPI_Irecv(ghostValuesBuffer[face].gradientValues.data(), nGradientValues, MPI_DOUBLE,
                                                neighbourRankNo, 0, currentRankSubset_->mpiCommunicator(), &receiveRequests.back()), "MPI_Irecv");

        // send values to neighbouring process, non-blocking
        // post non-blocking send call to send solution values
        sendRequests.push_back(0);
        MPIUtility::handleReturnValue(MPI_Isend(boundaryValues.nodePositionValues.data(), nNodePositionValues, MPI_DOUBLE,
                                                neighbourRankNo, 0, currentRankSubset_->mpiCommunicator(), &sendRequests.back()), "MPI_Isend");

        // post non-blocking send call to send solution values
        sendRequests.push_back(0);
        MPIUtility::handleReturnValue(MPI_Isend(boundaryValues.solutionValues.data(), nSolutionValues, MPI_DOUBLE,
                                                neighbourRankNo, 0, currentRankSubset_->mpiCommunicator(), &sendRequests.back()), "MPI_Isend");

        // post non-blocking send call to send gradient values
        sendRequests.push_back(0);
        MPIUtility::handleReturnValue(MPI_Isend(boundaryValues.gradientValues.data(), nGradientValues, MPI_DOUBLE,
                                                neighbourRankNo, 0, currentRankSubset_->mpiCommunicator(), &sendRequests.back()), "MPI_Isend");
      }
    }

    // wait for communication to finish
    MPIUtility::handleReturnValue(MPI_Waitall(sendRequests.size(), sendRequests.data(), MPI_STATUSES_IGNORE), "MPI_Waitall");
    MPIUtility::handleReturnValue(MPI_Waitall(receiveRequests.size(), receiveRequests.data(), MPI_STATUSES_IGNORE), "MPI_Waitall");

    // create ghost element meshes
    std::array<std::shared_ptr<FunctionSpaceType>,4> ghostMesh;
    std::array<std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,1>>,4> ghostSolution;
    std::array<std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,3>>,4> ghostGradient;

    for (int face = Mesh::face_t::face0Minus; face <= Mesh::face_t::face1Plus; face++)
    {
      if (!ghostValuesBuffer[face].nodePositionValues.empty())
      {
        std::stringstream meshName;
        meshName << this->functionSpace_->meshName() << "_ghost_" << Mesh::getString((Mesh::face_t)face);

        // transform the node Position data from vector of double to vector of Vec3
        int nNodes = ghostValuesBuffer[face].nodePositionValues.size() / 3;
        std::vector<Vec3> nodePositions(nNodes);

        for (node_no_t nodeIndex = 0; nodeIndex < nNodes; nodeIndex++)
        {
          nodePositions[nodeIndex][0] = ghostValuesBuffer[face].nodePositionValues[nodeIndex*3 + 0];
          nodePositions[nodeIndex][1] = ghostValuesBuffer[face].nodePositionValues[nodeIndex*3 + 1];
          nodePositions[nodeIndex][2] = ghostValuesBuffer[face].nodePositionValues[nodeIndex*3 + 2];
        }

        ghostMesh[face] = context_.meshManager()->template createFunctionSpace<FunctionSpaceType>(meshName.str(), nodePositions, ghostValuesBuffer[face].nElementsPerCoordinateDirection);
        ghostSolution[face] = ghostMesh[face]->template createFieldVariable<1>("solution");
        ghostSolution[face]->setValuesWithGhosts(ghostValuesBuffer[face].solutionValues);
        ghostGradient[face] = ghostMesh[face]->template createFieldVariable<3>("gradient");
        int nGradientValues = ghostValuesBuffer[face].gradientValues.size()/3;
        std::vector<Vec3> gradientValues(nGradientValues);
        for (int i = 0; i < nGradientValues; i++)
        {
          gradientValues[i][0] = ghostValuesBuffer[face].gradientValues[3*i + 0];
          gradientValues[i][1] = ghostValuesBuffer[face].gradientValues[3*i + 1];
          gradientValues[i][2] = ghostValuesBuffer[face].gradientValues[3*i + 2];
        }
        ghostGradient[face]->setValuesWithGhosts(gradientValues);
      }
    }

    // initialize values in base class
    this->functionSpace_ = problem_->data().functionSpace();
    this->solution_ = problem_->data().solution();
    this->gradient_ = data_.gradient();

    // determine seed points
    // nodePositions contains all node positions in the current 3D mesh
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

    int nNodesX = nElementsPerCoordinateDirectionLocal[0]+1;
    int nNodesY = nElementsPerCoordinateDirectionLocal[1]+1;
    int nNodesZ = nElementsPerCoordinateDirectionLocal[2]+1;

    LOG(DEBUG) << "nNodes: " << nNodesX << " x " << nNodesY << " x " << nNodesZ;

    assert(nNodesX == nBorderPointsXNew);
    assert(nNodesY == nBorderPointsXNew);
    assert(nNodesZ == nBorderPointsZ_);

    int rankZNo = meshPartition->ownRankPartitioningIndex(2);
    int nRanksZ = meshPartition->nRanks(2);
    bool streamlineDirectionUpwards = rankZNo > nRanksZ/2;

    int zIndex = 0;
    double streamlineDirection = 1.0;
    if (!streamlineDirectionUpwards)
    {
      zIndex = nNodesZ-1;
      streamlineDirection = -1.0;
    }

    if (nRanksZ == 1)
    {
      // if there is only one rank, start the streamlines at the center and trace towards top and bottom
      zIndex = nNodesZ/2;
    }

    // trace final fibers
    if (level == maxLevel_)
    {
      LOG(DEBUG) << "final level";

      std::vector<Vec3> seedPoints;

      // create seed points
      for (int j = 0; j < nBorderPointsXNew-1; j++)
      {
        for (int i = 0; i < nBorderPointsXNew-1; i++)
        {
          Vec3 p0 = nodePositions[zIndex*nNodesX*nNodesY + j*nNodesX + i];
          Vec3 p1 = nodePositions[zIndex*nNodesX*nNodesY + j*nNodesX + i+1];
          Vec3 p2 = nodePositions[zIndex*nNodesX*nNodesY + (j+1)*nNodesX + i];
          Vec3 p3 = nodePositions[zIndex*nNodesX*nNodesY + (j+1)*nNodesX + i+1];
          seedPoints.push_back(0.25*(p0 + p1 + p2 + p3));
        }
      }

#ifndef NDEBUG
      PyObject* functionOutputPoints = PyObject_GetAttrString(moduleStlDebugOutput_, "output_points");
      assert(functionOutputPoints);

      PyObject_CallFunction(functionOutputPoints, "s i O f", "03_final_seed_points", currentRankSubset_->ownRankNo(),
                            PythonUtility::convertToPython<std::vector<Vec3>>::get(seedPoints), 1.0);
      PythonUtility::checkForError();
#endif

      // trace streamlines from seed points
      int nStreamlines = seedPoints.size();
      std::vector<std::vector<Vec3>> streamlinePoints(nStreamlines);
      for (int i = 0; i < nStreamlines; i++)
      {
        Vec3 &startingPoint = seedPoints[i];
        streamlinePoints[i].push_back(startingPoint);

        for (int face = (int)Mesh::face_t::face0Minus; face <= (int)Mesh::face_t::face1Plus; face++)
        {
          this->functionSpace_->setGhostMesh(Mesh::face_t::face0Minus, ghostMesh[face]);
        }

        this->traceStreamline(startingPoint, streamlineDirection, streamlinePoints[i]);
      }

      // exchange end points of streamlines and adjust them

      // tracing done
      return;
    }

    // boundary indices for face0Minus and face0Plus (vertical direction)
    int iBeginVertical = 0;
    int iEndVertical = nBorderPointsXNew;

    if (subdomainIsAtBorder[(int)Mesh::face_t::face1Minus])
      iBeginVertical += 1;

    if (subdomainIsAtBorder[(int)Mesh::face_t::face1Plus])
      iEndVertical -= 1;

    // boundary indices for face1Minus and face1Plus (horizontal direction)
    int iBeginHorizontal = 0;
    int iEndHorizontal = nBorderPointsXNew;

    if (subdomainIsAtBorder[(int)Mesh::face_t::face0Minus])
      iBeginHorizontal += 1;

    if (subdomainIsAtBorder[(int)Mesh::face_t::face0Plus])
      iEndHorizontal -= 1;

    // face0Minus
    if (!subdomainIsAtBorder[(int)Mesh::face_t::face0Minus])
    {
      for (int i = iBeginVertical; i < iEndVertical; i++)
      {
        seedPoints.push_back(nodePositions[zIndex*nNodesX*nNodesY + i*nNodesX + 0]);
      }
    }

    // face0Plus
    if (!subdomainIsAtBorder[(int)Mesh::face_t::face0Plus])
    {
      for (int i = iBeginVertical; i < iEndVertical; i++)
      {
        seedPoints.push_back(nodePositions[zIndex*nNodesX*nNodesY + i*nNodesX + (nNodesX-1)]);
      }
    }

    // face1Minus (with corner points)
    if (!subdomainIsAtBorder[(int)Mesh::face_t::face1Minus])
    {
      for (int i = iBeginHorizontal; i < iEndHorizontal; i++)
      {
        seedPoints.push_back(nodePositions[zIndex*nNodesX*nNodesY + i]);
      }
    }

    // face1Plus (with corner points)
    if (!subdomainIsAtBorder[(int)Mesh::face_t::face1Plus])
    {
      for (int i = iBeginHorizontal; i < iEndHorizontal; i++)
      {
        seedPoints.push_back(nodePositions[zIndex*nNodesX*nNodesY + (nNodesY-1)*nNodesX + i]);
      }
    }

    // horizontal center line (with corner points)
    for (int i = iBeginHorizontal; i < iEndHorizontal; i++)
    {
      seedPoints.push_back(nodePositions[zIndex*nNodesX*nNodesY + int(nNodesY/2)*nNodesX + i]);
    }

    // vertical center line (with corner points and center point)
    for (int i = iBeginVertical; i < iEndVertical; i++)
    {
      seedPoints.push_back(nodePositions[zIndex*nNodesX*nNodesY + i*nNodesX + int(nNodesX/2)]);
    }

#ifndef NDEBUG
      PyObject* functionOutputPoints = PyObject_GetAttrString(moduleStlDebugOutput_, "output_points");
      assert(functionOutputPoints);

      PyObject_CallFunction(functionOutputPoints, "s i O f", "03_seed_points", currentRankSubset_->ownRankNo(),
                            PythonUtility::convertToPython<std::vector<Vec3>>::get(seedPoints), 0.6);
      PythonUtility::checkForError();
#endif

    // trace streamlines
    int nStreamlines = seedPoints.size();
    std::vector<std::vector<Vec3>> streamlinePoints(nStreamlines);

    if (nRanksZ == 1)
    {
      // if there is only one rank, trace streamline from the center in both directions
      for (int face = (int)Mesh::face_t::face0Minus; face <= (int)Mesh::face_t::face1Plus; face++)
      {
        this->functionSpace_->setGhostMesh(Mesh::face_t::face0Minus, ghostMesh[face]);
      }

      // trace streamlines from seed points
      int nStreamlines = seedPoints.size();
      LOG(DEBUG) << " on 1 rank, trace " << nStreamlines << " streamlines";
      for (int i = 0; i < nStreamlines; i++)
      {
        // get starting point
        Vec3 &startingPoint = seedPoints[i];

        // trace streamlines forwards
        std::vector<Vec3> forwardPoints;

        this->traceStreamline(startingPoint, 1.0, forwardPoints);

        if (forwardPoints.empty())  // if there was not even the first point found
        {
          LOG(ERROR) << "Seed point " << startingPoint << " is outside of domain.";
          continue;
        }

        // trace streamline backwards
        std::vector<Vec3> backwardPoints;
        this->traceStreamline(startingPoint, -1.0, backwardPoints);

        // copy collected points to result vector, note avoiding this additional copy-step is not really possible, since it would require a push_front which is only efficient with lists, but we need a vector here
        streamlinePoints[i].insert(streamlinePoints[i].begin(), backwardPoints.rbegin(), backwardPoints.rend());
        streamlinePoints[i].insert(streamlinePoints[i].end(), startingPoint);
        streamlinePoints[i].insert(streamlinePoints[i].end(), forwardPoints.begin(), forwardPoints.end());

        // now streamline is in order from bottom (Z-) to top (Z+)
        LOG(DEBUG) << "i=" << i << ", n points in streamline: " << streamlinePoints[i].size();
        VLOG(1) << "streamline: " << streamlinePoints[i];

  #ifndef NDEBUG
        std::stringstream name;
        name << "04_streamline_" << i << "_";
        PyObject_CallFunction(functionOutputPoints, "s i O f", name.str().c_str(), currentRankSubset_->ownRankNo(),
                              PythonUtility::convertToPython<std::vector<Vec3>>::get(streamlinePoints[i]), 0.5);
        PythonUtility::checkForError();
  #endif
      }
    }
    else
    {
      // multiple ranks

      // determine if previously set seedPoints are used or if they are received from neighbouring rank
      if (nRanksZ > 1 && rankZNo != int(nRanksZ/2) && rankZNo != int(nRanksZ/2+1))
      {
        int neighbourRankNo;
        if (streamlineDirectionUpwards)
        {
          neighbourRankNo = meshPartition->neighbourRank(Mesh::face_t::face2Minus);
        }
        else
        {
          neighbourRankNo = meshPartition->neighbourRank(Mesh::face_t::face2Plus);
        }

        // receive seed points
        MPIUtility::handleReturnValue(MPI_Recv(seedPoints.data(), seedPoints.size(), MPI_DOUBLE, neighbourRankNo,
                                              0, currentRankSubset_->mpiCommunicator(), MPI_STATUS_IGNORE), "MPI_Recv");
      }

      LOG(DEBUG) << " on " << nRanksZ << " ranks in Z direction, trace " << nStreamlines << " streamlines";

      // trace streamlines from seed points
      for (int i = 0; i < nStreamlines; i++)
      {
        Vec3 &startingPoint = seedPoints[i];
        streamlinePoints[i].push_back(startingPoint);

        for (int face = (int)Mesh::face_t::face0Minus; face <= (int)Mesh::face_t::face1Plus; face++)
        {
          this->functionSpace_->setGhostMesh(Mesh::face_t::face0Minus, ghostMesh[face]);
        }

        this->traceStreamline(startingPoint, streamlineDirection, streamlinePoints[i]);

  #ifndef NDEBUG
        std::stringstream name;
        name << "04_streamline_" << i << "_";
        PyObject_CallFunction(functionOutputPoints, "s i O f", name.str().c_str(), currentRankSubset_->ownRankNo(),
                              PythonUtility::convertToPython<std::vector<Vec3>>::get(streamlinePoints[i]), 1.0);
        PythonUtility::checkForError();
  #endif

      }

      // send end points of streamlines to next rank that continues the streamline
      if (nRanksZ > 1 && rankZNo != nRanksZ-1 && rankZNo != 0)
      {
        // fill send buffer
        std::vector<double> sendBuffer(nStreamlines);
        for (int streamlineIndex = 0; streamlineIndex < nStreamlines; streamlineIndex++)
        {
          for (int i = 0; i < 3; i++)
          {
            sendBuffer[streamlineIndex] = streamlinePoints[streamlineIndex].back()[i];
          }
        }

        int neighbourRankNo;
        if (streamlineDirectionUpwards)
        {
          neighbourRankNo = meshPartition->neighbourRank(Mesh::face_t::face2Plus);
        }
        else
        {
          neighbourRankNo = meshPartition->neighbourRank(Mesh::face_t::face2Minus);
        }

        // send end points of streamlines
        MPIUtility::handleReturnValue(MPI_Send(sendBuffer.data(), sendBuffer.size(), MPI_DOUBLE, neighbourRankNo,
                                              0, currentRankSubset_->mpiCommunicator()), "MPI_Send");
      }

      // reorder streamline points such that they go from bottom to top
      if (streamlineDirection < 0)
      {
        for (int streamlineIndex = 0; streamlineIndex < nStreamlines; streamlineIndex++)
        {
          std::reverse(streamlinePoints[streamlineIndex].begin(), streamlinePoints[streamlineIndex].end());
        }
      }
    }

    // sample streamlines at equidistant z points
    std::vector<std::vector<Vec3>> streamlineZPoints(nStreamlines);

    int nBorderPointsZNew = nBorderPointsZ_*2 - 1;
    double currentZ = bottomZClip_;
    double zIncrement = double(topZClip_ - bottomZClip_) / (nBorderPointsZNew - 1.0);

    LOG(DEBUG) << "z bounds: " << bottomZClip_ << ", " << topZClip_ << ", zIncrement: " << zIncrement;

    for (int i = 0; i < nStreamlines; i++)
    {
      LOG(DEBUG) << " streamline " << i << " has " << streamlinePoints[i].size() << " points.";

      assert(streamlinePoints[i].size() > 1);

      Vec3 previousPoint = streamlinePoints[i].front();
      streamlineZPoints[i].reserve(nBorderPointsZNew);
      streamlineZPoints[i].push_back(previousPoint);
      currentZ = bottomZClip_ + zIncrement;
      LOG(DEBUG) << "first point: " << streamlinePoints[i].front() << ", currentZ: " << currentZ;

      for (std::vector<Vec3>::const_iterator iter = streamlinePoints[i].begin()+1; iter != streamlinePoints[i].end(); iter++)
      {
        const Vec3 &currentPoint = *iter;
        if (currentPoint[2] > currentZ || iter+1 == streamlinePoints[i].end())
        {
          double alpha = (currentZ - previousPoint[2]) / (currentPoint[2] - previousPoint[2]);
          Vec3 point = (1.-alpha) * previousPoint + alpha * currentPoint;
          streamlineZPoints[i].push_back(point);

          currentZ += zIncrement;
        }
        previousPoint = currentPoint;
      }

#ifndef NDEBUG
      std::stringstream name;
      name << "05_sampled_streamline_" << i << "_";
      PyObject_CallFunction(functionOutputPoints, "s i O f", name.str().c_str(), currentRankSubset_->ownRankNo(),
                            PythonUtility::convertToPython<std::vector<Vec3>>::get(streamlineZPoints[i]), 1.0);
      PythonUtility::checkForError();
#endif

      LOG(DEBUG) << "n sampled points: " << streamlineZPoints[i].size() << ", nBorderPointsXNew: " << nBorderPointsXNew << ", nBorderPointsZNew: " << nBorderPointsZNew;

      if (streamlineZPoints[i].size() != nBorderPointsZNew)
      {
        LOG(ERROR) << "Streamline " << i << " is not complete, i.e. does not run from \"bottomZClip\" to \"topZClip\" "
          << std::endl << "Try adjusting \"bottomZClip\" and \"topZClip\" or mesh width.";
      }
      //assert(streamlineZPoints[i].size() == nBorderPointsZNew);
    }


    // save streamline points to the portions of the faces
    //std::array<std::array<std::vector<std::vector<Vec3>>,4>,8> borderPointsSubdomain;  // [subdomain index][face_t][z-level][point index]

    // assign sampled points to the data structure borderPointsSubdomain, which contains the points for each subdomain and face, as list of points for each z level
    reorganizeStreamlinePoints(streamlineZPoints, borderPointsSubdomain, subdomainIsAtBorder);

#ifndef NDEBUG
    PyObject_CallFunction(functionOutputBorderPoints, "s i O f", "06_border_points_subdomain", currentRankSubset_->ownRankNo(),
                          PythonUtility::convertToPython<std::array<std::array<std::vector<std::vector<Vec3>>,4>,8>>::get(borderPointsSubdomain), 1.0);
    PythonUtility::checkForError();
#endif

    // fill points that are on the border of the domain
    // loop over z levels
    for (int zLevelIndex = 0; zLevelIndex < nBorderPointsZNew; zLevelIndex++)
    {
      currentZ = bottomZClip_ + double(zLevelIndex) / (nBorderPointsZNew-1) * (topZClip_ - bottomZClip_);

      Vec3 startPoint;
      Vec3 endPoint;

      for (int face = Mesh::face_t::face0Minus; face <= Mesh::face_t::face1Plus; face++)
      {
        if (subdomainIsAtBorder[face])
        {
          // get start and end point of horizontal line on boundary that is to be found
          if (zLevelIndex % 2 == 0)
          {
            startPoint = borderPoints[face][zLevelIndex/2][0];
            endPoint = borderPoints[face][zLevelIndex/2][nBorderPointsXNew-1];
          }
          else
          {
            startPoint = 0.5 * (borderPoints[face][zLevelIndex/2][0] + borderPoints[face][zLevelIndex/2+1][0]);
            endPoint = 0.5 * (borderPoints[face][zLevelIndex/2][nBorderPointsXNew-1] + borderPoints[face][zLevelIndex/2][nBorderPointsXNew-1]);
          }

        // determine boundary points between startPoint and endPoint, with same zLevel value
        // call stl_create_rings.create_ring_section
        PyObject *functionCreateRingSection = PyObject_GetAttrString(moduleStlCreateRings_, "create_ring_section");
        assert(functionCreateRingSection);

        PyObject *startPointPy = PythonUtility::convertToPython<Vec3>::get(startPoint);
        PyObject *endPointPy = PythonUtility::convertToPython<Vec3>::get(endPoint);
        PyObject *loopSectionPy = PyObject_CallFunction(functionCreateRingSection, "s O O f i", stlFilename_.c_str(), startPointPy, endPointPy, currentZ, nBorderPointsXNew);
        PythonUtility::checkForError();
        //  create_ring_section(input_filename, start_point, end_point, z_value, n_points)
        assert(loopSectionPy);

        std::vector<Vec3> loopSection = PythonUtility::convertFromPython<std::vector<Vec3>>::get(loopSectionPy);

        // determine affected subdomains
        int subdomainIndex0 = 0;
        int subdomainIndex1 = 0;

        if (face == (int)Mesh::face_t::face0Minus)
        {
          subdomainIndex0 = 0;
          subdomainIndex1 = 2;
        }
        else if (face == (int)Mesh::face_t::face0Plus)
        {
          subdomainIndex0 = 1;
          subdomainIndex1 = 3;
        }
        else if (face == (int)Mesh::face_t::face1Minus)
        {
          subdomainIndex0 = 0;
          subdomainIndex1 = 1;
        }
        else if (face == (int)Mesh::face_t::face1Plus)
        {
          subdomainIndex0 = 2;
          subdomainIndex1 = 3;
        }

        // bottom subdomain
        if (zLevelIndex <= nBorderPointsZ_)
        {
          int zLevelIndexSubdomain = zLevelIndex;
          std::copy(loopSection.begin(), loopSection.begin()+nBorderPointsX_, borderPointsSubdomain[subdomainIndex0][face][zLevelIndexSubdomain].begin());
          std::copy(loopSection.begin()+nBorderPointsX_-1, loopSection.end(), borderPointsSubdomain[subdomainIndex1][face][zLevelIndexSubdomain].begin());
        }

        // top subdomain
        if (zLevelIndex >= nBorderPointsZ_)
        {
          int zLevelIndexSubdomain = zLevelIndex - nBorderPointsZ_;
          std::copy(loopSection.begin(), loopSection.begin()+nBorderPointsX_, borderPointsSubdomain[subdomainIndex0][face][zLevelIndexSubdomain].begin());
          std::copy(loopSection.begin()+nBorderPointsX_-1, loopSection.end(), borderPointsSubdomain[subdomainIndex1][face][zLevelIndexSubdomain].begin());
        }
      }
    }


#ifndef NDEBUG
    PyObject_CallFunction(functionOutputBorderPoints, "s i O f", "07_border_points_subdomain_filled", currentRankSubset_->ownRankNo(),
                          PythonUtility::convertToPython<std::array<std::array<std::vector<std::vector<Vec3>>,4>,8>>::get(borderPointsSubdomain), 1.0);
    PythonUtility::checkForError();
#endif

    // output points for debugging
#ifndef NDEBUG
    std::ofstream file("points.csv", std::ios::out | std::ios::trunc);
    assert (file.is_open());
    for (int subdomainIndex = 0; subdomainIndex < 8; subdomainIndex++)
    {
      for (int faceNo = (int)Mesh::face_t::face0Minus; faceNo <= (int)Mesh::face_t::face1Plus; faceNo++)
      {
        for (int zLevelIndex = 0; zLevelIndex < borderPointsSubdomain[subdomainIndex][faceNo].size(); zLevelIndex++)
        {
          for (int pointIndex = 0; pointIndex < borderPointsSubdomain[subdomainIndex][faceNo][zLevelIndex].size(); pointIndex++)
          {
            file << borderPointsSubdomain[subdomainIndex][faceNo][zLevelIndex][0] << ";" << borderPointsSubdomain[subdomainIndex][faceNo][zLevelIndex][1] << ";" << borderPointsSubdomain[subdomainIndex][faceNo][zLevelIndex][2] << ";";
          }
          file << std::endl;
        }
      }
    }
    file.close();
#endif
  }  // if own rank is part of this stage of the algorithm

  // create subdomains
  // create new rank subset i.e. a new MPI communicator
  int nRanksPerCoordinateDirectionPreviously = nRanksPerCoordinateDirection_[0];
  nRanksPerCoordinateDirection_.fill(nRanksPerCoordinateDirectionPreviously*2);

  int nRanks = nRanksPerCoordinateDirection_[0] * nRanksPerCoordinateDirection_[1] * nRanksPerCoordinateDirection_[2];
  std::vector<int> ranks(nRanks);
  std::iota(ranks.begin(), ranks.end(), 0);
  currentRankSubset_ = std::make_shared<Partition::RankSubset>(ranks.begin(), ranks.end());

  // send border points to ranks that will handle the new subdomains
  std::vector<MPI_Request> sendRequest;
  if (refineSubdomainsOnThisRank)
  {
    sendRequest.resize(8);

    // send border points to subdomains
    for (int subdomainIndex = 0; subdomainIndex < 8; subdomainIndex++)
    {
      // fill send buffer
      std::vector<double> sendBuffer;
      for (int faceNo = (int)Mesh::face_t::face0Minus; faceNo <= (int)Mesh::face_t::face1Plus; faceNo++)
      {
        for (int zLevelIndex = 0; zLevelIndex < borderPointsSubdomain[subdomainIndex][faceNo].size(); zLevelIndex++)
        {
          for (int pointIndex = 0; pointIndex < borderPointsSubdomain[subdomainIndex][faceNo][zLevelIndex].size(); pointIndex++)
          {
            for (int i = 0; i < 3; i++)
            {
              sendBuffer.push_back(borderPointsSubdomain[subdomainIndex][faceNo][zLevelIndex][pointIndex][i]);
            }
          }
        }
      }

      std::array<int,3> oldRankIndex;

      for (int i = 0; i < 3; i++)
      {
        oldRankIndex[i] = meshPartition->ownRankPartitioningIndex(i);
      }

      std::array<int,3> subdomainRankIndex;
      subdomainRankIndex[0] = oldRankIndex[0]*2 + subdomainIndex % 2;
      subdomainRankIndex[1] = oldRankIndex[1]*2 + int((subdomainIndex % 4) / 2);
      subdomainRankIndex[2] = oldRankIndex[2]*2 + int(subdomainIndex / 4);

      int subdomainRankNo = subdomainRankIndex[2]*(nRanksPerCoordinateDirection_[0]*nRanksPerCoordinateDirection_[1]) + subdomainRankIndex[1]*nRanksPerCoordinateDirection_[0] + subdomainRankIndex[2];

      MPIUtility::handleReturnValue(MPI_Isend(sendBuffer.data(), sendBuffer.size(), MPI_DOUBLE, subdomainRankNo, 0, currentRankSubset_->mpiCommunicator(), &sendRequest[subdomainIndex]), "MPI_Isend");
    }
  }

  // receive border points
  std::array<std::vector<std::vector<Vec3>>,4> borderPointsNew;
  std::array<bool,4> subdomainIsAtBorderNew;

  if (currentRankSubset_->ownRankIsContained())
  {
    int ownRankNo = currentRankSubset_->ownRankNo();
    std::array<int,3> rankIndex;
    rankIndex[0] = ownRankNo % nRanksPerCoordinateDirection_[0];
    rankIndex[1] = int((ownRankNo % (nRanksPerCoordinateDirection_[0]*nRanksPerCoordinateDirection_[1])) / nRanksPerCoordinateDirection_[0]);
    rankIndex[2] = int(ownRankNo / (nRanksPerCoordinateDirection_[0]*nRanksPerCoordinateDirection_[1]));

    subdomainIsAtBorderNew[(int)Mesh::face_t::face0Minus] = rankIndex[0] == 0;
    subdomainIsAtBorderNew[(int)Mesh::face_t::face0Plus] = rankIndex[0] == nRanksPerCoordinateDirection_[0]-1;
    subdomainIsAtBorderNew[(int)Mesh::face_t::face1Minus] = rankIndex[1] == 0;
    subdomainIsAtBorderNew[(int)Mesh::face_t::face1Plus] = rankIndex[1] == nRanksPerCoordinateDirection_[1]-1;

    int rankToReceiveFrom = int(rankIndex[2]/2)*(nRanksPerCoordinateDirectionPreviously*nRanksPerCoordinateDirectionPreviously) + int(rankIndex[1]/2)*nRanksPerCoordinateDirectionPreviously + int(rankIndex[0]/2);

    std::vector<double> recvBuffer(4*nBorderPointsX_*nBorderPointsX_);
    MPIUtility::handleReturnValue(MPI_Recv(recvBuffer.data(), recvBuffer.size(), MPI_DOUBLE, rankToReceiveFrom, 0, currentRankSubset_->mpiCommunicator(), MPI_STATUS_IGNORE), "MPI_Recv");

    // extract received values
    int recvBufferIndex = 0;
    for (int faceNo = (int)Mesh::face_t::face0Minus; faceNo <= (int)Mesh::face_t::face1Plus; faceNo++)
    {
      borderPointsNew[faceNo].resize(nBorderPointsX_);
      for (int zLevelIndex = 0; zLevelIndex < nBorderPointsX_; zLevelIndex++)
      {
        borderPointsNew[faceNo][zLevelIndex].resize(nBorderPointsX_);
        for (int pointIndex = 0; pointIndex < nBorderPointsX_; pointIndex++, recvBufferIndex+=3)
        {
          borderPointsNew[faceNo][zLevelIndex][pointIndex] = Vec3({recvBuffer[recvBufferIndex], recvBuffer[recvBufferIndex+1], recvBuffer[recvBufferIndex+2]});
        }
      }
    }
  }

  MPIUtility::handleReturnValue(MPI_Waitall(sendRequest.size(), sendRequest.data(), MPI_STATUSES_IGNORE), "MPI_Waitall");

  LOG(FATAL) << "done";

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

    LOG(FATAL) << "SUCCESS";
  /*
    std::vector<PyObject *> loopList = PythonUtility::convertFromPython<std::vector<PyObject *>>::get(returnValue);

    for (std::vector<PyObject *>::iterator iter = loopList.begin(); iter != loopList.end(); iter++)
    {
      std::vector<double> loop = PythonUtility::convertFromPython<std::vector<double>>::get(*iter);
      loops.push_back(loop);
    }


    LOG(DEBUG) << "loops: " << loops;*/
  }  // if own rank is part of this stage of the algorithm
}

template<typename BasisFunctionType>
void ParallelFiberEstimation<BasisFunctionType>::
reorganizeStreamlinePoints(std::vector<std::vector<Vec3>> &streamlineZPoints, std::array<std::array<std::vector<std::vector<Vec3>>,4>,8> &borderPointsSubdomain, std::array<bool,4> &subdomainIsAtBorder)
{
  // borderPointsSubdomain[subdomain index][face_t][z-level][point index]

  // the numbering of the subdomains from 0-7 is as expected (morton numbering)

  // allocate space for borderPointsSubdomain
  for (int subdomainIndex = 0; subdomainIndex < 8; subdomainIndex++)
  {
    for (int faceNo = (int)Mesh::face_t::face0Minus; faceNo <= (int)Mesh::face_t::face1Plus; faceNo++)
    {
      borderPointsSubdomain[subdomainIndex][faceNo].resize(nBorderPointsX_);   // resize to number of z levels
      for (int zLevelIndex = 0; zLevelIndex < nBorderPointsX_; zLevelIndex++)
      {
        borderPointsSubdomain[subdomainIndex][faceNo][zLevelIndex].resize(nBorderPointsX_);   // resize to number of points with same z level per face of subdomain
      }
    }
  }

  // assign sampled points to the data structure borderPointsSubdomain, which contains the points for each subdomain and face, as list of points for each z level

  int nBorderPointsXNew = nBorderPointsX_*2-1;
  int nBorderPointsZNew = nBorderPointsZ_*2 - 1;
  int streamlineIndex = 0;

  // boundary indices for face0Minus and face0Plus (vertical)
  int iBeginVertical = 0;
  int iEndVertical = nBorderPointsXNew;

  if (subdomainIsAtBorder[(int)Mesh::face_t::face1Minus])
    iBeginVertical += 1;

  if (subdomainIsAtBorder[(int)Mesh::face_t::face1Plus])
    iEndVertical -= 1;

  // boundary indices for face1Minus and face1Plus (horizontal)
  int iBeginHorizontal = 0;
  int iEndHorizontal = nBorderPointsXNew;

  if (subdomainIsAtBorder[(int)Mesh::face_t::face0Minus])
    iBeginHorizontal += 1;

  if (subdomainIsAtBorder[(int)Mesh::face_t::face0Plus])
    iEndHorizontal -= 1;

  // face0Minus
  if (!subdomainIsAtBorder[(int)Mesh::face_t::face0Minus])
  {
    // subdomains 0,4
    int pointIndex = iBeginVertical;
    for (int i = iBeginVertical; i < nBorderPointsX_; i++, streamlineIndex++, pointIndex++)
    {
      // loop over bottom half of the streamline points
      for (int zLevelIndex = 0; zLevelIndex <= nBorderPointsZ_; zLevelIndex++)
      {
        const int subdomainIndex = 0;
        borderPointsSubdomain[subdomainIndex][(int)Mesh::face_t::face0Minus][zLevelIndex][pointIndex] = streamlineZPoints[streamlineIndex][zLevelIndex];
      }

      // loop over top half of the streamline points
      for (int zLevelIndex = nBorderPointsZ_; zLevelIndex < nBorderPointsZNew; zLevelIndex++)
      {
        const int subdomainIndex = 4;
        borderPointsSubdomain[subdomainIndex][(int)Mesh::face_t::face0Minus][zLevelIndex][pointIndex] = streamlineZPoints[streamlineIndex][zLevelIndex];
      }
    }

    // subdomains 2,6
    pointIndex = 0;
    streamlineIndex--;  // consider center point twice
    for (int i = nBorderPointsX_-1; i < iEndVertical; i++, streamlineIndex++, pointIndex++)
    {
      // loop over bottom half of the streamline points
      for (int zLevelIndex = 0; zLevelIndex < nBorderPointsZ_; zLevelIndex++)
      {
        const int subdomainIndex = 2;
        borderPointsSubdomain[subdomainIndex][(int)Mesh::face_t::face0Minus][zLevelIndex][pointIndex] = streamlineZPoints[streamlineIndex][zLevelIndex];
      }

      // loop over top half of the streamline points
      for (int zLevelIndex = nBorderPointsZ_; zLevelIndex < nBorderPointsZNew; zLevelIndex++)
      {
        const int subdomainIndex = 6;
        borderPointsSubdomain[subdomainIndex][(int)Mesh::face_t::face0Minus][zLevelIndex][pointIndex] = streamlineZPoints[streamlineIndex][zLevelIndex];
      }
    }
  }

  // face0Plus
  if (!subdomainIsAtBorder[(int)Mesh::face_t::face0Plus])
  {
    // subdomains 1,5
    int pointIndex = iBeginVertical;
    for (int i = iBeginVertical; i < nBorderPointsX_; i++, streamlineIndex++, pointIndex++)
    {
      // loop over bottom half of the streamline points
      for (int zLevelIndex = 0; zLevelIndex <= nBorderPointsZ_; zLevelIndex++)
      {
        const int subdomainIndex = 1;
        borderPointsSubdomain[subdomainIndex][(int)Mesh::face_t::face0Plus][zLevelIndex][pointIndex] = streamlineZPoints[streamlineIndex][zLevelIndex];
      }

      // loop over top half of the streamline points
      for (int zLevelIndex = nBorderPointsZ_; zLevelIndex < nBorderPointsZNew; zLevelIndex++)
      {
        const int subdomainIndex = 5;
        borderPointsSubdomain[subdomainIndex][(int)Mesh::face_t::face0Plus][zLevelIndex][pointIndex] = streamlineZPoints[streamlineIndex][zLevelIndex];
      }
    }

    // subdomains 3,7
    pointIndex = 0;
    streamlineIndex--;  // consider center point twice
    for (int i = nBorderPointsX_-1; i < iEndVertical; i++, streamlineIndex++, pointIndex++)
    {
      // loop over bottom half of the streamline points
      for (int zLevelIndex = 0; zLevelIndex <= nBorderPointsZ_; zLevelIndex++)
      {
        const int subdomainIndex = 3;
        borderPointsSubdomain[subdomainIndex][(int)Mesh::face_t::face0Plus][zLevelIndex][pointIndex] = streamlineZPoints[streamlineIndex][zLevelIndex];
      }

      // loop over top half of the streamline points
      for (int zLevelIndex = nBorderPointsZ_; zLevelIndex < nBorderPointsZNew; zLevelIndex++)
      {
        const int subdomainIndex = 7;
        borderPointsSubdomain[subdomainIndex][(int)Mesh::face_t::face0Plus][zLevelIndex][pointIndex] = streamlineZPoints[streamlineIndex][zLevelIndex];
      }
    }
  }

  // face1Minus (without corner points)
  if (!subdomainIsAtBorder[(int)Mesh::face_t::face1Minus])
  {
    // subdomains 0,4
    int pointIndex = iBeginHorizontal;
    for (int i = iBeginHorizontal; i < nBorderPointsX_; i++, streamlineIndex++, pointIndex++)
    {
      // loop over bottom half of the streamline points
      for (int zLevelIndex = 0; zLevelIndex <= nBorderPointsZ_; zLevelIndex++)
      {
        const int subdomainIndex = 0;
        borderPointsSubdomain[subdomainIndex][(int)Mesh::face_t::face1Minus][zLevelIndex][pointIndex] = streamlineZPoints[streamlineIndex][zLevelIndex];
      }

      // loop over top half of the streamline points
      for (int zLevelIndex = nBorderPointsZ_; zLevelIndex < nBorderPointsZNew; zLevelIndex++)
      {
        const int subdomainIndex = 4;
        borderPointsSubdomain[subdomainIndex][(int)Mesh::face_t::face1Minus][zLevelIndex][pointIndex] = streamlineZPoints[streamlineIndex][zLevelIndex];
      }
    }

    // subdomains 1,5
    pointIndex = 0;
    streamlineIndex--;  // consider center point twice
    for (int i = nBorderPointsX_-1; i < iEndHorizontal; i++, streamlineIndex++, pointIndex++)
    {
      // loop over bottom half of the streamline points
      for (int zLevelIndex = 0; zLevelIndex <= nBorderPointsZ_; zLevelIndex++)
      {
        const int subdomainIndex = 1;
        borderPointsSubdomain[subdomainIndex][(int)Mesh::face_t::face1Minus][zLevelIndex][pointIndex] = streamlineZPoints[streamlineIndex][zLevelIndex];
      }

      // loop over top half of the streamline points
      for (int zLevelIndex = nBorderPointsZ_; zLevelIndex < nBorderPointsZNew; zLevelIndex++)
      {
        const int subdomainIndex = 5;
        borderPointsSubdomain[subdomainIndex][(int)Mesh::face_t::face1Minus][zLevelIndex][pointIndex] = streamlineZPoints[streamlineIndex][zLevelIndex];
      }
    }
  }

  // face1Plus (with corner points)
  if (!subdomainIsAtBorder[(int)Mesh::face_t::face1Plus])
  {
    // subdomains 2,6
    int pointIndex = iBeginHorizontal;
    for (int i = iBeginHorizontal; i < nBorderPointsX_; i++, streamlineIndex++, pointIndex++)
    {
      // loop over bottom half of the streamline points
      for (int zLevelIndex = 0; zLevelIndex <= nBorderPointsZ_; zLevelIndex++)
      {
        const int subdomainIndex = 2;
        borderPointsSubdomain[subdomainIndex][(int)Mesh::face_t::face1Plus][zLevelIndex][pointIndex] = streamlineZPoints[streamlineIndex][zLevelIndex];
      }

      // loop over top half of the streamline points
      for (int zLevelIndex = nBorderPointsZ_; zLevelIndex < nBorderPointsZNew; zLevelIndex++)
      {
        const int subdomainIndex = 6;
        borderPointsSubdomain[subdomainIndex][(int)Mesh::face_t::face1Plus][zLevelIndex][pointIndex] = streamlineZPoints[streamlineIndex][zLevelIndex];
      }
    }

    // subdomains 3,7
    pointIndex = 0;
    streamlineIndex--;  // consider center point twice
    for (int i = nBorderPointsX_-1; i < iEndHorizontal; i++, streamlineIndex++, pointIndex++)
    {
      // loop over bottom half of the streamline points
      for (int zLevelIndex = 0; zLevelIndex <= nBorderPointsZ_; zLevelIndex++)
      {
        const int subdomainIndex = 3;
        borderPointsSubdomain[subdomainIndex][(int)Mesh::face_t::face1Plus][zLevelIndex][pointIndex] = streamlineZPoints[streamlineIndex][zLevelIndex];
      }

      // loop over top half of the streamline points
      for (int zLevelIndex = nBorderPointsZ_; zLevelIndex < nBorderPointsZNew; zLevelIndex++)
      {
        const int subdomainIndex = 7;
        borderPointsSubdomain[subdomainIndex][(int)Mesh::face_t::face1Plus][zLevelIndex][pointIndex] = streamlineZPoints[streamlineIndex][zLevelIndex];
      }
    }
  }

  // horizontal center line (with corner points)
  // subdomains 0,4
  int pointIndex = iBeginHorizontal;
  for (int i = iBeginHorizontal; i < nBorderPointsX_; i++, streamlineIndex++, pointIndex++)
  {
    // loop over bottom half of the streamline points
    for (int zLevelIndex = 0; zLevelIndex <= nBorderPointsZ_; zLevelIndex++)
    {
      const int subdomainIndex = 0;
      LOG(DEBUG) << borderPointsSubdomain[subdomainIndex][(int)Mesh::face_t::face1Plus].size() << " levels ";
      LOG(DEBUG) << borderPointsSubdomain[subdomainIndex][(int)Mesh::face_t::face1Plus][zLevelIndex].size() << " points";
      LOG(DEBUG) << "z: " << zLevelIndex << ", i=" << i;
      borderPointsSubdomain[subdomainIndex][(int)Mesh::face_t::face1Plus][zLevelIndex][pointIndex] = streamlineZPoints[streamlineIndex][zLevelIndex];
    }

    // loop over top half of the streamline points
    for (int zLevelIndex = nBorderPointsZ_; zLevelIndex < nBorderPointsZNew; zLevelIndex++)
    {
      const int subdomainIndex = 4;
      borderPointsSubdomain[subdomainIndex][(int)Mesh::face_t::face1Plus][zLevelIndex][pointIndex] = streamlineZPoints[streamlineIndex][zLevelIndex];
    }
  }

  // subdomains 1,5
  pointIndex = 0;
  streamlineIndex--;  // consider center point twice
  for (int i = nBorderPointsX_-1; i < iEndHorizontal; i++, streamlineIndex++, pointIndex++)
  {
    // loop over bottom half of the streamline points
    for (int zLevelIndex = 0; zLevelIndex <= nBorderPointsZ_; zLevelIndex++)
    {
      const int subdomainIndex = 1;
      borderPointsSubdomain[subdomainIndex][(int)Mesh::face_t::face1Plus][zLevelIndex][pointIndex] = streamlineZPoints[streamlineIndex][zLevelIndex];
    }

    // loop over top half of the streamline points
    for (int zLevelIndex = nBorderPointsZ_; zLevelIndex < nBorderPointsZNew; zLevelIndex++)
    {
      const int subdomainIndex = 5;
      borderPointsSubdomain[subdomainIndex][(int)Mesh::face_t::face1Plus][zLevelIndex][pointIndex] = streamlineZPoints[streamlineIndex][zLevelIndex];
    }
  }

  // subdomains 2,6
  pointIndex = iBeginHorizontal;
  for (int i = iBeginHorizontal; i < nBorderPointsX_; i++, streamlineIndex++, pointIndex++)
  {
    // loop over bottom half of the streamline points
    for (int zLevelIndex = 0; zLevelIndex <= nBorderPointsZ_; zLevelIndex++)
    {
      const int subdomainIndex = 2;
      borderPointsSubdomain[subdomainIndex][(int)Mesh::face_t::face1Minus][zLevelIndex][pointIndex] = streamlineZPoints[streamlineIndex][zLevelIndex];
    }

    // loop over top half of the streamline points
    for (int zLevelIndex = nBorderPointsZ_; zLevelIndex < nBorderPointsZNew; zLevelIndex++)
    {
      const int subdomainIndex = 6;
      borderPointsSubdomain[subdomainIndex][(int)Mesh::face_t::face1Minus][zLevelIndex][pointIndex] = streamlineZPoints[streamlineIndex][zLevelIndex];
    }
  }

  // subdomains 3,7
  pointIndex = 0;
  streamlineIndex--;  // consider center point twice
  for (int i = nBorderPointsX_-1; i < iEndHorizontal; i++, streamlineIndex++, pointIndex++)
  {
    // loop over bottom half of the streamline points
    for (int zLevelIndex = 0; zLevelIndex <= nBorderPointsZ_; zLevelIndex++)
    {
      const int subdomainIndex = 3;
      borderPointsSubdomain[subdomainIndex][(int)Mesh::face_t::face1Minus][zLevelIndex][pointIndex] = streamlineZPoints[streamlineIndex][zLevelIndex];
    }

    // loop over top half of the streamline points
    for (int zLevelIndex = nBorderPointsZ_; zLevelIndex < nBorderPointsZNew; zLevelIndex++)
    {
      const int subdomainIndex = 7;
      borderPointsSubdomain[subdomainIndex][(int)Mesh::face_t::face1Minus][zLevelIndex][pointIndex] = streamlineZPoints[streamlineIndex][zLevelIndex];
    }
  }

  // vertical center line (with corner points and center point)
  // subdomains 0,4
  pointIndex = iBeginVertical;
  for (int i = iBeginVertical; i < nBorderPointsX_; i++, streamlineIndex++, pointIndex++)
  {
    // loop over bottom half of the streamline points
    for (int zLevelIndex = 0; zLevelIndex <= nBorderPointsZ_; zLevelIndex++)
    {
      const int subdomainIndex = 0;
      borderPointsSubdomain[subdomainIndex][(int)Mesh::face_t::face0Plus][zLevelIndex][pointIndex] = streamlineZPoints[streamlineIndex][zLevelIndex];
    }

    // loop over top half of the streamline points
    for (int zLevelIndex = nBorderPointsZ_; zLevelIndex < nBorderPointsZNew; zLevelIndex++)
    {
      const int subdomainIndex = 4;
      borderPointsSubdomain[subdomainIndex][(int)Mesh::face_t::face0Plus][zLevelIndex][pointIndex] = streamlineZPoints[streamlineIndex][zLevelIndex];
    }
  }

  // subdomains 1,5
  pointIndex = 0;
  streamlineIndex--;  // consider center point twice
  for (int i = nBorderPointsX_-1; i < iEndVertical; i++, streamlineIndex++, pointIndex++)
  {
    // loop over bottom half of the streamline points
    for (int zLevelIndex = 0; zLevelIndex <= nBorderPointsZ_; zLevelIndex++)
    {
      const int subdomainIndex = 1;
      borderPointsSubdomain[subdomainIndex][(int)Mesh::face_t::face0Minus][zLevelIndex][pointIndex] = streamlineZPoints[streamlineIndex][zLevelIndex];
    }

    // loop over top half of the streamline points
    for (int zLevelIndex = nBorderPointsZ_; zLevelIndex < nBorderPointsZNew; zLevelIndex++)
    {
      const int subdomainIndex = 5;
      borderPointsSubdomain[subdomainIndex][(int)Mesh::face_t::face0Minus][zLevelIndex][pointIndex] = streamlineZPoints[streamlineIndex][zLevelIndex];
    }
  }

  // subdomains 2,6
  pointIndex = iBeginVertical;
  for (int i = iBeginVertical; i < nBorderPointsX_; i++, streamlineIndex++, pointIndex++)
  {
    // loop over bottom half of the streamline points
    for (int zLevelIndex = 0; zLevelIndex <= nBorderPointsZ_; zLevelIndex++)
    {
      const int subdomainIndex = 2;
      borderPointsSubdomain[subdomainIndex][(int)Mesh::face_t::face0Plus][zLevelIndex][pointIndex] = streamlineZPoints[streamlineIndex][zLevelIndex];
    }

    // loop over top half of the streamline points
    for (int zLevelIndex = nBorderPointsZ_; zLevelIndex < nBorderPointsZNew; zLevelIndex++)
    {
      const int subdomainIndex = 6;
      borderPointsSubdomain[subdomainIndex][(int)Mesh::face_t::face0Plus][zLevelIndex][pointIndex] = streamlineZPoints[streamlineIndex][zLevelIndex];
    }
  }

  // subdomains 3,7
  pointIndex = 0;
  streamlineIndex--;  // consider center point twice
  for (int i = nBorderPointsX_-1; i < iEndVertical; i++, streamlineIndex++, pointIndex++)
  {
    // loop over bottom half of the streamline points
    for (int zLevelIndex = 0; zLevelIndex <= nBorderPointsZ_; zLevelIndex++)
    {
      const int subdomainIndex = 3;
      borderPointsSubdomain[subdomainIndex][(int)Mesh::face_t::face0Minus][zLevelIndex][pointIndex] = streamlineZPoints[streamlineIndex][zLevelIndex];
    }

    // loop over top half of the streamline points
    for (int zLevelIndex = nBorderPointsZ_; zLevelIndex < nBorderPointsZNew; zLevelIndex++)
    {
      const int subdomainIndex = 7;
      borderPointsSubdomain[subdomainIndex][(int)Mesh::face_t::face0Minus][zLevelIndex][pointIndex] = streamlineZPoints[streamlineIndex][zLevelIndex];
    }
  }

  // set index to beginning of next face
  assert(streamlineIndex == streamlineZPoints.size());

}

};  // namespace
