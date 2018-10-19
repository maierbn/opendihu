#include "postprocessing/parallel_fiber_estimation.h"

#include <algorithm>
#include <petscvec.h>

#include "utility/python_utility.h"
#include "mesh/face_t.h"
#include "partition/01_mesh_partition.h"

namespace Postprocessing
{

template<typename BasisFunctionType>
ParallelFiberEstimation<BasisFunctionType>::
ParallelFiberEstimation(DihuContext context) :
  context_(context["ParallelFiberEstimation"]), problem_(nullptr), data_(context_)
{
  LOG(TRACE) << "ParallelFiberEstimation::ParallelFiberEstimation()";

  specificSettings_ = context_.getPythonConfig();
  outputWriterManager_.initialize(context_, specificSettings_);

  stlFilename_ = PythonUtility::getOptionString(specificSettings_, "stlFilename", "");
  bottomZClip_ = PythonUtility::getOptionInt(specificSettings_, "bottomZClip", 0);
  topZClip_ = PythonUtility::getOptionInt(specificSettings_, "topZClip", 100);
  nElementsZPerSubdomain_ = PythonUtility::getOptionInt(specificSettings_, "nElementsZPerSubdomain", 13);
}

template<typename BasisFunctionType>
void ParallelFiberEstimation<BasisFunctionType>::
initialize()
{
  LOG(TRACE) << "ParallelFiberEstimation::initialize";

  data_.initialize();
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
  LOG(DEBUG) << "generateParallelMesh";

  int nBorderPointsX = 4;
  int nBorderPoints = 4*nBorderPointsX;
  int level = 0;

  // get loops of whole domain

  // run python script to generate loops for the whole volume
  PyObject* moduleStlCreateRings = PyImport_ImportModule("stl_create_rings");

  if (moduleStlCreateRings == NULL)
  {
    LOG(FATAL) << "Could not load python module stl_create_rings. Ensure that the PYTHONPATH environment variable contains the path to opendihu/scripts/geometry_manipulation" << std::endl
      << "Execute the following: " << std::endl << "export PYTHONPATH=$PYTHONPATH:" << OPENDIHU_HOME << "/scripts/geometry_manipulation";
  }

  PyObject* functionCreateRings = PyObject_GetAttrString(moduleStlCreateRings, "create_rings");
  assert(functionCreateRings);

  PyObject* loopsPy = PyObject_CallFunction(functionCreateRings, "s i i i O", stlFilename_.c_str(), bottomZClip_, topZClip_, nElementsZPerSubdomain_+1, Py_False);
  assert(loopsPy);

  PyObject* moduleStlCreateMesh = PyImport_ImportModule("stl_create_mesh");
  assert(moduleStlCreateMesh);

  PyObject* functionRingsToBorderPoints = PyObject_GetAttrString(moduleStlCreateMesh, "rings_to_border_points");
  assert(functionRingsToBorderPoints);

  PyObject* returnValue = PyObject_CallFunction(functionRingsToBorderPoints, "O i", loopsPy, nBorderPoints);
  if (returnValue == NULL)
  {
    LOG(FATAL) << "rings_to_border_points did not return a valid value";
  }

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

  PyObject* functionBorderPointLoopsToList = PyObject_GetAttrString(moduleStlCreateMesh, "border_point_loops_to_list");
  assert(functionRingsToBorderPoints);

  PyObject* borderPointLoops = PyObject_CallFunction(functionBorderPointLoopsToList, "O", borderPointsPy);


  std::vector<std::vector<Vec3>> loops = PythonUtility::convertFromPython<std::vector<std::vector<Vec3>>>::get(borderPointLoops);
  std::vector<double> lengths = PythonUtility::convertFromPython<std::vector<double>>::get(lengthsPy);
  //LOG(DEBUG) << "loops: " << loops << ", lengths: " << lengths;

  // rearrange the border points from the loops to the portitions of the faces
  std::array<std::vector<std::vector<Vec3>>,4> borderPoints;  // borderPoints[face_t][z-level][x]

  // loop over faces: face0Minus = 0, face0Plus, face1Minus, face1Plus
  for (int face = Mesh::face_t::face0Minus; face <= Mesh::face_t::face1Plus; face++)
  {
    borderPoints[face].resize(nElementsZPerSubdomain_+1);
    for (int zIndex = 0; zIndex < nElementsZPerSubdomain_+1; zIndex++)
    {
      borderPoints[face][zIndex].resize(nBorderPointsX);

      if (face == Mesh::face_t::face1Minus)
      {
        std::copy(loops[zIndex].begin(), loops[zIndex].begin() + nBorderPointsX, borderPoints[face][zIndex].begin());
      }
      else if (face == Mesh::face_t::face0Plus)
      {
        std::copy(loops[zIndex].begin() + nBorderPointsX, loops[zIndex].begin() + 2*nBorderPointsX, borderPoints[face][zIndex].begin());
      }
      else if (face == Mesh::face_t::face1Plus)
      {
        std::reverse_copy(loops[zIndex].begin() + 2*nBorderPointsX, loops[zIndex].begin() + 3*nBorderPointsX, borderPoints[face][zIndex].begin());
      }
      else if (face == Mesh::face_t::face0Minus)
      {
        std::reverse_copy(loops[zIndex].begin() + 3*nBorderPointsX, loops[zIndex].end(), borderPoints[face][zIndex].begin());
      }
    }
  }

  // create mesh in own domain, using python, harmonic maps

  PyObject *functionCreate3dMeshFromBorderPointsFaces = PyObject_GetAttrString(moduleStlCreateMesh, "create_3d_mesh_from_border_points_faces");
  assert(functionCreate3dMeshFromBorderPointsFaces);
  PythonUtility::checkForError();

  PyObject *borderPointsFacesPy = PythonUtility::convertToPython<std::array<std::vector<std::vector<Vec3>>,4>>::get(borderPoints);
  PythonUtility::checkForError();

  LOG(DEBUG) << PythonUtility::getString(borderPointsFacesPy);
  LOG(DEBUG) << "call function create_3d_mesh_from_border_points_faces";

  PyObject *meshData = PyObject_CallFunction(functionCreate3dMeshFromBorderPointsFaces, "(O)", borderPointsFacesPy);
  PythonUtility::checkForError();

  LOG(DEBUG) << PythonUtility::getString(meshData);

  std::vector<Vec3> nodePositions;
  PyObject *object = PythonUtility::getOptionPyObject(meshData, "node_positions");
  nodePositions = PythonUtility::convertFromPython<std::vector<Vec3>>::get(object);
  std::array<int,3> nElementsPerCoordinateDirectionLocal = PythonUtility::getOptionArray<int,3>(meshData, "n_linear_elements_per_coordinate_direction", std::array<int,3>({0,0,0}));

  LOG(DEBUG) << "nodePositions: " << nodePositions;
  LOG(DEBUG) << "nElementsPerCoordinateDirectionLocal: " << nElementsPerCoordinateDirectionLocal;

  std::stringstream meshName;
  meshName << "meshLevel" << level;


  std::shared_ptr<Partition::RankSubset> currentRankSubset = std::make_shared<Partition::RankSubset>(0);
  std::array<int,3> nRanksPerCoordinateDirection({1,1,1});

  std::array<global_no_t,3> nElementsPerCoordinateDirectionGlobal;

  // create meshPartition for function space
  std::shared_ptr<Partition::MeshPartition<FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<3>, BasisFunctionType>>> meshPartition
    = context_.partitionManager()->template createPartitioningStructuredLocal<FunctionSpaceType>(
      nElementsPerCoordinateDirectionGlobal, nElementsPerCoordinateDirectionLocal, nRanksPerCoordinateDirection);

  LOG(DEBUG) << "nElementsPerCoordinateDirectionGlobal: " << nElementsPerCoordinateDirectionGlobal;

  // create function space
  context_.partitionManager()->setRankSubsetForNextCreatedMesh(currentRankSubset);
  functionSpace_ = context_.meshManager()->template createFunctionSpaceWithGivenMeshPartition<FunctionSpaceType>(
    meshName.str(), meshPartition, nodePositions, nElementsPerCoordinateDirectionLocal);
  /*
   std::array<global_no_t,FunctionSpace::dim()> &nElementsGlobal,
    const std::array<element_no_t,FunctionSpace::dim()> nElementsLocal,
    const std::array<int,FunctionSpace::dim()> nRanks
   * */
  // solve laplace problem globally
  problem_ = std::make_shared<FiniteElementMethodType>(context_, functionSpace_);

  std::shared_ptr<SpatialDiscretization::DirichletBoundaryConditions<FunctionSpaceType,1>> dirichletBoundaryConditions
    = std::make_shared<SpatialDiscretization::DirichletBoundaryConditions<FunctionSpaceType,1>>();

  typedef typename SpatialDiscretization::DirichletBoundaryConditions<FunctionSpaceType,1>::ElementWithNodes ElementWithNodes;

  std::vector<ElementWithNodes> boundaryConditionElements;
  std::vector<dof_no_t> boundaryConditionNonGhostDofLocalNos;
  std::vector<std::array<double,1>> boundaryConditionValues;

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

          dof_no_t dofLocalNo = functionSpace_->getDofNo(elementNoLocal, elementalDofIndex);
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
  std::fill(boundaryConditionValues.begin(), boundaryConditionValues.end(),  std::array<double,1>({0.0}));
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

          dof_no_t dofLocalNo = functionSpace_->getDofNo(elementNoLocal, elementalDofIndex);
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

  dirichletBoundaryConditions->initialize(functionSpace_, boundaryConditionElements, boundaryConditionNonGhostDofLocalNos, boundaryConditionValues);

  // set boundary conditions to the problem
  problem_->setDirichletBoundaryConditions(dirichletBoundaryConditions);
  problem_->initialize();

  // solve the laplace problem, globally
  problem_->run();



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
}

};
