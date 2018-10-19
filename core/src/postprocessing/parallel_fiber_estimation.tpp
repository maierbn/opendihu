#include "postprocessing/parallel_fiber_estimation.h"

#include <algorithm>
#include <petscvec.h>

#include "utility/python_utility.h"
#include "mesh/face_t.h"

namespace Postprocessing
{

template<typename DiscretizableInTimeType>
ParallelFiberEstimation<DiscretizableInTimeType>::
ParallelFiberEstimation(DihuContext context) :
  context_(context["ParallelFiberEstimation"]), problem_(context_), data_(context_)
{
  LOG(TRACE) << "ParallelFiberEstimation::ParallelFiberEstimation()";

  specificSettings_ = context_.getPythonConfig();
  outputWriterManager_.initialize(context_, specificSettings_);

  stlFilename_ = PythonUtility::getOptionString(specificSettings_, "stlFilename", "");
  bottomZClip_ = PythonUtility::getOptionInt(specificSettings_, "bottomZClip", 0);
  topZClip_ = PythonUtility::getOptionInt(specificSettings_, "topZClip", 100);
  nElementsZPerSubdomain_ = PythonUtility::getOptionInt(specificSettings_, "nElementsZPerSubdomain", 13);
}

template<typename DiscretizableInTimeType>
void ParallelFiberEstimation<DiscretizableInTimeType>::
initialize()
{
  LOG(TRACE) << "ParallelFiberEstimation::initialize";

  // initialize the problem
  problem_.initialize();

  data_.initialize();
}

template<typename DiscretizableInTimeType>
void ParallelFiberEstimation<DiscretizableInTimeType>::
run()
{
  initialize();

  generateParallelMesh();

  // call the method of the underlying problem
  problem_.run();

  // output
  outputWriterManager_.writeOutput(data_);
}

template<typename DiscretizableInTimeType>
void ParallelFiberEstimation<DiscretizableInTimeType>::
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
  std::array<int,3> nElements = PythonUtility::getOptionArray<int,3>(meshData, "n_linear_elements_per_coordinate_direction", std::array<int,3>({0,0,0}));

  LOG(DEBUG) << "nodePositions: " << nodePositions;
  LOG(DEBUG) << "nElements: " << nElements;

  std::stringstream meshName;
  meshName << "meshLevel" << level;

  std::shared_ptr<Partition::RankSubset> currentRankSubset = std::make_shared<Partition::RankSubset>(0);
  context_.partitionManager()->setRankSubsetForNextCreatedMesh(currentRankSubset);
  functionSpace_ = context_.meshManager()->template createFunctionSpace<FunctionSpace>(meshName.str(), nodePositions, nElements);

  // solve laplace problem globally


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

  LOG(FATAL) << "done";
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
