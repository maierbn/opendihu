#include "postprocessing/parallel_fiber_estimation/parallel_fiber_estimation.h"

namespace Postprocessing
{

template<typename BasisFunctionType>
void ParallelFiberEstimation<BasisFunctionType>::
createMesh(std::array<std::vector<std::vector<Vec3>>,4> &borderPoints, std::vector<Vec3> &nodePositions, std::array<int,3> &nElementsPerCoordinateDirectionLocal)
{
  int subdomainNNodesX;
  int subdomainNNodesY;
  int subdomainNNodesZ;

  LOG(DEBUG) << "createMesh";

#ifdef USE_CHECKPOINT_MESH  // load from file
  std::stringstream filename;
  filename << "checkpoints/checkpoint_mesh_l" << level_ << "_" << currentRankSubset_->ownRankNo() << ".csv";
  std::ifstream file(filename.str().c_str());

  if (!file.is_open())
  {
    LOG(FATAL) << "Could not open file \"" << filename.str() << "\".";
  }

  file >> subdomainNNodesX >> subdomainNNodesY >> subdomainNNodesZ;
  nElementsPerCoordinateDirectionLocal[0] = subdomainNNodesX-1;
  nElementsPerCoordinateDirectionLocal[1] = subdomainNNodesY-1;
  nElementsPerCoordinateDirectionLocal[2] = subdomainNNodesZ-1;

  while(!file.eof())
  {
    Vec3 nodePosition;
    int i = 0;
    for (i = 0; i < 3; i++)
    {
      file >> nodePosition[i];
      if (file.eof())
        break;
    }
    if (i == 3)
    {
      nodePositions.push_back(nodePosition);
    }
  }

  file.close();

#else
  // call stl_create_mesh.create_3d_mesh_from_border_points_faces
  PyObject *borderPointsFacesPy = PythonUtility::convertToPython<std::array<std::vector<std::vector<Vec3>>,4>>::get(borderPoints);
  PythonUtility::checkForError();

  //LOG(DEBUG) << PythonUtility::getString(borderPointsFacesPy);
  LOG(DEBUG) << "call function create_3d_mesh_from_border_points_faces";

  PyObject *meshData = PyObject_CallFunction(functionCreate3dMeshFromBorderPointsFaces_, "(O,O,i)", borderPointsFacesPy, (improveMesh_? Py_True : Py_False), level_);
  PythonUtility::checkForError();

  if (meshData == Py_None)
  {
    LOG(FATAL) << "Python function create_3d_mesh_from_border_points_faces returned None!";
  }

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

  PyObject *object = PythonUtility::getOptionPyObject(meshData, "node_positions", "");
  nodePositions = PythonUtility::convertFromPython<std::vector<Vec3>>::get(object);
  nElementsPerCoordinateDirectionLocal = PythonUtility::getOptionArray<int,3>(meshData, "n_linear_elements_per_coordinate_direction", "", std::array<int,3>({0,0,0}));

  subdomainNNodesX = nElementsPerCoordinateDirectionLocal[0]+1;
  subdomainNNodesY = nElementsPerCoordinateDirectionLocal[1]+1;
  subdomainNNodesZ = nElementsPerCoordinateDirectionLocal[2]+1;

  LOG(DEBUG) << "subdomainNNodes: " << subdomainNNodesX << " x " << subdomainNNodesY << " x " << subdomainNNodesZ;

  if (subdomainNNodesX != nBorderPointsXNew_ || subdomainNNodesY != nBorderPointsXNew_ || subdomainNNodesZ != nBorderPointsZ_)
  {
    PyObject_CallFunction(functionOutputBorderPoints_, "s i i O f", "xx_failed_border_points", currentRankSubset_->ownRankNo(), level_,
                          PythonUtility::convertToPython<std::array<std::vector<std::vector<Vec3>>,4>>::get(borderPoints), 0.03);
    PythonUtility::checkForError();
  }
  assert(subdomainNNodesX == nBorderPointsXNew_);
  assert(subdomainNNodesY == nBorderPointsXNew_);
  assert(subdomainNNodesZ == nBorderPointsZ_);
/*
  // revert order of node positions
  std::vector<Vec3> nodePositions(nodePositionsOrderReversed.size());
  for (int z = 0; z < subdomainNNodesZ; z++)
  {
    for (int y = 0; y < subdomainNNodesY; y++)
    {
      for (int x = 0; x < subdomainNNodesX; x++)
      {
        nodePositions[z * subdomainNNodesX*subdomainNNodesY + y * subdomainNNodesY + x]
          = nodePositionsOrderReversed[z * subdomainNNodesX*subdomainNNodesY + (subdomainNNodesY-1-y) * subdomainNNodesY + (subdomainNNodesX-1-x)];
      }
    }
  }*/

  //std::vector<Vec3> &nodePositions = nodePositionsOrderReversed;

#ifdef WRITE_CHECKPOINT_MESH
  std::stringstream filename;
  filename << "checkpoints/checkpoint_mesh_l" << level_ << "_" << currentRankSubset_->ownRankNo() << ".csv";
  std::ofstream file(filename.str().c_str(), std::ios::out | std::ios::trunc);
  assert(file.is_open());

  file << subdomainNNodesX << " " << subdomainNNodesY << " " << subdomainNNodesZ << " " << std::endl;
  for (std::vector<Vec3>::iterator iter = nodePositions.begin(); iter != nodePositions.end(); iter++)
  {
    for (int i = 0; i < 3; i++)
    {
      file << (*iter)[i] << " ";
    }
  }

  file.close();

#endif

#endif
  //LOG(DEBUG) << "nodePositions: " << nodePositions;
  LOG(DEBUG) << "nElementsPerCoordinateDirectionLocal: " << nElementsPerCoordinateDirectionLocal;

#ifndef NDEBUG
#ifdef STL_OUTPUT
  PyObject_CallFunction(functionOutputPoints_, "s i i O f", "03_mesh_points", currentRankSubset_->ownRankNo(), level_,
                        PythonUtility::convertToPython<std::vector<Vec3>>::get(nodePositions), 0.05);
  PythonUtility::checkForError();
#endif
#endif
}

} // namespace
