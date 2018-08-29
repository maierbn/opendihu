#include "output_writer/python/python.h"

#include <Python.h>   // this has to be the first included header
#include <iostream>

#include "easylogging++.h"
#include "utility/python_utility.h"
#include "output_writer/python/python_base.h"

namespace OutputWriter
{

template<int D, typename BasisFunctionType, typename OutputFieldVariablesType>
PyObject *Python<BasisOnMesh::BasisOnMesh<Mesh::StructuredRegularFixedOfDimension<D>,BasisFunctionType>,OutputFieldVariablesType>::
buildPyDataObject(OutputFieldVariablesType fieldVariables,
                  std::string meshName, int timeStepNo, double currentTime, bool onlyNodalValues)
{
  // build python dict containing all information
  // data = {
  //   "meshType" : "RegularFixed",
  //   "dimension": dim,
  //   "nElements" : [x,y,z],
  //   "basisFunction" : "Lagrange",
  //   "basisOrder" : "1",
  //   "onlyNodalValues" : True,
  //   "data" : [
  //      {"name" : "fieldVariableName",
  //       "components" : [
  //           {"name" : "componentName", "values": data},
  //       ]
  //      },
  //   ]
  //   "timeStepNo" : timeStepNo,
  //   "currentTime" : currentTime
  // }

  // build python object for data
  std::shared_ptr<Mesh::Mesh> meshBase;
  PyObject *pyData = PythonBase<OutputFieldVariablesType>::buildPyFieldVariablesObject(fieldVariables, meshName, onlyNodalValues, meshBase);

  // cast mesh to its real type
  typedef BasisOnMesh::BasisOnMesh<Mesh::StructuredRegularFixedOfDimension<D>,BasisFunctionType> BasisOnMeshType;
  std::shared_ptr<BasisOnMeshType> mesh = std::static_pointer_cast<BasisOnMeshType>(meshBase);

  // convert nElementsPerCoordinateDirectionGlobal to python list
  std::array<global_no_t, BasisOnMeshType::dim()> nElementsPerCoordinateDirectionGlobal = mesh->nElementsPerCoordinateDirectionGlobal();
  std::array<long, BasisOnMeshType::dim()> nElementsPerCoordinateDirectionGlobalArray;

  std::copy(nElementsPerCoordinateDirectionGlobal.begin(), nElementsPerCoordinateDirectionGlobal.end(), nElementsPerCoordinateDirectionGlobalArray.begin());
  PyObject *pyNElementsGlobal = PythonUtility::convertToPythonList<BasisOnMeshType::dim()>(nElementsPerCoordinateDirectionGlobalArray);

  // convert nElementsPerCoordinateDirectionLocal to python list
  std::array<element_no_t, BasisOnMeshType::dim()> nElementsPerCoordinateDirectionLocal = mesh->nElementsPerCoordinateDirectionLocal();
  std::array<long, BasisOnMeshType::dim()> nElementsPerCoordinateDirectionLocalArray;

  std::copy(nElementsPerCoordinateDirectionLocal.begin(), nElementsPerCoordinateDirectionLocal.end(), nElementsPerCoordinateDirectionLocalArray.begin());
  PyObject *pyNElementsLocal = PythonUtility::convertToPythonList<BasisOnMeshType::dim()>(nElementsPerCoordinateDirectionLocalArray);

  std::array<long, BasisOnMeshType::dim()> beginNodeGlobal;
  std::array<long, BasisOnMeshType::dim()> hasFullNumberOfNodes;

  for (int i = 0; i < BasisOnMeshType::dim(); i++)
  {
    beginNodeGlobal[i] = mesh->meshPartition()->beginNodeGlobal(i);
    hasFullNumberOfNodes[i] = mesh->meshPartition()->hasFullNumberOfNodes(i);
  }

  PyObject *pyBeginNodeGlobal = PythonUtility::convertToPythonList<BasisOnMeshType::dim()>(beginNodeGlobal);
  PyObject *pyHasFullNumberOfNodes = PythonUtility::convertToPythonList<BasisOnMeshType::dim()>(hasFullNumberOfNodes);

  // convert basis function information
  std::string basisFunction = BasisOnMeshType::BasisFunction::getBasisFunctionString();
  int basisOrder = BasisOnMeshType::BasisFunction::getBasisOrder();

  int nRanks = mesh->meshPartition()->nRanks();
  int ownRankNo = mesh->meshPartition()->ownRankNo();

  LOG(DEBUG) << "PythonRegularFixed";

  // start critical section for python API calls
  PythonUtility::GlobalInterpreterLock lock;
  
  // build python dict that will contain all information and data
  PyObject *data = Py_BuildValue("{s s, s i, s O, s O, s O, s O, s s, s i, s O, s i, s i, s O, s i, s d}",
                                 "meshType", "StructuredRegularFixed",
                                 "dimension", D, "nElementsGlobal", pyNElementsGlobal, "nElementsLocal", pyNElementsLocal,
                                 "beginNodeGlobal", pyBeginNodeGlobal, "hasFullNumberOfNodes", pyHasFullNumberOfNodes,
                                 "basisFunction", basisFunction.c_str(), "basisOrder", basisOrder,
                                 "onlyNodalValues", onlyNodalValues ? Py_True: Py_False,
                                 "nRanks", nRanks, "ownRankNo", ownRankNo,
                                 "data", pyData,
                                 "timeStepNo", timeStepNo, "currentTime", currentTime);

  return data;
}

};
