#include "output_writer/python/python.h"

#include <Python.h>   // this has to be the first included header
#include <iostream>

#include "easylogging++.h"
#include "utility/python_utility.h"

namespace OutputWriter
{

template<int D, typename BasisFunctionType, typename OutputFieldVariablesType>
PyObject *Python<BasisOnMesh::BasisOnMesh<Mesh::StructuredDeformableOfDimension<D>,BasisFunctionType>,OutputFieldVariablesType>::
buildPyDataObject(OutputFieldVariablesType fieldVariables,
                  int timeStepNo, double currentTime, bool onlyNodalValues)
{
  // build python dict containing all information
  // data = {
  //   "meshType" : "StructuredDeformable",
  //   "dimension": dim,
  //   "nElements" : [x, y, z],
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
  
  PyObject *pyData = PythonBase<OutputFieldVariablesType>::buildPyFieldVariablesObject(fieldVariables, onlyNodalValues);
  
  // prepare number of elements in the dimensions
  std::array<element_no_t, BasisOnMeshType::dim()> nElementsPerCoordinateDirection = std::get<0>(fieldVariables)->nElementsPerCoordinateDirection();
  std::array<long, BasisOnMeshType::dim()> nElementsPerCoordinateDirectionArray;
  
  std::copy(nElementsPerCoordinateDirection.begin(), nElementsPerCoordinateDirection.end(), nElementsPerCoordinateDirectionArray.begin());
  PyObject *pyNElements = PythonUtility::convertToPythonList<BasisOnMeshType::dim()>(nElementsPerCoordinateDirectionArray);
  
  std::string basisFunction = BasisOnMeshType::BasisFunction::getBasisFunctionString();
  int basisOrder = BasisOnMeshType::BasisFunction::getBasisOrder();
  
  LOG(DEBUG) << "PythonStructuredDeformable";
  
  // build python dict that will contain all information and data
  PyObject *data = Py_BuildValue("{s s, s i, s O, s s, s i, s O, s O, s i, s d}", "meshType", "StructuredDeformable",
                                 "dimension", D, "nElements", pyNElements, 
                                 "basisFunction", basisFunction.c_str(), "basisOrder", basisOrder, 
                                 "onlyNodalValues", onlyNodalValues ? Py_True: Py_False,
                                 "data", pyData, 
                                 "timeStepNo", timeStepNo, "currentTime", currentTime);
  
  return data;
} 
 
};