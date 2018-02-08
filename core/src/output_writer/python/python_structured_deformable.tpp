#include "output_writer/python/python.h"

#include <Python.h>   // this has to be the first included header
#include <iostream>

#include "easylogging++.h"
#include "utility/python_utility.h"

namespace OutputWriter
{

template<int D, typename BasisFunctionType>
PyObject *Python<BasisOnMesh::BasisOnMesh<Mesh::StructuredDeformable<D>,BasisFunctionType>>::
buildPyDataObject(std::vector<std::shared_ptr<FieldVariable::FieldVariable<BasisOnMeshType>>> fieldVariables, int timeStepNo, double currentTime)
{
  // build python dict containing all information
  // data = {
  //   "meshType" : "StructuredDeformable",
  //   "dimension": dim,
  //   "nElements" : [x, y, z],
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
  
  PyObject *pyData = PythonBase<BasisOnMesh::BasisOnMesh<Mesh::StructuredDeformable<D>,BasisFunctionType>>::buildPyFieldVariablesObject(fieldVariables);
  
  // prepare number of elements in the dimensions
  std::array<element_no_t, BasisOnMeshType::dim()> nElementsPerDimension = fieldVariables.front()->nElementsPerDimension();
  std::array<long, BasisOnMeshType::dim()> nElementsPerDimensionArray;
  
  std::copy(nElementsPerDimension.begin(), nElementsPerDimension.end(), nElementsPerDimensionArray.begin());
  PyObject *pyNElements = PythonUtility::convertToPythonList<BasisOnMeshType::dim()>(nElementsPerDimensionArray);
  
  LOG(DEBUG) << "PythonStructuredDeformable";
  
  // build python dict that will contain all information and data
  PyObject *data = Py_BuildValue("{s s, s i, s O, s O, s i, s d}", "meshType", "StructuredDeformable",
                                 "dimension", D, "nElements", pyNElements, "data", pyData, 
                                 "timeStepNo", timeStepNo, "currentTime", currentTime);
  
  return data;
} 
 
};