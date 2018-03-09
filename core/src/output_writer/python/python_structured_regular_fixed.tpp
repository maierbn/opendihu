#include "output_writer/python/python.h"

#include <Python.h>   // this has to be the first included header
#include <iostream>

#include "easylogging++.h"
#include "utility/python_utility.h"
#include "output_writer/python/python_base.h"

namespace OutputWriter
{

template<int D, typename BasisFunctionType>
PyObject *Python<BasisOnMesh::BasisOnMesh<Mesh::StructuredRegularFixedOfDimension<D>,BasisFunctionType>>::
buildPyDataObject(std::vector<std::shared_ptr<FieldVariable::FieldVariable<BasisOnMeshType>>> fieldVariables, int timeStepNo, double currentTime)
{
  // build python dict containing all information
  // data = {
  //   "meshType" : "RegularFixed",
  //   "dimension": dim,
  //   "nElements" : [x,y,z],
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
  
  PyObject *pyData = PythonBase<BasisOnMesh::BasisOnMesh<Mesh::StructuredRegularFixedOfDimension<D>,BasisFunctionType>>::buildPyFieldVariablesObject(fieldVariables);
  
  // prepare number of elements in the dimensions
  std::array<element_no_t, BasisOnMeshType::Mesh::dim()> nElementsPerCoordinateDirection = fieldVariables.front()->nElementsPerCoordinateDirection();
  std::array<long, BasisOnMeshType::Mesh::dim()> nElementsPerCoordinateDirectionArray;
  
  std::copy(nElementsPerCoordinateDirection.begin(), nElementsPerCoordinateDirection.end(), nElementsPerCoordinateDirectionArray.begin());
  PyObject *pyNElements = PythonUtility::convertToPythonList<BasisOnMeshType::Mesh::dim()>(nElementsPerCoordinateDirectionArray);
  
  std::string basisFunction = BasisOnMeshType::BasisFunction::getBasisFunctionString();
  int basisOrder = BasisOnMeshType::BasisFunction::getBasisOrder();
  
  LOG(DEBUG) << "PythonRegularFixed";
  
  // build python dict that will contain all information and data
  PyObject *data = Py_BuildValue("{s s, s i, s O, s s, s i, s O, s i, s d}", "meshType", "StructuredRegularFixed",
                                 "dimension", D, "nElements", pyNElements, "basisFunction", basisFunction.c_str(), "basisOrder", basisOrder, "data", pyData, 
                                 "timeStepNo", timeStepNo, "currentTime", currentTime);
  
  return data;
} 
 
};