#include "output_writer/python/python.h"

#include <Python.h>   // this has to be the first included header
#include <iostream>

#include "easylogging++.h"
#include "utility/python_utility.h"

namespace OutputWriter
{

template<int D, typename BasisFunctionType>
PyObject *Python<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformable<D>,BasisFunctionType>>::
buildPyDataObject(std::vector<std::shared_ptr<FieldVariable::FieldVariable<BasisOnMeshType>>> fieldVariables, int timeStepNo, double currentTime)
{
  // build python dict containing all information
  // data = {
  //   "meshType" : "UnstructuredDeformable",
  //   "dimension": dim,
  //   "nElements" : nElements,
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
  
 
  LOG(DEBUG) << "build pyData, " << fieldVariables.size();
 
  // build python object for data  
  PyObject *pyData = PythonBase<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformable<D>,BasisFunctionType>>::buildPyFieldVariablesObject(fieldVariables);

  LOG(DEBUG) << "get mesh";
  
  std::shared_ptr<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformable<D>,BasisFunctionType>> mesh = fieldVariables.front()->mesh();
  
  LOG(DEBUG) << pyData << " build data";
  LOG(DEBUG) << "PythonUnStructuredDeformable";
  
  // build python dict that will contain all information and data
  PyObject *data = Py_BuildValue("{s s, s i, s i, s O, s i, s d}", "meshType", "UnstructuredDeformable",
                                 "dimension", D, "nElements", mesh->nElements(), "data", pyData, 
                                 "timeStepNo", timeStepNo, "currentTime", currentTime);
  
  LOG(DEBUG) << data << " done";
  
  return data;
} 
 
};