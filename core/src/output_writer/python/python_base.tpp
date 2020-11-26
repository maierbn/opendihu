#include "output_writer/python/python_base.h"

#include <Python.h>   // this has to be the first included header
#include <iostream>
#include <tuple>

#include "easylogging++.h"
#include "utility/python_utility.h"
#include "output_writer/python/loop_build_py_field_variable_object.h"
#include "output_writer/loop_count_n_field_variables_of_mesh.h"

namespace OutputWriter
{

template<typename FieldVariablesForOutputWriterType>
PyObject *PythonBase<FieldVariablesForOutputWriterType>::
buildPyFieldVariablesObject(FieldVariablesForOutputWriterType fieldVariables, std::string meshName, bool onlyNodalValues, std::shared_ptr<Mesh::Mesh> &mesh)
{
  // build python dict containing field variables
  // [
  //   {"name" : "fieldVariableName",
  //    "components" : [
  //      {"name" : "componentName", "values": data},
  //    ]
  //   },
  // ]


  int nFieldVariablesInMesh = 0;
  LoopOverTuple::loopCountNFieldVariablesOfMesh(fieldVariables, meshName, nFieldVariablesInMesh);
  VLOG(2) << "buildPyFieldVariablesObject for " << nFieldVariablesInMesh << " field variables: "
    << StringUtility::demangle(typeid(FieldVariablesForOutputWriterType).name());

  PyObject *pyData = PyList_New((Py_ssize_t)nFieldVariablesInMesh);

  int fieldVariableIndex = 0;
  PythonLoopOverTuple::loopBuildPyFieldVariableObject<FieldVariablesForOutputWriterType>(fieldVariables, fieldVariableIndex, meshName, pyData, onlyNodalValues, mesh);

  return pyData;
}

}  // namespace
