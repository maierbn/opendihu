#include "output_writer/python/python_base.h"

#include <Python.h>   // this has to be the first included header
#include <iostream>
#include <tuple>

#include "easylogging++.h"
#include "utility/python_utility.h"
#include "output_writer/python/loop_build_py_field_variable_object.h"

namespace OutputWriter
{

template<typename OutputFieldVariablesType>
PyObject *PythonBase<OutputFieldVariablesType>::
buildPyFieldVariablesObject(OutputFieldVariablesType fieldVariables, std::string meshName, bool onlyNodalValues, std::shared_ptr<Mesh::Mesh> &mesh)
{
  // build python dict containing field variables
  // [
  //   {"name" : "fieldVariableName",
  //    "components" : [
  //      {"name" : "componentName", "values": data},
  //    ]
  //   },
  // ]

  const int nFieldVariables = std::tuple_size<OutputFieldVariablesType>::value;
  VLOG(2) << "buildPyFieldVariablesObject for " << nFieldVariables << " field variables";

  PyObject *pyData = PyList_New((Py_ssize_t)nFieldVariables);

  int fieldVariableIndex = 0;
  PythonLoopOverTuple::loopBuildPyFieldVariableObject<OutputFieldVariablesType>(fieldVariables, fieldVariableIndex, meshName, pyData, onlyNodalValues, mesh);

  return pyData;
}

};