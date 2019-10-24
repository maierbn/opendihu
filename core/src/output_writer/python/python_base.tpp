#include "output_writer/python/python_base.h"

#include <Python.h>   // this has to be the first included header
#include <iostream>
#include <tuple>

#include "easylogging++.h"
#include "utility/python_utility.h"
#include "output_writer/python/loop_build_py_field_variable_object.h"

namespace OutputWriter
{

/** Helper struct that determines the number of field variables.
 */
template<typename T>
struct NFieldVariables
{
  static int get(T &fieldVariables)
  {
    VLOG(2) << "nFieldVariables of " << StringUtility::demangle(typeid(T).name()) << ": " << std::tuple_size<T>::value;
    return std::tuple_size<T>::value;
  }
};

template<typename T>
struct NFieldVariables<std::tuple<std::vector<T>>>
{
  static int get(std::tuple<std::vector<T>> &fieldVariables)
  {
    return NFieldVariables<std::vector<T>>::get(std::get<0>(fieldVariables));
  }
};

template<typename T>
struct NFieldVariables<std::vector<T>>
{
  static int get(std::vector<T> &fieldVariables)
  {
    VLOG(2) << "nFieldVariables of vector of " << StringUtility::demangle(typeid(T).name()) << ": " << fieldVariables.size() << " vector entries";

    int result = 0;
    for (int i = 0; i < fieldVariables.size(); i++)
    {
      int nVariables = NFieldVariables<T>::get(fieldVariables[i]);
      result += nVariables;
      VLOG(2) << "  entry " << i << ": " << nVariables;
    }
    VLOG(2) << "result: " << result;
    return result;
  }
};

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

  const int nFieldVariables = NFieldVariables<FieldVariablesForOutputWriterType>::get(fieldVariables);
  VLOG(2) << "buildPyFieldVariablesObject for " << nFieldVariables << " field variables: "
    << StringUtility::demangle(typeid(FieldVariablesForOutputWriterType).name());

  PyObject *pyData = PyList_New((Py_ssize_t)nFieldVariables);

  int fieldVariableIndex = 0;
  PythonLoopOverTuple::loopBuildPyFieldVariableObject<FieldVariablesForOutputWriterType>(fieldVariables, fieldVariableIndex, meshName, pyData, onlyNodalValues, mesh);

  return pyData;
}

}  // namespace
