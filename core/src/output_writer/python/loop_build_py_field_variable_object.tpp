#include "output_writer/python/loop_build_py_field_variable_object.h"

#include <cstdlib>

namespace OutputWriter
{

namespace PythonLoopOverTuple
{
 
 /** Static recursive loop from 0 to number of entries in the tuple
 * Loop body
 */
template<typename OutputFieldVariablesType, int i>
inline typename std::enable_if<i < std::tuple_size<OutputFieldVariablesType>::value, void>::type
loopBuildPyFieldVariableObject(const OutputFieldVariablesType &fieldVariables, int &fieldVariableIndex, std::string meshName, 
                               PyObject *pyData, bool onlyNodalValues, std::shared_ptr<Mesh::Mesh> &mesh)
{
  // call what to do in the loop body
  if (buildPyFieldVariableObject<typename std::tuple_element<i,OutputFieldVariablesType>::type>(
       std::get<i>(fieldVariables), fieldVariableIndex, meshName, pyData, onlyNodalValues, mesh))
    return;
  
  // advance iteration to next tuple element
  loopBuildPyFieldVariableObject<OutputFieldVariablesType, i+1>(fieldVariables, fieldVariableIndex, meshName, pyData, onlyNodalValues, mesh);
}
 
// current element is of pointer type (not vector)
template<typename CurrentFieldVariableType>
typename std::enable_if<!TypeUtility::isTuple<CurrentFieldVariableType>::value && !TypeUtility::isVector<CurrentFieldVariableType>::value, bool>::type
buildPyFieldVariableObject(CurrentFieldVariableType currentFieldVariable, int &fieldVariableIndex, std::string meshName, 
                           PyObject *pyData, bool onlyNodalValues, std::shared_ptr<Mesh::Mesh> &mesh)
{
  // if mesh name is not the specified meshName step over this field variable but do not exit the loop over field variables
  if (currentFieldVariable->functionSpace()->meshName() != meshName)
  {
    return false;  // do not break iteration
  }
 
  // the first time retrieve the mesh
  if (fieldVariableIndex == 0)
  {
    mesh = currentFieldVariable->functionSpace();
  }
 
  VLOG(2) << "field variable \"" << currentFieldVariable->name() << "\"";
  VLOG(2) << "has " << currentFieldVariable->nComponents() << " components: " << currentFieldVariable->componentNames();

  const int nComponents = currentFieldVariable->nComponents();
  PyObject *pyComponents = PyList_New((Py_ssize_t)nComponents);

  // loop over components of field variable
  for (int componentNo = 0; componentNo < nComponents; componentNo++)
  {
    std::string componentName = currentFieldVariable->componentNames()[componentNo];
    VLOG(2) << "  component " << componentNo << " " << componentName;

    std::vector<double> values;
    currentFieldVariable->getValuesWithoutGhosts(componentNo, values, onlyNodalValues);

    VLOG(2) << "  values: " << values << ", values.size(): " << values.size();

    PyObject *pyValues = PythonUtility::convertToPythonList(values);
    VLOG(2) << " create pyComponent";
    PyObject *pyComponent = Py_BuildValue("{s s, s O}", "name", componentName.c_str(), "values", pyValues);

    VLOG(2) << " add to list";

    // add to list
    PyList_SetItem(pyComponents, (Py_ssize_t)componentNo, pyComponent);    // steals reference to pyComponent
  }

  VLOG(2) << "create pyFieldVariable";

  PyObject *pyFieldVariable = Py_BuildValue("{s s, s O}", "name", currentFieldVariable->name().c_str(), "components", pyComponents);

  // add to list
  PyList_SetItem(pyData, (Py_ssize_t)fieldVariableIndex, pyFieldVariable);    // steals reference to pyFieldVariable
  fieldVariableIndex++;

  return false;  // do not break iteration
}

// element i is of vector type
template<typename VectorType>
typename std::enable_if<TypeUtility::isVector<VectorType>::value, bool>::type
buildPyFieldVariableObject(VectorType currentFieldVariableVector, int &fieldVariableIndex, std::string meshName, 
                           PyObject *pyData, bool onlyNodalValues, std::shared_ptr<Mesh::Mesh> &mesh)
{
  for (auto& currentFieldVariable : currentFieldVariableVector)
  {
    // call function on all vector entries
    if (buildPyFieldVariableObject<typename VectorType::value_type>(currentFieldVariable, fieldVariableIndex, meshName, pyData, onlyNodalValues, mesh))
      return true;
  }
  
  return false;  // do not break iteration
}

// element i is of tuple type
template<typename TupleType>
typename std::enable_if<TypeUtility::isTuple<TupleType>::value, bool>::type
buildPyFieldVariableObject(TupleType currentFieldVariableTuple, int &fieldVariableIndex, std::string meshName, 
                           PyObject *pyData, bool onlyNodalValues, std::shared_ptr<Mesh::Mesh> &mesh)
{
  // call for tuple element
  loopBuildPyFieldVariableObject<TupleType>(currentFieldVariableTuple, fieldVariableIndex, meshName,
                                            pyData, onlyNodalValues, mesh);
  
  return false;  // do not break iteration
}

}  // namespace ExfileLoopOverTuple
}  // namespace OutputWriter
