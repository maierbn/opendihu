#include "output_writer/python/python_base.h"

#include <Python.h>   // this has to be the first included header
#include <iostream>
#include <tuple>

#include "easylogging++.h"
#include "utility/python_utility.h"

namespace OutputWriter
{
 
namespace PythonLoopOverTuple
{
 
template<typename OutputFieldVariablesType, int i=0>
inline typename std::enable_if<i < std::tuple_size<OutputFieldVariablesType>::value, void>::type
loopOverTuple(const OutputFieldVariablesType &fieldVariables, PyObject *pyData, bool onlyNodalValues)
{
  // extract field variable
  typedef typename std::tuple_element<i,OutputFieldVariablesType>::type FieldVariableType;
  FieldVariableType fieldVariable = std::get<i>(fieldVariables);
  
  // create python object
  PyObject *pyFieldVariable = PythonBase<FieldVariableType>::
    buildPyFieldVariableObject(fieldVariable, onlyNodalValues);
  
  // add to list 
  PyList_SetItem(pyData, (Py_ssize_t)i, pyFieldVariable);    // steals reference to pyFieldVariable 
 
  loopOverTuple<OutputFieldVariablesType,i+1>(fieldVariables, pyData, onlyNodalValues);
}

}  // namespace

template<typename OutputFieldVariablesType>
PyObject *PythonBase<OutputFieldVariablesType>::
buildPyFieldVariablesObject(OutputFieldVariablesType fieldVariables, bool onlyNodalValues)
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
  
  PythonLoopOverTuple::loopOverTuple<OutputFieldVariablesType>(fieldVariables, pyData, onlyNodalValues);
  
  return pyData;
}
 
 
template<typename OutputFieldVariablesType>
template<typename FieldVariableType>
PyObject *PythonBase<OutputFieldVariablesType>::
buildPyFieldVariableObject(FieldVariableType fieldVariable, bool onlyNodalValues)
{
  VLOG(2) << "field variable";
  VLOG(2) << "has "<<fieldVariable->getNComponents()<<" components: "<<fieldVariable->componentNames();
  
  const int nComponents = fieldVariable->getNComponents();
  PyObject *pyComponents = PyList_New((Py_ssize_t)nComponents);
  
  // loop over components of field variable 
  for (int componentNo = 0; componentNo < nComponents; componentNo++)
  {
    std::string componentName = fieldVariable->componentNames()[componentNo];
    VLOG(2) << "  component " << componentNo << " " << componentName;
   
    std::vector<double> values;
    fieldVariable->getValues(componentNo, values, onlyNodalValues);
    
    VLOG(2) << "  values: " << values << ", values.size(): " << values.size();
    
    PyObject *pyValues = PythonUtility::convertToPythonList(values);
    VLOG(2) << " create pyComponent";
    PyObject *pyComponent = Py_BuildValue("{s s, s O}", "name", componentName.c_str(), "values", pyValues);
    
    VLOG(2) << " add to list";
    
    // add to list 
    PyList_SetItem(pyComponents, (Py_ssize_t)componentNo, pyComponent);    // steals reference to pyComponent 
  }
 
  VLOG(2) << "create pyFieldVariable";
 
  PyObject *pyFieldVariable = Py_BuildValue("{s s, s O}", "name", fieldVariable->name().c_str(), "components", pyComponents);
   
  return pyFieldVariable;
}

};