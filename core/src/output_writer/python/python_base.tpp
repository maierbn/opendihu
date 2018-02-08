#include "output_writer/python/python_base.h"

#include <Python.h>   // this has to be the first included header
#include <iostream>

#include "easylogging++.h"
#include "utility/python_utility.h"

namespace OutputWriter
{

template<typename BasisOnMeshType>
PyObject *PythonBase<BasisOnMeshType>::
buildPyFieldVariablesObject(std::vector<std::shared_ptr<FieldVariable::FieldVariable<BasisOnMeshType>>> fieldVariables)
{
  // build python dict containing field variables
  // [
  //   {"name" : "fieldVariableName",
  //    "components" : [
  //      {"name" : "componentName", "values": data},
  //    ]
  //   },
  // ]
  
 LOG(TRACE) << "buildPyFieldVariablesObject for " << fieldVariables.size() << " field variables";
 
  PyObject *pyData = PyList_New((Py_ssize_t)fieldVariables.size());
  
  // loop over field variables
  int fieldVariableNo = 0;
  for (typename std::vector<std::shared_ptr<FieldVariable::FieldVariable<BasisOnMeshType>>>::iterator fieldVariableIter = fieldVariables.begin();
       fieldVariableIter != fieldVariables.end(); fieldVariableIter++, fieldVariableNo++)
  {
    LOG(TRACE) << "field variable no " << fieldVariableNo;
    LOG(DEBUG) << "has "<<(*fieldVariableIter)->nComponents()<<" components: "<<(*fieldVariableIter)->componentNames();
   
    std::vector<std::string> componentNames = (*fieldVariableIter)->componentNames();
    
    PyObject *pyComponents = PyList_New((Py_ssize_t)(*fieldVariableIter)->nComponents());
    
    // loop over components of field variable 
    int componentNo = 0;
    for (std::vector<std::string>::iterator componentIter = componentNames.begin(); 
         componentIter != componentNames.end(); componentIter++, componentNo++)
    {
      LOG(DEBUG) << "  component " << componentNo << " [" << *componentIter << "]";
     
      std::vector<double> values;
      (*fieldVariableIter)->getValues(*componentIter, values);
      
      LOG(DEBUG) << "  values: " << values;
      
      PyObject *pyValues = PythonUtility::convertToPythonList(values);
      LOG(DEBUG) << " create pyComponent";
      PyObject *pyComponent = Py_BuildValue("{s s, s O}", "name", (*componentIter).c_str(), "values", pyValues);
      
      LOG(DEBUG) << " add to list";
      
      // add to list 
      PyList_SetItem(pyComponents, (Py_ssize_t)componentNo, pyComponent);    // steals reference to pyComponent 
    }
   
    LOG(DEBUG) << "create pyFieldVariable";
   
    PyObject *pyFieldVariable = Py_BuildValue("{s s, s O}", "name", (*fieldVariableIter)->name().c_str(), "components", pyComponents);
    
    LOG(DEBUG) << "add to list";
    
    // add to list 
    PyList_SetItem(pyData, (Py_ssize_t)fieldVariableNo, pyFieldVariable);    // steals reference to pyFieldVariable 
  }
  
  return pyData;
}
 
};