#include "output_writer/python_callback/python_callback_writer.h"

#include <Python.h>   // this has to be the first included header
#include <iostream>

#include "easylogging++.h"
#include "output_writer/python/python.h"

namespace OutputWriter
{

template<typename BasisOnMeshType>
void PythonCallbackWriter<BasisOnMeshType>::
callCallback(PyObject *callback, std::vector<std::shared_ptr<FieldVariable::FieldVariable<BasisOnMeshType>>> fieldVariables, int timeStepNo, double currentTime)
{
  LOG(TRACE) << "callCallback timeStepNo="<<timeStepNo<<", currentTime="<<currentTime;
 
  if (callback == NULL || fieldVariables.empty())
  {
    if (callback == NULL)
      LOG(DEBUG) << "PythonCallbackWriter: no callback specified";
    return;
  }
  
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
  PyObject *data = Python<BasisOnMeshType>::buildPyDataObject(fieldVariables, timeStepNo, currentTime);
  //old signature: def callback(data, shape, nEntries, dimension, timeStepNo, currentTime)
  //PyObject *arglist = Py_BuildValue("(O,O,i,i,i,d)", dataList, nEntriesList, data.size(), nEntries.size(), timeStepNo_, currentTime_);
  
  //new signature def callback(data)
  PyObject *arglist = Py_BuildValue("(O)", data);
  
  // call callback function
  PyObject *returnValue = PyObject_CallObject(callback, arglist);
 
  // if there was an error while executing the function, print the error message
  if (returnValue == NULL)
    PyErr_Print();
  
  // decrement reference counters for python objects
  Py_DECREF(returnValue);
  Py_DECREF(arglist);
}

};