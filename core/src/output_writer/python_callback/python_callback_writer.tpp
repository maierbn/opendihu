#include "output_writer/python_callback/python_callback_writer.h"

#include <Python.h>   // this has to be the first included header
#include <iostream>

#include "easylogging++.h"
#include "output_writer/python/python.h"

namespace OutputWriter
{

template<typename BasisOnMeshType, typename OutputFieldVariablesType>
void PythonCallbackWriter<BasisOnMeshType,OutputFieldVariablesType>::
callCallback(PyObject *callback, OutputFieldVariablesType fieldVariables,
             int timeStepNo, double currentTime, bool onlyNodalValues)
{
  LOG(TRACE) << "callCallback timeStepNo="<<timeStepNo<<", currentTime="<<currentTime;

  if (callback == NULL)
  {
    LOG(DEBUG) << "PythonCallbackWriter: no callback specified";
    return;
  }

  // collect all available meshes
  std::set<std::string> meshNames;
  LoopOverTuple::loopCollectMeshNames<OutputFieldVariablesType>(fieldVariables, meshNames);

  PyObject *pyDataList = PyList_New((Py_ssize_t)meshNames.size());
  
  int meshIndex = 0;
  // loop over meshes and create an output file for each
  for (std::set<std::string>::iterator iter = meshNames.begin(); iter != meshNames.end(); iter++, meshIndex++)
  {
    std::string meshName = *iter;
    
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
    PyObject *pyData = Python<BasisOnMeshType,OutputFieldVariablesType>::buildPyDataObject(fieldVariables, meshName, timeStepNo, currentTime, onlyNodalValues);
    
    // set entry in list
    PyList_SetItem(pyDataList, (Py_ssize_t)meshIndex, pyData);    // steals reference to pyData
  }
  
  //old signature: def callback(data, shape, nEntries, dimension, timeStepNo, currentTime)
  //PyObject *arglist = Py_BuildValue("(O,O,i,i,i,d)", dataList, nEntriesList, data.size(), nEntries.size(), timeStepNo_, currentTime_);

  //new signature def callback([data0,data1,...])
  PyObject *pyArglist = Py_BuildValue("(O)", pyDataList);

  // call callback function
  PyObject *pyReturnValue = PyObject_CallObject(callback, pyArglist);

  // if there was an error while executing the function, print the error message
  if (pyReturnValue == NULL)
    PyErr_Print();

  // decrement reference counters for python objects
  Py_DECREF(pyReturnValue);
  Py_DECREF(pyArglist);
}

};