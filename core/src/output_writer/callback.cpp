#include "output_writer/callback.h"

#include <iostream>

#include "easylogging++.h"
#include <Python.h>

#include <utility/python_utility.h>
#include <utility/petsc_utility.h>
#include <mesh/regular_fixed.h>
#include <mesh/structured_deformable.h>
#include <mesh/mesh.h>

namespace OutputWriter
{

Callback::Callback(PyObject *settings) : Generic(settings)
{
  callback_ = PythonUtility::getOptionPyObject(settings, "callback");
}

void Callback::callCallback(std::vector<double> &data, std::vector<long> &nEntries, int timeStepNo, double currentTime)
{
  if (callback_ == NULL)
    return;
  
  // compose callback function
  PyObject *dataList = PythonUtility::convertToPythonList(data);
  PyObject *nEntriesList = PythonUtility::convertToPythonList(nEntries);
  
  // signature: def callback(data, shape, nEntries, dimension, timeStepNo, currentTime)
  PyObject *arglist = Py_BuildValue("(O,O,i,i,i,d)", dataList, nEntriesList, data.size(), nEntries.size(), timeStepNo, currentTime);
  PyObject *returnValue = PyObject_CallObject(callback_, arglist);
  
  // if there was an error while executing the function, print the error message
  if (returnValue == NULL)
    PyErr_Print();
  
  // decrement reference counters for python objects
  Py_DECREF(dataList);
  Py_DECREF(nEntriesList);
  Py_DECREF(returnValue);
  Py_DECREF(arglist);
}

};