#include "output_writer/python_callback/python_callback.h"

#include <Python.h>  // has to be the first included header
#include <iostream>

#include "easylogging++.h"
#include "utility/python_utility.h"

namespace OutputWriter
{

PythonCallback::PythonCallback(DihuContext context, PyObject *settings) : Generic(context, settings)
{
  callback_ = PythonUtility::getOptionPyObject(settings, "callback");
  onlyNodalValues_ = PythonUtility::getOptionBool(settings, "onlyNodalValues", true);
}

};
