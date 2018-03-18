#include "output_writer/python_file/python_file.h"

#include <Python.h>  // has to be the first included header
#include <iostream>

#include "easylogging++.h"
#include "utility/python_utility.h"

namespace OutputWriter
{

PythonFile::PythonFile(PyObject *settings) : Generic(settings)
{
  onlyNodalValues_ = PythonUtility::getOptionBool(settings, "onlyNodalValues", true);
}

};
