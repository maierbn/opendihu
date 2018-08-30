#pragma once

#include <Python.h>  // has to be the first included header
#include <iostream>
#include <vector>
#include <memory>

#include "function_space/function_space.h"

namespace OutputWriter
{

template<typename FunctionSpaceType, typename OutputFieldVariablesType>
class PythonCallbackWriter
{
public:
  //! call python callback
  static void callCallback(PyObject *callback, OutputFieldVariablesType fieldVariables,
                           int timeStepNo, double currentTime, bool onlyNodalValues);
};

};  // namespace

#include "output_writer/python_callback/python_callback_writer.tpp"