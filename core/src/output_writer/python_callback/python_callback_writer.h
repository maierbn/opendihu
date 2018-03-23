#pragma once

#include <Python.h>  // has to be the first included header
#include <iostream>
#include <vector>
#include <memory>

#include "basis_on_mesh/05_basis_on_mesh.h"

namespace OutputWriter
{
 
template<typename BasisOnMeshType, typename OutputFieldVariablesType>
class PythonCallbackWriter
{
public:
  //! call python callback
  static void callCallback(PyObject *callback, OutputFieldVariablesType fieldVariables, 
                           int timeStepNo, double currentTime, bool onlyNodalValues);  
};

};  // namespace

#include "output_writer/python_callback/python_callback_writer.tpp"