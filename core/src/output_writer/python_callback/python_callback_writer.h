#pragma once

#include <Python.h>  // has to be the first included header
#include <iostream>
#include <vector>
#include <memory>

#include "basis_on_mesh/05_basis_on_mesh.h"

namespace OutputWriter
{
 
template<typename BasisOnMeshType>
class PythonCallbackWriter
{
public:
  //! call python callback
  static void callCallback(PyObject *callback, std::vector<std::shared_ptr<FieldVariable::FieldVariable<BasisOnMeshType>>> fieldVariables, int timeStepNo, double currentTime);  
};

};  // namespace

#include "output_writer/python_callback/python_callback_writer.tpp"