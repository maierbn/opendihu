#pragma once

#include <Python.h>  // has to be the first included header
#include <iostream>
#include <vector>

namespace OutputWriter
{
 
template<typename BasisOnMeshType>
class PythonBase
{
public:
  //! create a python dict that contains data and meta data of field variables
  static PyObject *buildPyFieldVariablesObject(std::vector<std::shared_ptr<FieldVariable::FieldVariable<BasisOnMeshType>>> fieldVariables);
};

};  // namespace

#include "output_writer/python/python_base.tpp"
