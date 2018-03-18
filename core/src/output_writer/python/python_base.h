#pragma once

#include <Python.h>  // has to be the first included header
#include <iostream>
#include <vector>

namespace OutputWriter
{
 
/** This is a base class for python output writer and contains functionality to create a python dictionary from some field variables
 * that will then be output.
 */
template<typename BasisOnMeshType>
class PythonBase
{
public:
  //! create a python dict that contains data and meta data of field variables
  //! @param onlyNodalValues: if only values at nodes should be contained, this discards the derivative values for Hermite
  static PyObject *buildPyFieldVariablesObject(std::vector<std::shared_ptr<FieldVariable::FieldVariable<BasisOnMeshType>>> fieldVariables, 
                                               bool onlyNodalValues);
};

};  // namespace

#include "output_writer/python/python_base.tpp"
