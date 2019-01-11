#pragma once

#include <Python.h>  // has to be the first included header
#include <iostream>
#include <vector>

namespace OutputWriter
{

/** This is a base class for python output writer and contains functionality to create a python dictionary from some field variables
 * that will then be output.
 */
template<typename OutputFieldVariablesType>
class PythonBase
{
public:
  //! create a python dict that contains data and meta data of field variables
  //! @param onlyNodalValues: if only values at nodes should be contained, this discards the derivative values for Hermite
  static PyObject *buildPyFieldVariablesObject(OutputFieldVariablesType fieldVariables, std::string meshName, bool onlyNodalValues, std::shared_ptr<Mesh::Mesh> &mesh);
};

} // namespace

#include "output_writer/python/python_base.tpp"
