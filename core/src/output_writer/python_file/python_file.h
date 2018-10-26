#pragma once

#include <Python.h>  // has to be the first included header
#include <iostream>
#include <vector>

#include "control/types.h"
#include "output_writer/generic.h"
#include "data_management/finite_element_method/finite_elements.h"

namespace OutputWriter
{

class PythonFile : public Generic
{
public:
  //! constructor
  PythonFile(DihuContext context, PyObject *specificSettings);

  //! write out solution to file, if timeStepNo is not -1, this value will be part of the filename
  template<typename DataType>
  void write(DataType &data, int timeStepNo = -1, double currentTime = -1);

private:

  //! open a python file stream using io module (python 3) or with low-level api (python 2.7)
  PyObject *openPythonFileStream(std::string filename, std::string writeFlag);

  //! write a python object to an already opened python file stream
  void outputPyObject(PyObject *file, PyObject *pyData);

  bool onlyNodalValues_;  ///< if only nodal values should be output, this omits the derivative values for Hermite ansatz functions, for Lagrange functions it has no effect
};

};  // namespace

#include "output_writer/python_file/python_file.tpp"
