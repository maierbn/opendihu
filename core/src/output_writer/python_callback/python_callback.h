#pragma once

#include <Python.h>  // has to be the first included header
#include <iostream>
#include <vector>

#include "control/types.h"
#include "output_writer/generic.h"
#include "data_management/finite_element_method/finite_elements.h"

namespace OutputWriter
{

class PythonCallback : public Generic
{
public:

  //! constructor
  PythonCallback(DihuContext context, PythonConfig specificSettings);

  //! write out solution i.e. call the callback function in this case
  template<typename DataType>
  void write(DataType &data, int timeStepNo = -1, double currentTime = -1);

private:

  PyObject *callback_;    ///< the python callback function object
  bool onlyNodalValues_;  ///< if only nodal values should be output, this omits the derivative values for Hermite ansatz functions, for Lagrange functions it has no effect
};

};  // namespace

#include "output_writer/python_callback/python_callback.tpp"
