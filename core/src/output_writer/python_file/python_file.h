#pragma once

#include <Python.h>  // has to be the first included header
#include <iostream>
#include <vector>

#include "control/types.h"
#include "output_writer/generic.h"
#include "data_management/finite_elements.h"

namespace OutputWriter
{

class PythonFile : public Generic
{
public:
  using Generic::Generic;
 
  //! write out solution to file, if timeStepNo is not -1, this value will be part of the filename
  template<typename DataType>
  void write(DataType &data, int timeStepNo = -1, double currentTime = -1);
  
private:
};

};  // namespace

#include "output_writer/python_file/python_file.tpp"
