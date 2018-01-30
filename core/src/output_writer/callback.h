#pragma once

#include <Python.h>  // has to be the first included header
#include <iostream>
#include <vector>

#include "control/types.h"
#include "output_writer/generic.h"
#include "data_management/finite_elements.h"

namespace OutputWriter
{

class Callback : public Generic
{
public:
 
  //! constructor
  Callback(PyObject *specificSettings);
 
  //! write out solution to given filename, if timeStepNo is not -1, this value will be part of the filename
  template<typename DataType>
  void write(DataType &data, int timeStepNo = -1, double currentTime = -1);
  
private:
 
  //! write out solution to given filename
  template<typename DataType>
  void writeSolution(DataType &data);
  
  //! write out solution to given filename
  template<int D, typename DataType>
  void writeSolutionDim(DataType &data);
  
  //! call python callback function callback_ if it is set
  void callCallback(std::vector<double> &data, std::vector<long> &nEntries);
  
  PyObject *callback_;
};

};  // namespace

#include "output_writer/callback.tpp"
