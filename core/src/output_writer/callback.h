#pragma once

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
 
private:
 
  //! write out solution templated by dimension 
  template <int dimension>
  void writeSolutionDim(Data::Data &data, int timeStepNo, double currentTime);
  
  //! write out solution to given filename
  void writeSolution(Data::Data &data, int timeStepNo, double currentTime);
  
  //! call python callback function callback_ if it is set
  void callCallback(std::vector<double> &data, std::vector<long> &nEntries, int timeStepNo, double currentTime);
  
  PyObject *callback_;
};

};