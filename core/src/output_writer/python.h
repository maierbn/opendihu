#pragma once

#include <iostream>
#include <vector>

#include "control/types.h"
#include "output_writer/generic.h"
#include "data_management/finite_elements.h"

namespace OutputWriter
{

class Python : public Generic
{
public:
 
  //! constructor
  Python(PyObject *specificSettings);
 
private:
 
  //! write out solution to given filename
  void writeSolution(Data::Data &data, int timeStepNo, double currentTime);
 
  //! write out solution templated by dimension 
  template <int dimension>
  void writeSolutionDim(Data::Data &data, int timeStepNo, double currentTime);
 
  //! write rhs matrix of FiniteElements solution
  template <int dimension> 
  void writeRhsMatrix(Data::FiniteElements &data);
 
  //! write data vector to a numpy file, data layout has shape given by nEntries and dimension
  void writeToNumpyFile(std::vector<double> &data, std::string filename, int dimension, std::vector<long> &nEntries);
  
  std::string filenameBase_;
};

};