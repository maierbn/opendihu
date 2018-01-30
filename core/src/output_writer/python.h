#pragma once

#include <Python.h>  // has to be the first included header
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
 
  //! write out solution to given filename, if timeStepNo is not -1, this value will be part of the filename
  template<typename DataType>
  void write(DataType &data, int timeStepNo = -1, double currentTime = -1);
 
private:
 
  //! write out solution templated by dimension 
  template <int dimension, typename DataType>
  void writeSolutionDim(DataType &data);
 
  //! write rhs matrix of FiniteElements solution
  template <int dimension, typename DataType> 
  void writeRhsMatrix(DataType &data);
 
  //! write data vector to a numpy file, data layout has shape given by nEntries and dimension
  void writeToNumpyFile(std::vector<double> &data, std::string filename, int dimension, std::vector<long> &nEntries);
  
  std::string filenameBase_;
};

};  // namespace

#include "output_writer/python.tpp"
