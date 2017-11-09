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
  Python(PyObject *settings);
 
private:
 
  void writeSolution(Data::Data &data);
 
  template <int dimension>
  void writeSolutionDim(Data::Data &data);
 
  template <int dimension> 
  void writeRhsMatrix(Data::FiniteElements &data);
 
  void writeToNumpyFile(std::vector<double> &data, std::string filename, int dimension, std::vector<long> &nEntries);
  
  std::string filenameBase_;
  PyObject *settings_;
};

};