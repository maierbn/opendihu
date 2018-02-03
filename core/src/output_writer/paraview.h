#pragma once

#include <Python.h>  // has to be the first included header
#include <iostream>
#include <vector>

#include "control/types.h"
#include "output_writer/generic.h"

namespace OutputWriter
{   

class Paraview : public Generic
{
public:
 
  //! constructor
  Paraview(PyObject *specificSettings);
 
  //! write out solution to given filename, if timeStepNo is not -1, this value will be part of the filename
  template<typename DataType>
  void write(DataType &data, int timeStepNo = -1, double currentTime = -1);
  
private:
 
  //! write out solution templated by dimension 
  template <int dimension, typename DataType>
  void writeSolutionDim(DataType &data);
  
  //! write serial vtkRectilinearGrid file (structured, suffix *.vtr)
  template <typename Mesh, typename DataType>
  void writeRectilinearGrid(DataType& data);
 
  //! write serial vtkStructuredGrid file (structured, suffix *.vts)
  template <int dimension, typename DataType>
  void writeStructuredGrid(DataType& data);
  
  //! write serial vtkUnstructuredGrid file (unstructured, suffix *.vtu)
  template <int dimension, typename DataType>
  void writeUnstructuredGrid(DataType& data);
  
  //! encode a Petsc vector in Base64
  std::string encodeBase64(Vec &vector);
  
  //! encode a std::vector as base64
  std::string encodeBase64(std::vector<double> &vector);
  
  //! convert to a string with space separated values
  std::string convertToAscii(Vec &vector, bool humanReadable);
  std::string convertToAscii(std::vector<double> &vector, bool humanReadable);
};

};  // namespace

#include "output_writer/paraview.tpp"
