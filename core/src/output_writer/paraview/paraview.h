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

  //! write the given field variable as VTK <DataArray> element to file, if onlyParallelDatasetElement write the <PDataArray> element
  template<typename FieldVariableType>
  static void writeParaviewFieldVariable(FieldVariableType &fieldVariable, std::ofstream &file, bool binaryOutput, bool fixedFormat, bool onlyParallelDatasetElement);
  
  //! encode a Petsc vector in Base64
  static std::string encodeBase64(const Vec &vector);

  //! encode a std::vector as base64
  static std::string encodeBase64(const std::vector<double> &vector);

  //! encode a std::vector as base64
  static std::string encodeBase64(const std::vector<element_no_t> &vector);

  //! convert to a string with space separated values
  static std::string convertToAscii(const Vec &vector, bool humanReadable);
  
  //! convert to a string with space separated values
  static std::string convertToAscii(const std::vector<double> &vector, bool humanReadable);
  
  //! convert to a string with space separated values
  static std::string convertToAscii(const std::vector<element_no_t> &vector, bool humanReadable);
};

};  // namespace

#include "output_writer/paraview/paraview.tpp"
