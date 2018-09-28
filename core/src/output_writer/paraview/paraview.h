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
  Paraview(DihuContext context, PyObject *specificSettings);

  //! write out solution to given filename, if timeStepNo is not -1, this value will be part of the filename
  template<typename DataType>
  void write(DataType &data, int timeStepNo = -1, double currentTime = -1);

  //! set the rank subset of all ranks that collectively call write afterwards. This is needed for the "combineFiles" option.
  void setRankSubset(Partition::RankSubset &rankSubset);

  //! write the given field variable as VTK <DataArray> element to file, if onlyParallelDatasetElement write the <PDataArray> element
  template<typename FieldVariableType>
  static void writeParaviewFieldVariable(FieldVariableType &fieldVariable, std::ofstream &file,
                                         bool binaryOutput, bool fixedFormat, bool onlyParallelDatasetElement=false);


  //! write the a field variable indicating which ranks own which portion of the domain as VTK <DataArray> element to file, if onlyParallelDatasetElement write the <PDataArray> element
  template<typename FieldVariableType>
  static void writeParaviewPartitionFieldVariable(FieldVariableType &geometryField, std::ofstream &file,
                                                  bool binaryOutput, bool fixedFormat, bool onlyParallelDatasetElement=false);
  
  //! write a single *.vtp file that contains all data of all field variables. This is uses MPI I/O. It can be enabled with the "combineFiles" option.
  template<typename OutputFieldVariablesType>
  void writePolyDataFile(const OutputFieldVariablesType &fieldVariables);

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

protected:
  Partition::RankSubset rankSubset_; ///< the ranks that collectively call Paraview::write

  bool binaryOutput_;  ///< if the data output should be binary encoded using base64
  bool fixedFormat_;   ///< if non-binary output is selected, if the ascii values should be written with a fixed precision, like 1.000000e5

  bool combineFiles_;   ///< if the output data should be combined for 1D meshes into a single PolyData output file (*.vtp) and for 2D and 3D meshes to normal *.vtu,*.vts or *.vtr files. This is needed when the number of output files should be reduced.
};

};  // namespace

#include "output_writer/paraview/paraview.tpp"
