 #pragma once

#include <Python.h>  // has to be the first included header
#include <iostream>
#include <vector>

#include "control/types.h"
#include "output_writer/generic.h"
#include "data_management/finite_elements.h"

namespace OutputWriter
{

/**
 * Encapsulates method to write a vector to a numpy file *.npy
 */
class NumpyFileWriter
{
public:

  //! write data vector to a numpy file, data layout has shape given by nEntries
  static void writeToNumpyFile(std::vector<double> &data, std::string filename, std::vector<long> &nEntries);
};

/* fall-back implementation for not RegularFixed meshes or not Finite Elements data
 */
template<typename DataType>
class PythonStiffnessMatrixWriter : public NumpyFileWriter
{
public:
  static void writeNumpySolution(DataType &data, std::string filename);
};

/** specialization for finite elements, regular fixed mesh. For that write also stiffness matrix and rhs to numpy files
 */
template<int D, typename BasisFunctionType, typename Term>
class PythonStiffnessMatrixWriter<
  Data::FiniteElements<FunctionSpace::FunctionSpace<Mesh::StructuredRegularFixedOfDimension<D>,BasisFunctionType>,Term>
> :
  public NumpyFileWriter
{
public:
  typedef Data::FiniteElements<FunctionSpace::FunctionSpace<Mesh::StructuredRegularFixedOfDimension<D>,BasisFunctionType>,Term> DataType;

  //! write out solution to file filename as numpy array, as 1D data and in correct shape, also write rhs matrix
  static void writeNumpySolution(DataType &data, std::string filename);

private:
  //! write stiffness and rhs matrix
  static void writeMatrices(DataType &data, std::string filename);

};

};   // namespace

#include "output_writer/python_file/python_stiffness_matrix_writer.tpp"
