#include "output_writer/python_file/python_stiffness_matrix_writer.h"

#include <Python.h>  // has to be the first included header
#include <memory>
#include <sstream>

#include "function_space/00_function_space_base_dim.h"

namespace OutputWriter
{

template<int D,typename BasisFunctionType,typename Term>
void PythonStiffnessMatrixWriter<Data::FiniteElements<FunctionSpace::FunctionSpace<Mesh::StructuredRegularFixedOfDimension<D>,BasisFunctionType>,Term>>::
writeNumpySolution(Data::FiniteElements<FunctionSpace::FunctionSpace<Mesh::StructuredRegularFixedOfDimension<D>,BasisFunctionType>,Term> &data, std::string filename)
{
  // This should dump debugging data, e.g. matrices in a numpy array format.
  // But because there are issues with linking to the numpy library, it is left out.

 /*
  LOG(TRACE) << "writeNumpySolution RegularFixed, D=" <<D;

  // solution and rhs vectors in mesh shape

  // determine file names
  std::stringstream s[2];
  s[0] << filename << "_solution.npy";
  s[1] << filename << "_solution_shaped.npy";
  std::string filenameSolution = s[0].str();
  std::string filenameSolutionShaped = s[1].str();

  // get data of solution vector
  std::vector<double> vectorValues;
  PetscUtility::getVectorEntries(data.solution()->values(), vectorValues);

  // determine number of entries in each D
  std::vector<long int> nEntries(D);
  for (int i=0; i<D; i++)
  {
    int averageNDofsPerElement1D = FunctionSpace::FunctionSpaceBaseDim<1,BasisFunctionType>::averageNDofsPerElement();
    dof_no_t dofsPerRow = (averageNDofsPerElement1D * data.mesh()->nElementsPerCoordinateDirectionLocal(i) + BasisFunctionType::nDofsPerNode());

    nEntries[i] = dofsPerRow * 1;
  }
  std::vector<long int> singleEntry({(long)vectorValues.size()});

  // write as numpy file
  writeToNumpyFile(vectorValues, filenameSolution, singleEntry);
  writeToNumpyFile(vectorValues, filenameSolutionShaped, nEntries);

  writeMatrices(data, filename);
  */
}

template<int D,typename BasisFunctionType,typename Term>
void PythonStiffnessMatrixWriter<Data::FiniteElements<FunctionSpace::FunctionSpace<Mesh::StructuredRegularFixedOfDimension<D>,BasisFunctionType>,Term>>::
writeMatrices(Data::FiniteElements<FunctionSpace::FunctionSpace<Mesh::StructuredRegularFixedOfDimension<D>,BasisFunctionType>,Term> &data, std::string filename)
{
  // solution and rhs vectors in mesh shape
  // determine file names
  std::stringstream s[3];
  s[0] << filename << "_rhs.npy";
  s[1] << filename << "_rhs_shaped.npy";
  s[2] << filename << "_stiffness.npy";
  std::string filenameRhs = s[0].str();
  std::string filenameRhsShaped = s[1].str();
  std::string filenameStiffness = s[2].str();

  // get data of rhs vector
  int vectorSize = 0;
  VecGetSize(data.rightHandSide()->valuesLocal(), &vectorSize);

  std::vector<int> indices(vectorSize);
  std::iota(indices.begin(), indices.end(), 0);
  std::vector<double> vectorValues(vectorSize);

  // determine number of entries in each D
  std::vector<long int> nEntries(D);
  for (int i=0; i<D; i++)
  {
    int averageNDofsPerElement1D = FunctionSpace::FunctionSpaceBaseDim<1,BasisFunctionType>::averageNDofsPerElement();
    dof_no_t dofsPerRow = (averageNDofsPerElement1D * data.functionSpace()->nElementsPerCoordinateDirectionLocal(i) + BasisFunctionType::nDofsPerNode());

    nEntries[i] = dofsPerRow;
  }
  std::vector<long int> singleEntry({(long)vectorValues.size()});

  VecGetValues(data.rightHandSide()->valuesLocal(), vectorSize, indices.data(), vectorValues.data());

  // write as numpy file
  writeToNumpyFile(vectorValues, filenameRhs, singleEntry);
  writeToNumpyFile(vectorValues, filenameRhsShaped, nEntries);

  // get stiffness matrix
  int nRows, nColumns;
  MatGetSize(data.stiffnessMatrix(), &nRows, &nColumns);
  std::vector<int> rowIndices(nRows);
  std::iota(rowIndices.begin(), rowIndices.end(), 0);
  std::vector<int> columnIndices(nColumns);
  std::iota(columnIndices.begin(), columnIndices.end(), 0);
  std::vector<double> matrixValues(nRows*nColumns);

  nEntries = {nRows, nColumns};

  MatGetValues(data.stiffnessMatrix().valuesLocal(), nRows, rowIndices.data(), nColumns, columnIndices.data(), matrixValues.data());

  // write as numpy file
  writeToNumpyFile(matrixValues, filenameStiffness, nEntries);

}

// implementation for not RegularFixed mesh
template<typename DataType>
void PythonStiffnessMatrixWriter<DataType>::
writeNumpySolution(DataType &data, std::string filename)
{
 // TODO: change to data.OutputFieldVariables instead of data.solution

 /*
  // determine file names
  std::stringstream s;
  s << filename << "_solution.npy";
  std::string filenameSolution = s.str();

  // get PETSc vector of values
  std::vector<double> values;
  PetscUtility::getVectorEntries(data.solution()->values(), values);

  // get number of entries
  std::vector<long> nEntries(1, values.size());

  // write to a file
  writeToNumpyFile(values, filenameSolution, nEntries);
  */
}

}  // namespace
