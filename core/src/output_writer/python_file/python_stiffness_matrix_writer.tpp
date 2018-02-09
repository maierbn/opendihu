#include "output_writer/python_file/python_stiffness_matrix_writer.h"

#include <Python.h>  // has to be the first included header
#include <memory>
#include <sstream>

namespace OutputWriter
{
 
template<int D, typename BasisFunctionType>
void PythonStiffnessMatrixWriter<Data::FiniteElements<BasisOnMesh::BasisOnMesh<Mesh::RegularFixed<D>,BasisFunctionType>>>::
writeNumpySolution(Data::FiniteElements<BasisOnMesh::BasisOnMesh<Mesh::RegularFixed<D>,BasisFunctionType>> &data, std::string filename)
{
  LOG(TRACE) << "writeNumpySolution RegularFixed, D="<<D;
  
  // solution and rhs vectors in mesh shape

  // determine file names
  std::stringstream s[2];
  s[0] << filename << "_solution.npy";
  s[1] << filename << "_solution_shaped.npy";
  std::string filenameSolution = s[0].str();
  std::string filenameSolutionShaped = s[1].str();
  
  // get data of solution vector
  std::vector<double> vectorValues;
  PetscUtility::getVectorEntries(data.solution().values(), vectorValues);
  
  // determine number of entries in each D
  std::vector<long int> nEntries(D);
  for (int i=0; i<D; i++)
  {
    nEntries[i] = (data.mesh()->nElementsPerDimension(i) + 1) * data.nComponentsPerNode();
  }
  std::vector<long int> singleEntry({(long)vectorValues.size()});
  
  // write as numpy file
  writeToNumpyFile(vectorValues, filenameSolution, singleEntry);
  writeToNumpyFile(vectorValues, filenameSolutionShaped, nEntries);
  
  writeMatrices(data, filename);
}

template<int D, typename BasisFunctionType>
void PythonStiffnessMatrixWriter<Data::FiniteElements<BasisOnMesh::BasisOnMesh<Mesh::RegularFixed<D>,BasisFunctionType>>>::
writeMatrices(Data::FiniteElements<BasisOnMesh::BasisOnMesh<Mesh::RegularFixed<D>,BasisFunctionType>> &data, std::string filename)
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
  VecGetSize(data.rightHandSide().values(), &vectorSize);

  std::vector<int> indices(vectorSize);
  std::iota(indices.begin(), indices.end(), 0);
  std::vector<double> vectorValues(vectorSize);

  // determine number of entries in each D
  std::vector<long int> nEntries(D);
  for (int i=0; i<D; i++)
  {
    nEntries[i] = (data.mesh()->nElementsPerDimension(i) + 1);
  }
  std::vector<long int> singleEntry({(long)vectorValues.size()});
  
  VecGetValues(data.rightHandSide().values(), vectorSize, indices.data(), vectorValues.data());
  
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
  
  MatGetValues(data.stiffnessMatrix(), nRows, rowIndices.data(), nColumns, columnIndices.data(), matrixValues.data());
  
  // write as numpy file
  writeToNumpyFile(matrixValues, filenameStiffness, nEntries);
    
}

// implementation for not RegularFixed mesh
template<typename DataType>
void PythonStiffnessMatrixWriter<DataType>::
writeNumpySolution(DataType &data, std::string filename)
{
  // get PETSc vector of values
  std::vector<double> values;
  PetscUtility::getVectorEntries(data.solution().values(), values);
  
  // get number of entries
  std::vector<long> nEntries(1, values.size());
  
  // write to a file
  writeToNumpyFile(values, filename, nEntries);
}

};   // namespace
