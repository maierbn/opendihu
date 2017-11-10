#include "control/petsc_utility.h"

#include <numeric>

void PetscUtility::getMatrixEntries(Mat &matrix, std::vector<double> &matrixValues)
{
  int nRows, nColumns;
  MatGetSize(matrix, &nRows, &nColumns);
  
  std::vector<int> rowIndices(nRows);
  std::iota(rowIndices.begin(), rowIndices.end(), 0);
  std::vector<int> columnIndices(nColumns);
  std::iota(columnIndices.begin(), columnIndices.end(), 0);
  matrixValues.resize(nRows*nColumns);
  
  MatGetValues(matrix, nRows, rowIndices.data(), nColumns, columnIndices.data(), matrixValues.data());
}

void PetscUtility::getVectorEntries(Vec &vector, std::vector<double> &vectorValues)
{
  int nEntries;
  VecGetSize(vector, &nEntries);
  
  std::vector<int> indices(nEntries);
  std::iota(indices.begin(), indices.end(), 0);
  vectorValues.resize(nEntries);
  
  VecGetValues(vector, nEntries, indices.data(), vectorValues.data());
}