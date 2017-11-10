#include "control/petsc_utility.h"

#include <numeric>
#include <sstream>
#include <iomanip>

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

void PetscUtility::setVector(std::vector< double >& vectorValues, Vec& vector)
{
  std::vector<int> indices(vectorValues.size());
  std::iota(indices.begin(), indices.end(), 0);
  VecSetValues(vector, vectorValues.size(), indices.data(), vectorValues.data(), INSERT_VALUES);
  
  VecAssemblyBegin(vector);
  VecAssemblyEnd(vector);
}

void PetscUtility::createVector(Vec& vector, int nEntries, std::string name)
{
  PetscErrorCode ierr;
  // create PETSc vector object
  ierr = VecCreate(PETSC_COMM_WORLD, &vector);  CHKERRV(ierr);
  
  if (name != "")
  {
    ierr = PetscObjectSetName((PetscObject)vector, name.c_str()); CHKERRV(ierr);
  }
  
  // initialize size of vector
  ierr = VecSetSizes(vector, PETSC_DECIDE, nEntries); CHKERRV(ierr);
  
  // set sparsity type and other options
  ierr = VecSetFromOptions(vector);  CHKERRV(ierr);
}

std::string PetscUtility::getStringMatrixVector(Mat& matrix, Vec& vector)
{
  std::string name;
  char *cName;
  PetscObjectGetName((PetscObject)vector, (const char **)&cName);
  name = cName;
  
  int nRows, nColumns;
  MatGetSize(matrix, &nRows, &nColumns);
  
  std::vector<double> matrixValues, vectorValues;
  PetscUtility::getMatrixEntries(matrix, matrixValues);
  PetscUtility::getVectorEntries(vector, vectorValues);
  
  std::stringstream s;
  s<<std::endl<<"    ";
  for (int j=0; j<nColumns; j++)
  {
    s<<std::setw(5)<<std::setfill('_')<<j;
  }
  s<<std::string(5,'_')<<" | "<<name;
  s<<std::endl;
  for (int i=0; i<nRows; i++)
  {
    s<<std::setw(3)<<std::setfill(' ')<<i<<"| ";
    for (int j=0; j<nColumns; j++)
    {
      if(matrixValues[i*nRows + j] == 0.0)
        s<<std::string(5, ' ');
      else
        s<<std::setw(4)<<std::setfill(' ')<<matrixValues[i*nRows + j]<<" ";
    }
    s<<std::string(5, ' ')<<"| "<<vectorValues[i];
    s<<std::endl;
  }
  s<<std::endl;
  
  return s.str();
}

std::string PetscUtility::getStringVector(Vec& vector)
{
  std::vector<double> vectorValues;
  PetscUtility::getVectorEntries(vector, vectorValues);
  
  int nEntries;
  VecGetSize(vector, &nEntries);
  
  std::stringstream s;
  for (int i=0; i<nEntries; i++)
  {
    s<<vectorValues[i]<<" ";
  }
  
  return s.str();
}

std::string PetscUtility::getStringSparsityPattern(Mat& matrix)
{
  int nRows, nColumns;
  MatGetSize(matrix, &nRows, &nColumns);
  
  std::vector<double> matrixValues;
  PetscUtility::getMatrixEntries(matrix, matrixValues);

  std::stringstream s;
  s<<" ";
  for (int j=0; j<nColumns; j++)
  {
    if (j%10 == 0)
      s<<"|";
    else if (j%2 == 0)
      s<<".";
    else 
      s<<" ";
  }
  s<<std::endl;
  for (int i=0; i<nRows; i++)
  {
    s<<" ";
    for (int j=0; j<nColumns; j++)
    {
      if(matrixValues[i*nRows + j] == 0.0)
        s<<" ";
      else
        s<<"*";
    }
    s<<std::endl;
  }
  s<<std::endl;
  return s.str();
}


