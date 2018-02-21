#include "utility/petsc_utility.h"

#include <numeric>
#include <sstream>
#include <iomanip>

#include "easylogging++.h"
#include "petscksp.h"

void PetscUtility::getMatrixEntries(const Mat &matrix, std::vector<double> &matrixValues)
{
  int nRows, nColumns;
  MatGetSize(matrix, &nRows, &nColumns);
  
  std::vector<int> rowIndices(nRows);
  std::iota(rowIndices.begin(), rowIndices.end(), 0);
  std::vector<int> columnIndices(nColumns);
  std::iota(columnIndices.begin(), columnIndices.end(), 0);
  matrixValues.resize(nRows*nColumns);
  LOG(DEBUG) << "matrixValues contains " << nRows*nColumns << " entries for the " << nRows << "x" << nColumns << " matrix";
  
  MatGetValues(matrix, nRows, rowIndices.data(), nColumns, columnIndices.data(), matrixValues.data());
}

void PetscUtility::getVectorEntries(const Vec &vector, std::vector<double> &vectorValues)
{
  int nEntries;
  VecGetSize(vector, &nEntries);
  
  std::vector<int> indices(nEntries);
  std::iota(indices.begin(), indices.end(), 0);
  vectorValues.resize(nEntries);
  
  VecGetValues(vector, nEntries, indices.data(), vectorValues.data());
}

void PetscUtility::setVector(const std::vector<double> &vectorValues, Vec& vector)
{
  std::vector<int> indices(vectorValues.size());
  std::iota(indices.begin(), indices.end(), 0);
  VecSetValues(vector, vectorValues.size(), indices.data(), vectorValues.data(), INSERT_VALUES);
  
  VecAssemblyBegin(vector);
  VecAssemblyEnd(vector);
  
  int size;
  VecGetSize(vector, &size);
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

std::string PetscUtility::getStringMatrixVector(const Mat& matrix, const Vec& vector)
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

std::string PetscUtility::getStringMatrix(const Mat& matrix)
{
  int nRows, nColumns;
  MatGetSize(matrix, &nRows, &nColumns);
  
  std::vector<double> matrixValues, vectorValues;
  PetscUtility::getMatrixEntries(matrix, matrixValues);
  
  std::stringstream s;
  s<<std::endl<<"    ";
  for (int j=0; j<nColumns; j++)
  {
    s<<std::setw(5)<<std::setfill('_')<<j;
  }
  s<<std::string(5,'_');
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
    s<<std::string(5, ' ');
    s<<std::endl;
  }
  s<<std::endl;
  
  return s.str();
}

std::string PetscUtility::getStringVector(const Vec& vector)
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

std::string PetscUtility::getStringSparsityPattern(const Mat& matrix)
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

std::string PetscUtility::getStringConvergedReason(KSPConvergedReason convergedReason)
{
  
  // source: http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/KSP/KSPGetConvergedReason.html
  switch(convergedReason)
  {
  case KSP_CONVERGED_RTOL:
    return "KSP_CONVERGED_RTOL: residual 2-norm decreased by a factor of rtol, from 2-norm of right hand side";

  case KSP_CONVERGED_ATOL:
    return "KSP_CONVERGED_ATOL: residual 2-norm less than abstol";

  case KSP_CONVERGED_ITS:
    return "KSP_CONVERGED_ITS: used by the preonly preconditioner that always uses ONE iteration, or when the KSPConvergedSkip() convergence test routine is set.";

  case KSP_CONVERGED_CG_NEG_CURVE:
    return "KSP_CONVERGED_CG_NEG_CURVE";

  case KSP_CONVERGED_CG_CONSTRAINED:
    return "KSP_CONVERGED_CG_CONSTRAINED";

  case KSP_CONVERGED_STEP_LENGTH:
    return "KSP_CONVERGED_STEP_LENGTH";

  case KSP_CONVERGED_ITERATING:
    return "KSP_CONVERGED_ITERATING: returned if the solver is not yet finished";

  case KSP_DIVERGED_ITS:
    return "KSP_DIVERGED_ITS: required more than its to reach convergence";

  case KSP_DIVERGED_DTOL:
    return "KSP_DIVERGED_DTOL: residual norm increased by a factor of divtol";

  case KSP_DIVERGED_NANORINF:
    return "KSP_DIVERGED_NANORINF: residual norm became Not-a-number or Inf likely due to 0/0";

  case KSP_DIVERGED_BREAKDOWN:
    return "KSP_DIVERGED_BREAKDOWN: generic breakdown in method";

  case KSP_DIVERGED_BREAKDOWN_BICG:
    return "KSP_DIVERGED_BREAKDOWN_BICG: Initial residual is orthogonal to preconditioned initial residual. Try a different preconditioner, or a different initial Level.";
      
  case KSP_DIVERGED_NULL:
    return "KSP_DIVERGED_NULL";
      
  case KSP_DIVERGED_NONSYMMETRIC:
    return "KSP_DIVERGED_NONSYMMETRIC";
      
  case KSP_DIVERGED_INDEFINITE_PC:
    return "KSP_DIVERGED_INDEFINITE_PC";
      
  case KSP_DIVERGED_INDEFINITE_MAT:
    return "KSP_DIVERGED_INDEFINITE_MAT";
      
  case KSP_DIVERGED_PCSETUP_FAILED:
    return "KSP_DIVERGED_PCSETUP_FAILED";
    
  default:
    break;
  }
  
  std::stringstream s;
  if (convergedReason < 0)
    s << "divergence, ";
  else
    s << "converged, ";
  s << "unknown reason (" << int(convergedReason) << ")";
  return s.str();
}


