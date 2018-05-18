#include "utility/petsc_utility.h"

#include <numeric>
#include <sstream>
#include <iomanip>

#include "easylogging++.h"
#include "petscksp.h"

// color codes: https://github.com/shiena/ansicolor/blob/master/README.md
#define ANSI_COLOR_RED     "\x1b[31m"
#define ANSI_COLOR_GREEN   "\x1b[32m"
#define ANSI_COLOR_YELLOW  "\x1b[33m"
#define ANSI_COLOR_BLUE    "\x1b[34m"
#define ANSI_COLOR_MAGENTA "\x1b[35m"
#define ANSI_COLOR_CYAN    "\x1b[36m"
#define ANSI_COLOR_LIGHT_GRAY    "\x1b[90m"
#define ANSI_COLOR_LIGHT_WHITE    "\x1b[97m"
#define ANSI_COLOR_RESET   "\x1b[0m"

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

  // get values in row-major format
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

  const double zeroTolerance = 1e-15;

  std::stringstream s;
  s<<std::endl<<"    ";
  for (int j=0; j<nColumns; j++)
  {
    s<<std::setw(6)<<std::setfill('_')<<j;
  }
  s<<std::string(6,'_');
  s<<std::endl;
  for (int i=0; i<nRows; i++)
  {
    s<<std::setw(4)<<std::setfill(' ')<<i<<"| ";
    for (int j=0; j<nColumns; j++)
    {
      if(fabs(matrixValues[i*nRows + j]) <= zeroTolerance)
        s<<std::string(6, ' ');
      else
        s<<std::showpos<<std::setw(5)<<std::setfill(' ')<<std::setprecision(3)<<matrixValues[i*nRows + j]<<" ";
    }
    s<<std::string(6, ' ');
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

  LOG(DEBUG) << " getStringVector: " << nEntries;

  const double zeroTolerance = 1e-15;

  std::stringstream s;
  for (int i=0; i<nEntries; i++)
  {
    s << (fabs(vectorValues[i]) < zeroTolerance? 0.0 : vectorValues[i]) << " ";
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

std::string PetscUtility::getStringLinearConvergedReason(KSPConvergedReason convergedReason)
{


  // source: http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/KSP/KSPGetConvergedReason.html
  switch(convergedReason)
  {
  case KSP_CONVERGED_RTOL:
    return ANSI_COLOR_GREEN "KSP_CONVERGED_RTOL" ANSI_COLOR_RESET ": residual 2-norm decreased by a factor of rtol, from 2-norm of right hand side";

  case KSP_CONVERGED_ATOL:
    return ANSI_COLOR_GREEN "KSP_CONVERGED_ATOL" ANSI_COLOR_RESET ": residual 2-norm less than abstol";

  case KSP_CONVERGED_ITS:
    return ANSI_COLOR_GREEN "KSP_CONVERGED_ITS" ANSI_COLOR_RESET ": used by the preonly preconditioner that always uses ONE iteration, or when the KSPConvergedSkip() convergence test routine is set.";

  case KSP_CONVERGED_CG_NEG_CURVE:
    return ANSI_COLOR_GREEN "KSP_CONVERGED_CG_NEG_CURVE" ANSI_COLOR_RESET;

  case KSP_CONVERGED_CG_CONSTRAINED:
    return ANSI_COLOR_GREEN "KSP_CONVERGED_CG_CONSTRAINED" ANSI_COLOR_RESET;

  case KSP_CONVERGED_STEP_LENGTH:
    return ANSI_COLOR_GREEN "KSP_CONVERGED_STEP_LENGTH" ANSI_COLOR_RESET;

  case KSP_CONVERGED_ITERATING:
    return ANSI_COLOR_GREEN "KSP_CONVERGED_ITERATING" ANSI_COLOR_RESET ": returned if the solver is not yet finished";

  case KSP_DIVERGED_ITS:
    return ANSI_COLOR_RED "KSP_DIVERGED_ITS" ANSI_COLOR_RESET ": required more than its to reach convergence";

  case KSP_DIVERGED_DTOL:
    return ANSI_COLOR_RED "KSP_DIVERGED_DTOL" ANSI_COLOR_RESET ": residual norm increased by a factor of divtol";

  case KSP_DIVERGED_NANORINF:
    return ANSI_COLOR_RED "KSP_DIVERGED_NANORINF" ANSI_COLOR_RESET ": residual norm became Not-a-number or Inf likely due to 0/0";

  case KSP_DIVERGED_BREAKDOWN:
    return ANSI_COLOR_RED "KSP_DIVERGED_BREAKDOWN" ANSI_COLOR_RESET ": generic breakdown in method";

  case KSP_DIVERGED_BREAKDOWN_BICG:
    return ANSI_COLOR_RED "KSP_DIVERGED_BREAKDOWN_BICG" ANSI_COLOR_RESET ": Initial residual is orthogonal to preconditioned initial residual. Try a different preconditioner, or a different initial Level.";

  case KSP_DIVERGED_NULL:
    return ANSI_COLOR_RED "KSP_DIVERGED_NULL" ANSI_COLOR_RESET;

  case KSP_DIVERGED_NONSYMMETRIC:
    return ANSI_COLOR_RED "KSP_DIVERGED_NONSYMMETRIC" ANSI_COLOR_RESET;

  case KSP_DIVERGED_INDEFINITE_PC:
    return ANSI_COLOR_RED "KSP_DIVERGED_INDEFINITE_PC" ANSI_COLOR_RESET;

  case KSP_DIVERGED_INDEFINITE_MAT:
    return ANSI_COLOR_RED "KSP_DIVERGED_INDEFINITE_MAT" ANSI_COLOR_RESET;

  case KSP_DIVERGED_PCSETUP_FAILED:
    return ANSI_COLOR_RED "KSP_DIVERGED_PCSETUP_FAILED" ANSI_COLOR_RESET;

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

std::string PetscUtility::getStringNonlinearConvergedReason(SNESConvergedReason convergedReason)
{

  // source: http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/SNES/SNESConvergedReason.html#SNESConvergedReason
  switch(convergedReason)
  {
  case SNES_CONVERGED_FNORM_ABS:
    return ANSI_COLOR_GREEN "SNES_CONVERGED_FNORM_ABS" ANSI_COLOR_RESET ": ||F|| < atol";

  case SNES_CONVERGED_FNORM_RELATIVE:
    return ANSI_COLOR_GREEN "SNES_CONVERGED_FNORM_RELATIVE" ANSI_COLOR_RESET ": ||F|| < rtol*||F_initial||";

  case SNES_CONVERGED_SNORM_RELATIVE:
    return ANSI_COLOR_GREEN "SNES_CONVERGED_SNORM_RELATIVE" ANSI_COLOR_RESET ": Newton computed step size small; || delta x || < stol || x ||";

  case SNES_CONVERGED_ITS:
    return ANSI_COLOR_GREEN "SNES_CONVERGED_ITS" ANSI_COLOR_RESET ": maximum iterations reached";

  case SNES_CONVERGED_TR_DELTA:
    return ANSI_COLOR_GREEN " SNES_CONVERGED_TR_DELTA" ANSI_COLOR_RESET;

  case SNES_DIVERGED_FUNCTION_DOMAIN:
    return ANSI_COLOR_RED "SNES_DIVERGED_FUNCTION_DOMAIN:" ANSI_COLOR_RESET " the new x location passed the function is not in the domain of F";

  case SNES_DIVERGED_FUNCTION_COUNT:
    return ANSI_COLOR_RED "SNES_DIVERGED_FUNCTION_COUNT:" ANSI_COLOR_RESET " returned if the solver is not yet finished";

  case SNES_DIVERGED_LINEAR_SOLVE:
    return ANSI_COLOR_RED "SNES_DIVERGED_LINEAR_SOLVE:" ANSI_COLOR_RESET " the linear solve failed";

  case SNES_DIVERGED_FNORM_NAN:
    return ANSI_COLOR_RED "SNES_DIVERGED_FNORM_NAN" ANSI_COLOR_RESET;

  case SNES_DIVERGED_MAX_IT:
    return ANSI_COLOR_RED "SNES_DIVERGED_MAX_IT" ANSI_COLOR_RESET;

  case SNES_DIVERGED_LINE_SEARCH :
    return ANSI_COLOR_RED "SNES_DIVERGED_LINE_SEARCH" ANSI_COLOR_RESET ": the line search failed";

  case SNES_DIVERGED_INNER:
    return ANSI_COLOR_RED "SNES_DIVERGED_INNER" ANSI_COLOR_RESET ": inner solve failed";

  case SNES_DIVERGED_LOCAL_MIN:
    return ANSI_COLOR_RED "SNES_DIVERGED_LOCAL_MIN" ANSI_COLOR_RESET ": || J^T b || is small, implies converged to local minimum of F()";

  case SNES_CONVERGED_ITERATING:
    return "SNES_CONVERGED_ITERATING";

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


