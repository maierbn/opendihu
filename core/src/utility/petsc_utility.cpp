#include "utility/petsc_utility.h"

#include <numeric>
#include <sstream>
#include <iomanip>
#include <mpi.h>

#include "easylogging++.h"
#include "petscksp.h"
#include "utility/mpi_utility.h"
#include "control/dihu_context.h"   // for DihuContext::nRanksCommWorld
#include "output_writer/generic.h"

namespace PetscUtility
{

void getMatrixEntries(const Mat &matrix, std::vector<double> &matrixValues)
{
  MatType matrixType = MATSEQMAIJ;
  MatGetType(matrix, &matrixType);

  // if matrix is nested, only consider first submatrix
  if (std::string(matrixType) == std::string(MATNEST))
  {
    PetscInt nNestedRows = 1;
    PetscInt nNestedColumns = 1;
    Mat **nestedMats;
    MatNestGetSubMats(matrix,&nNestedRows,&nNestedColumns,&nestedMats);

    PetscInt nRows, nColumns;
    MatGetLocalSize(nestedMats[0][0], &nRows, &nColumns);

    std::vector<PetscInt> rowIndices(nRows);
    std::iota(rowIndices.begin(), rowIndices.end(), 0);
    std::vector<PetscInt> columnIndices(nColumns);
    std::iota(columnIndices.begin(), columnIndices.end(), 0);
    matrixValues.resize(nRows*nColumns, 0.0);

    // get values in row-major format with global indexing, note, there is no "MatGetValuesLocal"
    MatGetValues(nestedMats[0][0], nRows, rowIndices.data(), nColumns, columnIndices.data(), matrixValues.data());
  }
  else
  {
    // normal matrix
    PetscInt nRows, nColumns;
    MatGetLocalSize(matrix, &nRows, &nColumns);

    PetscInt rowBegin, rowEnd;
    MatGetOwnershipRange(matrix, &rowBegin, &rowEnd);
    std::vector<PetscInt> rowIndices(nRows);
    std::iota(rowIndices.begin(), rowIndices.end(), rowBegin);
    std::vector<PetscInt> columnIndices(nColumns);
    std::iota(columnIndices.begin(), columnIndices.end(), 0);
    matrixValues.resize(nRows*nColumns, 0.0);
    //LOG(DEBUG) << "matrixValues contains " << nRows*nColumns << " local entries for the " << nRows << "x" << nColumns << " local matrix";

    // get values in row-major format with global indexing, note, there is no "MatGetValuesLocal"
    MatGetValues(matrix, nRows, rowIndices.data(), nColumns, columnIndices.data(), matrixValues.data());
  }
}

void getVectorEntries(const Vec &vector, std::vector<double> &vectorValues)
{
  PetscInt ownershipBegin = 0;
  PetscInt ownershipEnd = 0;

  PetscErrorCode ierr;
  ierr = VecGetOwnershipRange(vector, &ownershipBegin, &ownershipEnd); CHKERRV(ierr);
  PetscInt nValues = ownershipEnd - ownershipBegin;

  std::vector<PetscInt> dofNosGlobal(nValues);
  std::iota(dofNosGlobal.begin(), dofNosGlobal.end(), ownershipBegin);
  vectorValues.resize(nValues);

  // get the values
  ierr = VecGetValues(vector, nValues, dofNosGlobal.data(), vectorValues.data()); CHKERRV(ierr);
}

void setVector(const std::vector<double> &vectorValues, Vec& vector)
{
  PetscInt ownershipBegin = 0;
  PetscInt ownershipEnd = 0;

  PetscErrorCode ierr;
  ierr = VecGetOwnershipRange(vector, &ownershipBegin, &ownershipEnd); CHKERRV(ierr);
  PetscInt nValues = ownershipEnd - ownershipBegin;
  assert (nValues >= vectorValues.size());

  std::vector<PetscInt> dofNosGlobal(nValues);
  std::iota(dofNosGlobal.begin(), dofNosGlobal.end(), ownershipBegin);

  VecSetValues(vector, vectorValues.size(), dofNosGlobal.data(), vectorValues.data(), INSERT_VALUES);

  VecAssemblyBegin(vector);
  VecAssemblyEnd(vector);
}

void createVector(Vec& vector, int nEntries, std::string name)
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

std::string getStringMatrixVector(const Mat& matrix, const Vec& vector)
{
#ifdef NDEBUG
  return std::string("");
#else
  std::string name;
  char *cName;
  PetscObjectGetName((PetscObject)vector, (const char **)&cName);
  name = cName;

  PetscInt nRows, nColumns;
  MatGetLocalSize(matrix, &nRows, &nColumns);

  std::vector<double> matrixValues, vectorValues;
  getMatrixEntries(matrix, matrixValues);
  getVectorEntries(vector, vectorValues);

  std::stringstream s;
  s << std::endl<< "    ";
  for (PetscInt j=0; j<nColumns; j++)
  {
    s << std::setw(5) << std::setfill('_') <<j;
  }
  s << std::string(5,'_') << " | " <<name;
  s << std::endl;
  for (PetscInt i=0; i<nRows; i++)
  {
    s << std::setw(3) << std::setfill(' ') <<i<< "| ";
    for (PetscInt j=0; j<nColumns; j++)
    {
      if (matrixValues[i*nRows + j] == 0.0)
        s << std::string(5, ' ');
      else
        s << std::setw(4) << std::setfill(' ') << matrixValues[i*nRows + j]<< " ";
    }
    s << std::string(5, ' ') << "| " <<vectorValues[i];
    s << std::endl;
  }
  s << std::endl;

  return s.str();
#endif
}

std::string getStringMatrix(const Mat& matrix)
{
#ifdef NDEBUG
  return std::string("");
#else

  PetscInt nNestedRows = 1;
  PetscInt nNestedColumns = 1;
  Mat **nestedMats;

  std::stringstream s;

  MatType matrixType = MATSEQMAIJ;
  MatGetType(matrix, &matrixType);

  if (std::string(matrixType) == std::string(MATNEST))
  {
    MatNestGetSubMats(matrix,&nNestedRows,&nNestedColumns,&nestedMats);
  }

  for (PetscInt i = 0; i < nNestedRows; i++)
  {
    for (PetscInt j = 0; j < nNestedColumns; j++)
    {
      const char *name;
      PetscObjectGetName((PetscObject)matrix, &name);
      s << "\"" << name << "\" (" << i << "," << j << ")/(" << nNestedRows << "," << nNestedColumns << "): ";

      Mat submatrix;
      if (std::string(matrixType) != std::string(MATNEST))
      {
        submatrix = matrix;
      }
      else
      {
        submatrix = nestedMats[i][j];
      }

      PetscInt nRows, nColumns;
      MatGetLocalSize(submatrix, &nRows, &nColumns);
      PetscInt nRowsGlobal, nColumnsGlobal;
      MatGetSize(submatrix, &nRowsGlobal, &nColumnsGlobal);

      std::vector<double> matrixValues;
      getMatrixEntries(submatrix, matrixValues);

      s << getStringMatrix(matrixValues, nRows, nColumns, nRowsGlobal, nColumnsGlobal);
      s << std::endl;
    }
  }
  return s.str();
#endif
}

std::string getStringMatrix(std::vector<double> &matrixValues, int nRows, int nColumns, int nRowsGlobal, int nColumnsGlobal)
{
#ifdef NDEBUG
  return std::string("");
#else
  const double zeroTolerance = 1e-15;

  std::stringstream s;

  PetscMPIInt nRanks = DihuContext::nRanksCommWorld();
  PetscMPIInt ownRank = DihuContext::ownRankNoCommWorld();
  //MPIUtility::handleReturnValue(MPI_Comm_size(MPI_COMM_WORLD, &nRanks), "MPI_Comm_size");

  if (nRanks > 1)
    s << ownRank << "/" << nRanks << ": ";
  s << nRows << "x" << nColumns;
  if (nRows != nRowsGlobal || nColumns != nColumnsGlobal)
    s << " (global: " << nRowsGlobal << "x" << nColumnsGlobal << ")";
  s << std::endl
    << "    ";
  for (PetscInt j=0; j<nColumns; j++)
  {
    s << std::setw(6) << std::setfill('_') <<j;
  }
  s << std::string(6,'_');
  s << std::endl;
  for (PetscInt i=0; i<nRows; i++)
  {
    s << std::setw(4) << std::setfill(' ') <<i<< "| ";
    for (PetscInt j=0; j<nColumns; j++)
    {
      if (fabs(matrixValues[i*nColumns + j]) <= zeroTolerance)
        s << std::string(6, ' ');
      else
        s << std::showpos << std::setw(5) << std::setfill(' ') << std::setprecision(3) << matrixValues[i*nRows + j]<< " " << std::noshowpos;
    }
    s << std::string(6, ' ');
    s << std::endl;
  }
  s << std::endl;

  return s.str();
#endif
}

std::string getStringVector(const Vec& vector)
{
#ifdef NDEBUG
  return std::string("");
#else
  std::vector<double> vectorValues;
  getVectorEntries(vector, vectorValues);

  PetscInt nEntries;
  VecGetLocalSize(vector, &nEntries);

  const double zeroTolerance = 1e-15;

  std::stringstream s;
  s << "(" << nEntries << " entries) ";
  for (PetscInt i=0; i<nEntries; i++)
  {
    s << std::setprecision(9) << (fabs(vectorValues[i]) < zeroTolerance? 0.0 : vectorValues[i]) << " ";
  }

  return s.str();
#endif
}

std::string getStringSparsityPattern(const Mat& matrix)
{
#ifdef NDEBUG
  return std::string("");
#else
  PetscInt nNestedRows = 1;
  PetscInt nNestedColumns = 1;
  Mat **nestedMats;

  std::stringstream s;

  MatType matrixType = MATSEQMAIJ;
  MatGetType(matrix, &matrixType);

  if (std::string(matrixType) == std::string(MATNEST))
  {
    MatNestGetSubMats(matrix,&nNestedRows,&nNestedColumns,&nestedMats);
  }

  for (PetscInt i = 0; i < nNestedRows; i++)
  {
    for (PetscInt j = 0; j < nNestedColumns; j++)
    {
      const char *name;
      PetscObjectGetName((PetscObject)matrix, &name);
      s << "\"" << name << "\" (" << i << "," << j << ")/(" << nNestedRows << "," << nNestedColumns << "): ";

      Mat submatrix;
      if (std::string(matrixType) != std::string(MATNEST))
      {
        submatrix = matrix;
      }
      else
      {
        submatrix = nestedMats[i][j];
      }

      PetscInt nRows, nColumns;
      MatGetLocalSize(submatrix, &nRows, &nColumns);
      PetscInt nRowsGlobal, nColumnsGlobal;
      MatGetSize(submatrix, &nRowsGlobal, &nColumnsGlobal);

      std::vector<double> matrixValues;
      getMatrixEntries(submatrix, matrixValues);

      s << std::endl << " ";
      for (PetscInt columnNo=0; columnNo<nColumns; columnNo++)
      {
        if (columnNo%10 == 0)
          s << "|";
        else if (columnNo%2 == 0)
          s << ".";
        else
          s << " ";
      }
      s << std::endl;
      for (PetscInt rowNo=0; rowNo<nRows; rowNo++)
      {
        s << " ";
        for (PetscInt columnNo=0; columnNo<nColumns; columnNo++)
        {
          if (fabs(matrixValues[rowNo*nColumns + columnNo]) < 1e-7)
            s << ".";
          else if (fabs(matrixValues[rowNo*nColumns + columnNo] - 1.0) < 1e-7)
            s << "1";
          else if (fabs(matrixValues[rowNo*nColumns + columnNo] - matrixValues[rowNo*nColumns + columnNo]) < 2e-4)
            s << "s";
          else
            s << "x";
        }
        s << std::endl;
      }
      s << std::endl;
    }
  }
  return s.str();
#endif
}

void checkDimensionsMatrixVector(Mat &matrix, Vec &input)
{
#ifndef NDEBUG
  PetscInt nRows, nColumns;
  MatGetSize(matrix, &nRows, &nColumns);
  PetscInt nEntries;
  VecGetSize(input, &nEntries);
  if (nColumns != nEntries)
  {
    LOG(ERROR) << "Matrix dimension " << nRows << "x" << nColumns << " does not match input vector (" << nEntries << ")!";
  }
  assert(nColumns == nEntries);
#endif
}

std::string getStringLinearConvergedReason(KSPConvergedReason convergedReason)
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

#ifndef __PGI

#if PETSC_VERSION_GE(3,11,0)
  case KSP_DIVERGED_PC_FAILED:
#else
  case KSP_DIVERGED_PCSETUP_FAILED:
#endif
    return ANSI_COLOR_RED "KSP_DIVERGED_PCSETUP_FAILED" ANSI_COLOR_RESET;
#endif

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

std::string getStringNonlinearConvergedReason(SNESConvergedReason convergedReason)
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

#if PETSC_VERSION_GE(3,12,0)
  case SNES_DIVERGED_TR_DELTA:
    return ANSI_COLOR_RED "SNES_DIVERGED_TR_DELTA" ANSI_COLOR_RESET;    
#else 
  case SNES_CONVERGED_TR_DELTA:
    return ANSI_COLOR_GREEN " SNES_CONVERGED_TR_DELTA" ANSI_COLOR_RESET;
#endif

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

void dumpVector(std::string filename, std::string format, Vec &vector, MPI_Comm mpiCommunicator, int componentNo, int nComponents)
{
  PetscViewer viewer;
  static int counter = 0;

  // set viewer format and suffix depending on format string
  std::string suffix = ".txt";
  PetscViewerFormat petscViewFormat = PETSC_VIEWER_DEFAULT;

  if (format == "matlab")
  {
    suffix = ".m";
    petscViewFormat = PETSC_VIEWER_ASCII_MATLAB;
  }
  else if (format == "ascii")
  {
    suffix = ".txt";
    petscViewFormat = PETSC_VIEWER_ASCII_INDEX;
  }
  else if (format != "default")
  {
    LOG(WARNING) << "In dumpVector: Format \"" << format << "\" is not allowed, possible values: default, ascii, matlab";
  }

  // compose output filename
  std::stringstream vectorOutputFilename;
  vectorOutputFilename << filename << "_" << std::setw(3) << std::setfill('0') << counter++ << suffix;

  std::stringstream filenameStream;
  filenameStream << vectorOutputFilename.str();

  // if there are multiple components, add component no to filename
  if (nComponents > 1)
  {
    filenameStream << componentNo;
  }

  // open file given by filename, create directory if necessary
  std::ofstream file;
  OutputWriter::Generic::openFile(file, filenameStream.str());
  file.close();

  // if the vector is nested, created single vector
  VecType vectorType = VECSEQ;
  VecGetType(vector, &vectorType);

  if (vectorType != NULL)
  {
    if (std::string(vectorType) == std::string(VECNEST))
    {
      PetscInt nNestedVecs = 1;
      Vec *nestedVecs;
      VecNestGetSubVecs(vector,&nNestedVecs,&nestedVecs);

      PetscInt nEntriesGlobal = 0;

      // collect global size
      PetscInt nEntries;
      for (PetscInt i = 0; i < nNestedVecs; i++)
      {
        VecGetSize(nestedVecs[i], &nEntries);

        nEntriesGlobal += nEntries;
      }

      int ownRankNo = 0;
      int nRanks = 0;
      MPI_Comm_rank(mpiCommunicator, &ownRankNo);
      MPI_Comm_size(mpiCommunicator, &nRanks);

      VLOG(1) << "ownRankNo: " << ownRankNo << ", nRanks: " << nRanks;

      Vec globalVector;

      if (ownRankNo == 0)
      {
        VecCreateSeq(MPI_COMM_SELF, nEntriesGlobal, &globalVector);
        const char *name = "v";
        //PetscObjectGetName((PetscObject)vector, &name);
        PetscObjectSetName((PetscObject)globalVector, name);
      }

      PetscErrorCode ierr;

      // gather all data from processes and set values in globalMatrix
      for (int j = 0; j < nNestedVecs; j++)
      {
        VecGetSize(nestedVecs[j], &nEntries);

        const PetscInt *ranges;
        VecGetOwnershipRanges(nestedVecs[j], &ranges);

        std::vector<double> values;

        for (int rankNo = 0; rankNo < nRanks; rankNo++)
        {
          PetscInt nRows = ranges[rankNo+1] - ranges[rankNo];
          values.resize(nRows, 0.0);

          //LOG(DEBUG) << "(" << j << "," << i << ") rank " << rankNo << " has " << nRows << " rows";

          if (ownRankNo == rankNo)
          {
            // get own data
            std::vector<PetscInt> rowIndices(nRows);
            std::iota(rowIndices.begin(), rowIndices.end(), ranges[rankNo]);

            ierr = VecGetValues(nestedVecs[j], nRows, rowIndices.data(), values.data()); CHKERRV(ierr);

            // send data to root
            if (rankNo != 0)
            {
              MPI_Send(values.data(), nRows, MPI_DOUBLE, 0, 0, mpiCommunicator);
            }

            VLOG(1) << nRows << " values [" << j << "] from rank " << rankNo << ", starting at " << ranges[rankNo] << ": " << values;
            //LOG(DEBUG) << "send " << nColumns*nRows << " values to 0";
          }

          if (ownRankNo == 0)
          {
            // receive data
            if (rankNo != 0)
            {
              MPI_Recv(values.data(), nRows, MPI_DOUBLE, rankNo, 0, mpiCommunicator, MPI_STATUS_IGNORE);
            }
            //LOG(DEBUG) << "receive " << nColumns*nRows << " values from " << rankNo;

            // set data
            std::vector<PetscInt> rowIndices(nRows);
            std::iota(rowIndices.begin(), rowIndices.end(), ranges[rankNo]);

            VLOG(1) << "insert at (" << ranges[rankNo] << ") " << nRows << " values [" << j << "] from rank " << rankNo << ": " << values;
            ierr = VecSetValues(globalVector, nRows, rowIndices.data(), values.data(), INSERT_VALUES); CHKERRV(ierr);
          }
        }
      }

      vector = globalVector;

      if (ownRankNo == 0)
      {
        mpiCommunicator = MPI_COMM_SELF;

        VecAssemblyBegin(vector);
        VecAssemblyEnd(vector);
      }
      else
      {
        mpiCommunicator = MPI_COMM_NULL;
      }
    }
  }

  // dump the data using a PetscViewer
  if (mpiCommunicator != MPI_COMM_NULL)
  {
    PetscErrorCode ierr;
    ierr = PetscViewerASCIIOpen(mpiCommunicator, filenameStream.str().c_str(), &viewer); CHKERRV(ierr);
    ierr = PetscViewerPushFormat(viewer, petscViewFormat); CHKERRV(ierr);
    ierr = VecView(vector, viewer); CHKERRV(ierr);
  }

  LOG(DEBUG) << "Vector written to \"" << filenameStream.str() << "\".";
}

//! dump the matrix to a file using PetscViewer, format is "default", "ascii" or "matlab"
void dumpMatrix(std::string filename, std::string format, Mat &matrix, MPI_Comm mpiCommunicator)
{
  PetscViewer viewer;
  static int counter = 0;

  // set viewer format and suffix depending on format string
  std::string suffix = ".txt";
  PetscViewerFormat petscViewFormat = PETSC_VIEWER_DEFAULT;

  if (format == "matlab")
  {
    suffix = ".m";
    petscViewFormat = PETSC_VIEWER_ASCII_MATLAB;
  }
  else if (format == "ascii")
  {
    suffix = ".txt";
    petscViewFormat = PETSC_VIEWER_ASCII_DENSE;
  }
  else if (format != "default")
  {
    LOG(WARNING) << "In dumpMatrix: Format \"" << format << "\" is not allowed, possible values: default, matlab";
  }

  // compose output filename
  std::stringstream matrixOutputFilename;
  matrixOutputFilename << filename << "_" << std::setw(3) << std::setfill('0') << counter++ << suffix;

  // open file given by filename, create directory if necessary
  std::ofstream file;
  OutputWriter::Generic::openFile(file, matrixOutputFilename.str());
  file.close();

  // if the matrix is nested, created single matrix
  MatType matrixType = MATSEQMAIJ;
  MatGetType(matrix, &matrixType);

  if (matrixType != NULL)
  {
    if (std::string(matrixType) == std::string(MATNEST))
    {
      PetscInt nNestedRows = 1;
      PetscInt nNestedColumns = 1;
      Mat **nestedMats;
      MatNestGetSubMats(matrix,&nNestedRows,&nNestedColumns,&nestedMats);

      PetscInt nRowsGlobal = 0;
      PetscInt nColumnsGlobal = 0;

      std::vector<PetscInt> nColumnsSubMats(nNestedColumns);
      std::vector<PetscInt> nRowsSubMats(nNestedRows);

      // collect global size
      // number of columns
      PetscInt nRows, nColumns;
      for (PetscInt i = 0; i < nNestedColumns; i++)
      {
        for (PetscInt j = 0; j < nNestedRows; j++)
        {
          if (nestedMats[j][i])
          {
            MatGetSize(nestedMats[j][i], &nRows, &nColumns);

            nColumnsSubMats[i] = nColumns;
            nColumnsGlobal += nColumns;
            break;
          }
        }
      }

      // number of rows
      for (PetscInt j = 0; j < nNestedRows; j++)
      {
        for (PetscInt i = 0; i < nNestedColumns; i++)
        {
          if (nestedMats[j][i])
          {
            MatGetSize(nestedMats[j][i], &nRows, &nColumns);

            nRowsSubMats[j] = nRows;
            nRowsGlobal += nRows;
            break;
          }
        }
      }

      int ownRankNo = 0;
      int nRanks = 0;
      MPI_Comm_rank(mpiCommunicator, &ownRankNo);
      MPI_Comm_size(mpiCommunicator, &nRanks);

      VLOG(1) << "ownRankNo: " << ownRankNo << ", nRanks: " << nRanks;

      Mat globalMatrix;

      if (ownRankNo == 0)
      {
        MatCreateSeqAIJ(MPI_COMM_SELF, nRowsGlobal, nColumnsGlobal, PETSC_DEFAULT, NULL, &globalMatrix);
        const char *name = "matrix";
        PetscObjectSetName((PetscObject)globalMatrix, name);
      }

      PetscErrorCode ierr;

      PetscInt globalRowOffset = 0;
      PetscInt globalColumnOffset = 0;

      // gather all data from processes and set values in globalMatrix
      for (PetscInt j = 0; j < nNestedRows; j++)
      {
        globalColumnOffset = 0;
        for (PetscInt i = 0; i < nNestedColumns; i++)
        {
          if (nestedMats[j][i])
          {

            MatGetSize(nestedMats[j][i], &nRows, &nColumns);

            const PetscInt *ranges;
            MatGetOwnershipRanges(nestedMats[j][i], &ranges);

            std::vector<double> values;

            for (int rankNo = 0; rankNo < nRanks; rankNo++)
            {
              PetscInt nRows = ranges[rankNo+1] - ranges[rankNo];
              values.resize(nColumns * nRows, 0.0);

              //LOG(DEBUG) << "(" << j << "," << i << ") rank " << rankNo << " has " << nRows << " rows";

              if (ownRankNo == rankNo)
              {
                // get own data
                std::vector<PetscInt> rowIndices(nRows);
                std::iota(rowIndices.begin(), rowIndices.end(), ranges[rankNo]);

                std::vector<PetscInt> columnIndices(nColumns);
                std::iota(columnIndices.begin(), columnIndices.end(), 0);

                ierr = MatGetValues(nestedMats[j][i], nRows, rowIndices.data(), nColumns, columnIndices.data(), values.data()); CHKERRV(ierr);

                // send data to root
                if (rankNo != 0)
                {
                  MPI_Send(values.data(), nColumns * nRows, MPI_DOUBLE, 0, 0, mpiCommunicator);
                }

                VLOG(1) << nColumns << "x" << nRows << " values [" << j << "," << i << "] from rank " << rankNo << ", starting at (" << ranges[rankNo] << ", 0): " << values;
                //LOG(DEBUG) << "send " << nColumns*nRows << " values to 0";
              }

              if (ownRankNo == 0)
              {
                // receive data
                if (rankNo != 0)
                {
                  MPI_Recv(values.data(), nColumns * nRows, MPI_DOUBLE, rankNo, 0, mpiCommunicator, MPI_STATUS_IGNORE);
                }
                //LOG(DEBUG) << "receive " << nColumns*nRows << " values from " << rankNo;

                // set data
                std::vector<PetscInt> rowIndices(nRows);
                std::iota(rowIndices.begin(), rowIndices.end(), globalRowOffset + ranges[rankNo]);

                std::vector<PetscInt> columnIndices(nColumns);
                std::iota(columnIndices.begin(), columnIndices.end(), globalColumnOffset);

                std::vector<double> nonZeroValues;
                nonZeroValues.reserve(nColumns);
                std::vector<PetscInt> nonZeroIndices;
                nonZeroIndices.reserve(nColumns);


                // extract non-zero entries
                for (int rowNo = 0; rowNo < nRows; rowNo++)
                {
                  int rowNoGlobal = globalRowOffset + ranges[rankNo] + rowNo;

                  nonZeroValues.clear();
                  nonZeroIndices.clear();
                  for (int columnNo = 0; columnNo < nColumns; columnNo++)
                  {
                    if (values[rowNo*nColumns + columnNo] != 0.0)
                    {
                      nonZeroValues.push_back(values[columnNo*nRows + rowNo]);
                      nonZeroIndices.push_back(globalColumnOffset + columnNo);
                    }
                  }

                  ierr = MatSetValues(globalMatrix, 1, &rowNoGlobal, nonZeroIndices.size(), nonZeroIndices.data(), nonZeroValues.data(), INSERT_VALUES); CHKERRV(ierr);
                }

                //VLOG(1) << "insert at (" << globalRowOffset + ranges[rankNo] << "," << globalColumnOffset << ") " << nRows << "x" << nColumns << " values [" << j << "," << i << "] from rank " << rankNo << ": " << values;
                //ierr = MatSetValues(globalMatrix, nRows, rowIndices.data(), nColumns, columnIndices.data(), values.data(), INSERT_VALUES); CHKERRV(ierr);
              }
            }
          }
          globalColumnOffset += nColumnsSubMats[i];
        }

        globalRowOffset += nRowsSubMats[j];
      }

      matrix = globalMatrix;

      if (ownRankNo == 0)
      {
        mpiCommunicator = MPI_COMM_SELF;

        MatAssemblyBegin(matrix, MAT_FINAL_ASSEMBLY);
        MatAssemblyEnd(matrix, MAT_FINAL_ASSEMBLY);
      }
      else
      {
        mpiCommunicator = MPI_COMM_NULL;
      }
    }
  }

  // dump the data using a PetscViewer
  if (mpiCommunicator != MPI_COMM_NULL)
  {
    PetscErrorCode ierr;
    ierr = PetscViewerASCIIOpen(mpiCommunicator, matrixOutputFilename.str().c_str(), &viewer); CHKERRV(ierr);
    ierr = PetscViewerPushFormat(viewer, petscViewFormat); CHKERRV(ierr);

    MatAssemblyBegin(matrix, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(matrix, MAT_FINAL_ASSEMBLY);

    ierr = MatView(matrix, viewer); CHKERRV(ierr);
  }

  LOG(DEBUG) << "Matrix written to \"" << matrixOutputFilename.str() << "\".";
}


}  // namespace
