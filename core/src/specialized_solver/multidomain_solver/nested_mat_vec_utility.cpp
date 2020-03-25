#include "specialized_solver/multidomain_solver/nested_mat_vec_utility.h"

#include <Python.h>  // has to be the first included header
#include <vector>
#include <algorithm>
#include <numeric>

#include "utility/vector_operators.h"
#include "utility/petsc_utility.h"
#include "easylogging++.h"

//#define USE_NESTED_MAT      // this disables this utility and directly uses the nested Vecs and Mats. This is faster in the handling of the Petsc variables but only GMRES solver is possible, no direct solvers.

namespace TimeSteppingScheme
{

namespace NestedMatVecUtility
{

//! from a Petsc Vec with nested type (nestedVec) create a new Petsc Vec (singleVec) that contains all values at once. If the singleVec already exists, do not create again, only copy the values.
void createVecFromNestedVec(Vec nestedVec, Vec &singleVec, std::shared_ptr<Partition::RankSubset> rankSubset)
{
#ifdef USE_NESTED_MAT
  singleVec = nestedVec;
  return;
#endif

  MPI_Comm mpiCommunicator = rankSubset->mpiCommunicator();
  PetscErrorCode ierr;

  // get the nested Vecs
  int nNestedVecs;
  Vec *nestedVecs;
  ierr = VecNestGetSubVecs(nestedVec, &nNestedVecs, &nestedVecs); CHKERRV(ierr);

  // if Vec object does not yet exist, create new one
  if (singleVec == PETSC_NULL)
  {

    // get the total sizes of all nested vecs
    int nEntriesGlobal = 0;
    int nEntriesLocal = 0;
    for (int nestedVecNo = 0; nestedVecNo < nNestedVecs; nestedVecNo++)
    {
      // sum up global size
      PetscInt nEntriesGlobalNestedVec = 0;
      ierr = VecGetSize(nestedVecs[nestedVecNo], &nEntriesGlobalNestedVec); CHKERRV(ierr);
      nEntriesGlobal += nEntriesGlobalNestedVec;
      
      // sum up local size
      PetscInt nEntriesLocalNestedVec = 0;
      ierr = VecGetLocalSize(nestedVecs[nestedVecNo], &nEntriesLocalNestedVec); CHKERRV(ierr);
      nEntriesLocal += nEntriesLocalNestedVec;
    }

    // get name of first subvec
    std::string name;
    char *cName;
    ierr = PetscObjectGetName((PetscObject)nestedVecs[0], (const char **)&cName); CHKERRV(ierr);
    name = cName;

    name += std::string("_singleVec");

    LOG(DEBUG) << "create single Vec \"" << name << "\" from " << nNestedVecs << " nested vecs, n entries global: " << nEntriesGlobal << ", local: " << nEntriesLocal;

    // create new Vec and assign name
    ierr = VecCreate(mpiCommunicator, &singleVec); CHKERRV(ierr);
    ierr = PetscObjectSetName((PetscObject) singleVec, name.c_str()); CHKERRV(ierr);

    // initialize size of vector
    ierr = VecSetSizes(singleVec, nEntriesLocal, nEntriesGlobal); CHKERRV(ierr);

    // set sparsity type and other options
    ierr = VecSetFromOptions(singleVec); CHKERRV(ierr);

  }

  // transfer the values from nestedVec to singleVec
  std::vector<int> indices;
  std::vector<double> values;

  // get global offset of the local portion in the single vec
  PetscInt singleVecGlobalBegin = 0;
  ierr = VecGetOwnershipRange(singleVec, &singleVecGlobalBegin, NULL); CHKERRV(ierr);

  // loop over subvecs
  for (int nestedVecNo = 0; nestedVecNo < nNestedVecs; nestedVecNo++)
  {
    Vec currentVec = nestedVecs[nestedVecNo];

    // get local range in the nested vec
    PetscInt valueNoGlobalBegin = 0;
    PetscInt valueNoGlobalEnd = 0;
    ierr = VecGetOwnershipRange(currentVec, &valueNoGlobalBegin, &valueNoGlobalEnd); CHKERRV(ierr);
    PetscInt nEntriesLocal = valueNoGlobalEnd - valueNoGlobalBegin;

    // get local values of the nested vec
    indices.resize(nEntriesLocal);
    std::iota(indices.begin(), indices.end(), valueNoGlobalBegin);

    values.resize(nEntriesLocal);
    ierr = VecGetValues(currentVec, indices.size(), indices.data(), values.data()); CHKERRV(ierr);
    
    // set values in the singleVec
    std::iota(indices.begin(), indices.end(), singleVecGlobalBegin);
    ierr = VecSetValues(singleVec, indices.size(), indices.data(), values.data(), INSERT_VALUES);

    singleVecGlobalBegin += nEntriesLocal;
  }

  // parallel assembly of single vec
  VecAssemblyBegin(singleVec);
  VecAssemblyEnd(singleVec);

  // dump vector singleVec to file
  if (VLOG_IS_ON(1))
  {
    // determine file name
    std::stringstream filename;
    // get name of first subvec
    std::string name;
    char *cName;
    ierr = PetscObjectGetName((PetscObject)nestedVecs[0], (const char **)&cName); CHKERRV(ierr);
    name = cName;

    filename << "out/" << name << "_createVecFromNestedVec";
    PetscUtility::dumpVector(filename.str(), std::string("matlab"), singleVec, mpiCommunicator);
  }
}

//! copy the values from a singleVec back to the nested Petsc Vec (nestedVec)
void fillNestedVec(Vec singleVec, Vec nestedVec)
{
#ifdef USE_NESTED_MAT
  singleVec = nestedVec;
  return;
#endif
  
  PetscErrorCode ierr;

  // get the nested Vecs
  int nNestedVecs;
  Vec *nestedVecs;
  ierr = VecNestGetSubVecs(nestedVec, &nNestedVecs, &nestedVecs); CHKERRV(ierr);

  // get all local values

  // get global offset of the local portion in the single vec
  PetscInt singleVecGlobalBegin = 0;
  PetscInt singleVecGlobalEnd = 0;
  ierr = VecGetOwnershipRange(singleVec, &singleVecGlobalBegin, &singleVecGlobalEnd); CHKERRV(ierr);

  PetscInt nEntriesLocal = singleVecGlobalEnd - singleVecGlobalBegin;

  std::vector<int> indices(nEntriesLocal);
  std::vector<double> values(nEntriesLocal);
  std::iota(indices.begin(), indices.end(), singleVecGlobalBegin);

  ierr = VecGetValues(singleVec, indices.size(), indices.data(), values.data()); CHKERRV(ierr);

  PetscInt valuesOffset = 0;

  // loop over subvecs
  for (int nestedVecNo = 0; nestedVecNo < nNestedVecs; nestedVecNo++)
  {
    Vec currentVec = nestedVecs[nestedVecNo];

    // get local range in the nested vec
    PetscInt valueNoGlobalBegin = 0;
    PetscInt valueNoGlobalEnd = 0;
    ierr = VecGetOwnershipRange(currentVec, &valueNoGlobalBegin, &valueNoGlobalEnd); CHKERRV(ierr);
    PetscInt nEntriesLocal = valueNoGlobalEnd - valueNoGlobalBegin;

    // set local values of the nested vec
    indices.resize(nEntriesLocal);
    std::iota(indices.begin(), indices.end(), valueNoGlobalBegin);

    ierr = VecSetValues(currentVec, indices.size(), indices.data(), values.data() + valuesOffset, INSERT_VALUES); CHKERRV(ierr);

    // parallel assembly of current vec
    VecAssemblyBegin(currentVec);
    VecAssemblyEnd(currentVec);

    valuesOffset += nEntriesLocal;
  }

  // dump vector nestedVec to file
  if (VLOG_IS_ON(1))
  {
    // determine file name
    std::stringstream filename;
    // get name of first subvec
    std::string name;
    char *cName;
    ierr = PetscObjectGetName((PetscObject)nestedVecs[0], (const char **)&cName); CHKERRV(ierr);
    name = cName;

    filename << "out/" << name << "_fillNestedVec";
    PetscUtility::dumpVector(filename.str(), std::string("matlab"), nestedVec, MPI_COMM_WORLD);
  }
}

//! from a Petsc Mat with nested type (nestedMat) create a new Petsc Mat (singleMat) that contains all values at once. If the singleMat already exists, do not create again, only copy the values.
void createMatFromNestedMat(Mat nestedMat, Mat &singleMat, std::shared_ptr<Partition::RankSubset> rankSubset)
{
#ifdef USE_NESTED_MAT
  singleMat = nestedMat;
  return;
#endif
  
  MPI_Comm mpiCommunicator = rankSubset->mpiCommunicator();
  int nRanks = rankSubset->size();
  int ownRankNo = rankSubset->ownRankNo();
  
  // collect sizes of nested Mats
  int nNestedMatRows;
  int nNestedMatColumns;
  PetscErrorCode ierr;

  // get the nested Mats as two-dimensional array, nestedMats[rowNo][columnNo]
  Mat **nestedMats;
  ierr = MatNestGetSubMats(nestedMat, &nNestedMatRows, &nNestedMatColumns, &nestedMats); CHKERRV(ierr);

  // get the sizes of all nested Mats and the global number of rows and columns
  std::vector<PetscInt> nRowsGlobalNestedMats(nNestedMatRows);
  std::vector<PetscInt> nRowsLocalNestedMats(nNestedMatRows);
  std::vector<PetscInt> nColumnsGlobalNestedMats(nNestedMatColumns);
  std::vector<PetscInt> nColumnsLocalNestedMats(nNestedMatColumns);

  int nRowsGlobal = 0;
  int nRowsLocal = 0;
  int nColumnsGlobal = 0;
  int nColumnsLocal = 0;

  // loop over rows of nested mat
  for (int nestedMatRowNo = 0; nestedMatRowNo < nNestedMatRows; nestedMatRowNo++)
  {
    for (int nestedMatColumnNo = 0; nestedMatColumnNo < nNestedMatColumns; nestedMatColumnNo++)
    {
      Mat currentMat = nestedMats[nestedMatRowNo][nestedMatColumnNo];

      if (currentMat != PETSC_NULL)
      {
        PetscInt nColumns = 0;
        ierr = MatGetSize(currentMat, &nRowsGlobalNestedMats[nestedMatRowNo], &nColumns); CHKERRV(ierr);
        ierr = MatGetLocalSize(currentMat, &nRowsLocalNestedMats[nestedMatRowNo], &nColumns); CHKERRV(ierr);
        
        nRowsGlobal += nRowsGlobalNestedMats[nestedMatRowNo];
        nRowsLocal += nRowsLocalNestedMats[nestedMatRowNo];

        break;
      }
    }
  }

  // loop over columns of nested mat
  for (int nestedMatColumnNo = 0; nestedMatColumnNo < nNestedMatColumns; nestedMatColumnNo++)
  {
    for (int nestedMatRowNo = 0; nestedMatRowNo < nNestedMatRows; nestedMatRowNo++)
    {
      Mat currentMat = nestedMats[nestedMatRowNo][nestedMatColumnNo];

      if (currentMat != PETSC_NULL)
      {

        PetscInt nRows = 0;
        ierr = MatGetSize(currentMat, &nRows, &nColumnsGlobalNestedMats[nestedMatColumnNo]); CHKERRV(ierr);
        ierr = MatGetLocalSize(currentMat, &nRows, &nColumnsLocalNestedMats[nestedMatColumnNo]); CHKERRV(ierr);
        
        nColumnsGlobal += nColumnsGlobalNestedMats[nestedMatColumnNo];
        nColumnsLocal += nColumnsLocalNestedMats[nestedMatColumnNo];

        break;
      }
    }
  }

  LOG(DEBUG) << "createMatFromNestedMat, nested mat has " << nNestedMatRows << "x" << nNestedMatColumns
    << " sub matrices, sizes global: " << nRowsGlobalNestedMats << " x " << nColumnsGlobalNestedMats
    << ", local: " << nRowsLocalNestedMats << " x " << nColumnsLocalNestedMats;

  // determine maximum number of nonzeros per row
  // loop over rows of sub matrices
  int nMaximumNonzerosPerRow = 0;
  for (int nestedMatRowNo = 0; nestedMatRowNo < nNestedMatRows; nestedMatRowNo++)
  {
    std::vector<PetscInt> nNonzeroEntriesRows;
    for (int nestedMatColumnNo = 0; nestedMatColumnNo < nNestedMatColumns; nestedMatColumnNo++)
    {
      Mat currentMat = nestedMats[nestedMatRowNo][nestedMatColumnNo];

      if (currentMat != PETSC_NULL)
      {
        PetscInt rowNoGlobalBegin = 0;
        PetscInt rowNoGlobalEnd = 0;
        ierr = MatGetOwnershipRange(currentMat, &rowNoGlobalBegin, &rowNoGlobalEnd); CHKERRV(ierr);

        if (nNonzeroEntriesRows.empty())
          nNonzeroEntriesRows.resize(rowNoGlobalEnd-rowNoGlobalBegin, 0);

        // loop over rows of mat
        for (int rowNoGlobal = rowNoGlobalBegin; rowNoGlobal < rowNoGlobalEnd; rowNoGlobal++)
        {
          PetscInt nNonzeroEntriesInRow;
          const PetscInt *columnIndices;
          const double *values;
          ierr = MatGetRow(currentMat, rowNoGlobal, &nNonzeroEntriesInRow, &columnIndices, &values); CHKERRV(ierr);
          ierr = MatRestoreRow(currentMat, rowNoGlobal, NULL, &columnIndices, &values); CHKERRV(ierr);

          int rowNoLocal = rowNoGlobal - rowNoGlobalBegin;
          nNonzeroEntriesRows[rowNoLocal] += nNonzeroEntriesInRow;
        }
      }
    }

    PetscInt nMaximumNonzerosPerRowCurrentMatRow = *std::max_element(nNonzeroEntriesRows.begin(), nNonzeroEntriesRows.end());
    nMaximumNonzerosPerRow = std::max(nMaximumNonzerosPerRow, nMaximumNonzerosPerRowCurrentMatRow);
  }
  
  LOG(DEBUG) << "nMaximumNonzerosPerRow: " << nMaximumNonzerosPerRow;

  // if Mat object does not yet exist, create new one
  if (singleMat == PETSC_NULL)
  {
    ierr = MatCreate(mpiCommunicator, &singleMat); CHKERRV(ierr);
    ierr = MatSetSizes(singleMat, nRowsLocal, nColumnsLocal, nRowsGlobal, nColumnsGlobal); CHKERRV(ierr);
    ierr = MatSetType(singleMat, MATAIJ); CHKERRV(ierr);
    
    // sparse matrix: preallocation of internal data structure
    // http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/Mat/MATAIJ.html#MATAIJ
    // MATAIJ = "aij" - A matrix type to be used for sparse matrices. This matrix type is identical to MATSEQAIJ when constructed with a single process communicator, and MATMPIAIJ otherwise.
    // As a result, for single process communicators, MatSeqAIJSetPreallocation is supported, and similarly MatMPIAIJSetPreallocation is supported for communicators controlling multiple processes.
    // It is recommended that you call both of the above preallocation routines for simplicity.
    PetscInt nNonZerosDiagonal = nMaximumNonzerosPerRow;
    PetscInt nNonZerosOffdiagonal = nMaximumNonzerosPerRow;
    ierr = MatSeqAIJSetPreallocation(singleMat, nNonZerosDiagonal, NULL); CHKERRV(ierr);
    ierr = MatMPIAIJSetPreallocation(singleMat, nNonZerosDiagonal, NULL, nNonZerosOffdiagonal, NULL); CHKERRV(ierr);
    
    // get name of first sub Mat
    std::string name;
    char *cName;
    ierr = PetscObjectGetName((PetscObject)nestedMats[0][0], (const char **)&cName); CHKERRV(ierr);
    name = cName;

    name += std::string("_singleMat");

    //ierr = MatSetFromOptions(singleMat); CHKERRV(ierr);         // either use MatSetFromOptions or MatSetUp to allocate internal data structures
    ierr = PetscObjectSetName((PetscObject)singleMat, name.c_str()); CHKERRV(ierr);
    
    LOG(DEBUG) << "create single Mat \"" << name << "\" with size local " << nRowsLocal << "x" << nColumnsLocal
      << ", global " << nRowsGlobal << "x" << nColumnsGlobal << ", nNonzeros: " << nMaximumNonzerosPerRow;
  }

  // set values
  PetscInt singleMatRowNoGlobalBegin = 0;
  ierr = MatGetOwnershipRange(singleMat, &singleMatRowNoGlobalBegin, NULL); CHKERRV(ierr);

  std::vector<PetscInt> indices;

  // determine column ranges of ranks of all nested mats
  std::vector<const PetscInt *> columnRangesNested(nNestedMatColumns);

  // loop over all nested mats
  for (int nestedMatColumnNo = 0; nestedMatColumnNo < nNestedMatColumns; nestedMatColumnNo++)
  {
    for (int nestedMatRowNo = 0; nestedMatRowNo < nNestedMatRows; nestedMatRowNo++)
    {
      Mat currentMat = nestedMats[nestedMatRowNo][nestedMatColumnNo];

      if (currentMat != PETSC_NULL)
      {
        ierr = MatGetOwnershipRangesColumn(currentMat, &columnRangesNested[nestedMatColumnNo]); CHKERRV(ierr);
        break;
      }
    }
  }

  // determine column ranges for the single mat
  const PetscInt *columnRangesGlobal;
  ierr = MatGetOwnershipRangesColumn(singleMat, &columnRangesGlobal); CHKERRV(ierr);

#ifndef NDEBUG
  std::stringstream s;
  for (int rankNo = 0; rankNo < nRanks; rankNo++)
  {
    s << columnRangesGlobal[rankNo] << ",";
  }
  LOG(DEBUG) << "columnRangesGlobal: " << s.str();
#endif

  PetscInt rowOffsetInSingleMatOwnDomain = 0;

  // loop over all nested mats
  for (int nestedMatRowNo = 0; nestedMatRowNo < nNestedMatRows; nestedMatRowNo++)
  {
    std::vector<PetscInt> columnOffsetInSingleMatOwnDomain(nRanks,0);
    for (int nestedMatColumnNo = 0; nestedMatColumnNo < nNestedMatColumns; nestedMatColumnNo++)
    {
      Mat currentMat = nestedMats[nestedMatRowNo][nestedMatColumnNo];

      if (currentMat != PETSC_NULL)
      {
        PetscInt rowNoGlobalBegin = 0;
        PetscInt rowNoGlobalEnd = 0;
        ierr = MatGetOwnershipRange(currentMat, &rowNoGlobalBegin, &rowNoGlobalEnd); CHKERRV(ierr);

        PetscInt nRowsLocal = rowNoGlobalEnd - rowNoGlobalBegin;

        // create index set indicating all rows of currentMat that are stored locally, in global numbering
        IS indexSetRows[1];
        ISCreateStride(mpiCommunicator, nRowsLocal, rowNoGlobalBegin, 1, &indexSetRows[0]);

        // create index set indicating all columns of currentMat that are stored locally (which are also all global columns), in global numbering
        IS indexSetColumns[1];
        ISCreateStride(mpiCommunicator, nColumnsGlobalNestedMats[nestedMatColumnNo], 0, 1, &indexSetColumns[0]);
        
        Mat *localSubMatrix;
        ierr = MatCreateSubMatrices(currentMat, 1, indexSetRows, indexSetColumns, MAT_INITIAL_MATRIX, &localSubMatrix); CHKERRV(ierr);

        PetscInt nRowsLocalSubMatrix = 0;
        PetscInt nColumnsLocalSubMatrix = 0;
        ierr = MatGetSize(localSubMatrix[0], &nRowsLocalSubMatrix, &nColumnsLocalSubMatrix); CHKERRV(ierr);
        
        if (VLOG_IS_ON(1))
        {
          std::stringstream filename;
          filename << "out/localSubMatrix" << nestedMatRowNo << "_" << nestedMatColumnNo << "_" << ownRankNo;
          PetscUtility::dumpMatrix(filename.str(), std::string("matlab"), localSubMatrix[0], MPI_COMM_SELF);
        }

        // loop over rows of mat
        for (PetscInt rowNoGlobal = rowNoGlobalBegin; rowNoGlobal < rowNoGlobalEnd; rowNoGlobal++)
        {
          PetscInt rowNoLocal = rowNoGlobal - rowNoGlobalBegin;
          
          // transfer row from nested mat to singleMat
          // get values of current row, this row will be split and the parts set at different locations
          PetscInt nNonzeroEntriesInRow;
          const PetscInt *columnIndices;
          const double *values;
          ierr = MatGetRow(localSubMatrix[0], rowNoLocal, &nNonzeroEntriesInRow, &columnIndices, &values); CHKERRV(ierr);

          PetscInt singleMatRowNoGlobal = singleMatRowNoGlobalBegin + rowOffsetInSingleMatOwnDomain + rowNoLocal;

#ifndef NDEBUG
          PetscInt begin = 0;
          PetscInt end = 0;
          ierr = MatGetOwnershipRange(singleMat, &begin, &end); CHKERRV(ierr);

          if (!(begin <= singleMatRowNoGlobal && singleMatRowNoGlobal < end))
          {
            LOG(FATAL) << "(" << nestedMatRowNo << "," << nestedMatColumnNo << ") set in row " << singleMatRowNoGlobal << " cols " << indices 
              << ", ownership: [" << begin << "," << end << "], nestedMatRowNo: " << nestedMatRowNo
              << ", nestedMatColumnNo: " << nestedMatColumnNo << ", rowNoGlobal: " << rowNoGlobal;

          }
          assert(begin <= singleMatRowNoGlobal);
          assert(singleMatRowNoGlobal < end);
#endif
          
          // loop over the ranges of ranks of the current retrieved row of currentMat, i.e. columnRangesNested[nestedMatColumnNo]
          // at the same time find out global column indices within this range, they are in [currentColumnIndexBegin,currentColumnIndexEnd]
          int currentColumnIndexBegin = 0;  
          int currentColumnIndexEnd = 0;
          int columnIndexNo = 0;
            
          for (int columnRangeNo = 0; columnRangeNo < nRanks; columnRangeNo++)
          {
            int currentColumnIndexRangeBegin = columnRangesNested[nestedMatColumnNo][columnRangeNo];
            int currentColumnIndexRangeEnd = columnRangesNested[nestedMatColumnNo][columnRangeNo+1];
          
            // find out end of current column range in columnIndices array
            for (; columnIndexNo < nNonzeroEntriesInRow; columnIndexNo++)
            {
              if (columnIndices[columnIndexNo] >= currentColumnIndexRangeEnd)
              {
                currentColumnIndexEnd = columnIndexNo;
                break;
              }
            }
            if (columnIndexNo == nNonzeroEntriesInRow)
              currentColumnIndexEnd = columnIndexNo;
            
            PetscInt nColumnIndicesInRange = currentColumnIndexEnd - currentColumnIndexBegin;

            // here, the columnIndices[i] are considered for i in [currentColumnIndexBegin,currentColumnIndexEnd]
            // the corresponding values will be stored to the global column block of the rank = columnRangeNo
            // and within this block to the sub part corresponding to the own part in the own block (columnRangesNested[columnRangeNo])

            PetscInt singleMatColumnNoGlobalBegin = columnRangesGlobal[columnRangeNo];

            // create global column indices for the singleMat
            std::stringstream s;
            indices.resize(nColumnIndicesInRange);
            for (int i = 0; i < nColumnIndicesInRange; i++)
            {
              PetscInt columnNoLocal = columnIndices[currentColumnIndexBegin + i] - currentColumnIndexRangeBegin;
              
              if (VLOG_IS_ON(1))
                s << columnNoLocal << ",";

              indices[i] = singleMatColumnNoGlobalBegin + columnOffsetInSingleMatOwnDomain[columnRangeNo] + columnNoLocal;
            }

            if (VLOG_IS_ON(1))
            {
              if (nestedMatRowNo == 1 && nestedMatColumnNo == 1)
              {
                VLOG(1) << "from (" << nestedMatRowNo << "," << nestedMatColumnNo << ")/(" << nNestedMatRows << "," << nNestedMatColumns << ") row " 
                  << rowNoGlobal << "/[" << rowNoGlobalBegin << "," << rowNoGlobalEnd << "]"
                  << " range " << columnRangeNo << "/" << nRanks << " (cols " << currentColumnIndexRangeBegin << "-" << currentColumnIndexRangeEnd-1 << ", "
                  << " indices from columnIndices[" << currentColumnIndexBegin << "-" << currentColumnIndexEnd-1 << "]/" << nNonzeroEntriesInRow 
                  << ": " << columnIndices[currentColumnIndexBegin] << "-" << columnIndices[currentColumnIndexEnd-1] << "), local: [" << s.str() << "] " 
                  << "set at global block " << columnRangeNo << "/" << nRanks << " starting at " << columnRangesGlobal[columnRangeNo] << " + local block " << nestedMatColumnNo 
                  << " (+" << columnOffsetInSingleMatOwnDomain[columnRangeNo] << "), begin: " << singleMatColumnNoGlobalBegin + columnOffsetInSingleMatOwnDomain[columnRangeNo]
                  << ", " << nColumnIndicesInRange << " indices: " << indices;
              }
            }

            ierr = MatSetValues(singleMat, 1, &singleMatRowNoGlobal, nColumnIndicesInRange, indices.data(), values + currentColumnIndexBegin, INSERT_VALUES);

            // columns in next columnRange start from end point of current columnRange
            currentColumnIndexBegin = currentColumnIndexEnd;
          }

          ierr = MatRestoreRow(localSubMatrix[0], rowNoLocal, NULL, &columnIndices, &values); CHKERRV(ierr);
        }

        ierr = MatDestroySubMatrices(1, &localSubMatrix); CHKERRV(ierr);
      }

      for (int rankNo = 0; rankNo < nRanks; rankNo++)
      {
        columnOffsetInSingleMatOwnDomain[rankNo] += columnRangesNested[nestedMatColumnNo][rankNo+1] - columnRangesNested[nestedMatColumnNo][rankNo];
      }
    }
    rowOffsetInSingleMatOwnDomain += nRowsLocalNestedMats[nestedMatRowNo];
  }

  LOG(DEBUG) << "final assembly";
 
  ierr = MatAssemblyBegin(singleMat, MAT_FINAL_ASSEMBLY); CHKERRV(ierr);
  ierr = MatAssemblyEnd(singleMat, MAT_FINAL_ASSEMBLY); CHKERRV(ierr);
}

} 
}  // namespace