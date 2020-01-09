#include "partition/partitioned_petsc_mat/partitioned_petsc_mat_for_hyperelasticity.h"

template<typename DisplacementsFunctionSpaceType, typename PressureFunctionSpaceType, int nDisplacementComponents>
PartitionedPetscMatForHyperelasticity<DisplacementsFunctionSpaceType,PressureFunctionSpaceType,nDisplacementComponents>::
PartitionedPetscMatForHyperelasticity(
    std::shared_ptr<PartitionedPetscVecForHyperelasticity<DisplacementsFunctionSpaceType,PressureFunctionSpaceType,nDisplacementComponents>> partitionedPetscVecForHyperelasticity,
    int diagonalNonZeros, int offdiagonalNonZeros,
    std::string name
) : PartitionedPetscMatOneComponent<FunctionSpace::Generic>(

  // create generic function space with number of entries as in the given vector
  DihuContext::meshManager()->createGenericFunctionSpace(
    partitionedPetscVecForHyperelasticity->nEntriesLocal(), partitionedPetscVecForHyperelasticity->meshPartition(), std::string("genericMeshForMatrix")+name)->meshPartition(),

  diagonalNonZeros, offdiagonalNonZeros, name
),
partitionedPetscVecForHyperelasticity_(partitionedPetscVecForHyperelasticity)
{

/*
 PartitionedPetscMatOneComponent(std::shared_ptr<Partition::MeshPartition<FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>>> meshPartition,
                                  int diagonalNonZeros, int offdiagonalNonZeros, std::string name);
 */
}

template<typename DisplacementsFunctionSpaceType, typename PressureFunctionSpaceType, int nDisplacementComponents>
void PartitionedPetscMatForHyperelasticity<DisplacementsFunctionSpaceType,PressureFunctionSpaceType,nDisplacementComponents>::
setValue(int componentNoRow, PetscInt row, int componentNoColumn, PetscInt column, PetscScalar value, InsertMode mode)
{
  if (VLOG_IS_ON(2))
  {
    std::stringstream stream;
    stream << "\"" << this->name_ << "\" setValue " << (mode==INSERT_VALUES? "(insert)" : "(add)") << ", row " << row << ", column " << column << ", value " << value;
    VLOG(2) << stream.str();
  }

  // replace dirichlet BC values with the prescribed values
  if (partitionedPetscVecForHyperelasticity_->isPrescribed(componentNoRow, row)
      || partitionedPetscVecForHyperelasticity_->isPrescribed(componentNoColumn, column))
  {
    return;
  }

  assert(componentNoRow < nDisplacementComponents+1);
  assert(componentNoColumn < nDisplacementComponents+1);  //  < 4

  //assert(row < this->meshPartitionRows_->nDofsLocalWithGhosts());
  //assert(column < this->meshPartitionColumns_->nDofsLocalWithGhosts());

  // determine new indices
  row    = partitionedPetscVecForHyperelasticity_->nonBCDofNoGlobal(componentNoRow,    row);
  column = partitionedPetscVecForHyperelasticity_->nonBCDofNoGlobal(componentNoColumn, column);

  // this wraps the standard PETSc MatSetValue on the global matrix
  PetscErrorCode ierr;
  ierr = MatSetValues(this->globalMatrix_, 1, &row, 1, &column, &value, mode); CHKERRV(ierr);
}

template<typename DisplacementsFunctionSpaceType, typename PressureFunctionSpaceType, int nDisplacementComponents>
void PartitionedPetscMatForHyperelasticity<DisplacementsFunctionSpaceType,PressureFunctionSpaceType,nDisplacementComponents>::
dumpMatrixGlobalNatural(std::string filename)
{
  VLOG(1) << "dumpMatrixGlobalNatural, name: " << this->name_;

  // the local matrices have ordering
  // r0         r1         r2
  // ux uy uz p ux uy uz p ...

  // the final global matrix has ordering:
  // ux_0 ux_1 ux_2 ... uy_0 uy_1 uy_2 ... ... p_0 p_1 ...
  // \-----v-----/
  //    global natural

  std::shared_ptr<Partition::MeshPartition<DisplacementsFunctionSpaceType>> meshPartition =
    partitionedPetscVecForHyperelasticity_->meshPartition();

  std::shared_ptr<Partition::MeshPartition<PressureFunctionSpaceType>> meshPartitionPressure =
    partitionedPetscVecForHyperelasticity_->meshPartitionPressure();

  Mat matrix = this->globalMatrix_;

  std::stringstream result;
  int ownRankNo = meshPartition->ownRankNo();
  int nRanks = meshPartition->nRanks();

  int nRowsGlobal = 0;
  int nColumnsGlobal = 0;
  PetscErrorCode ierr;

  int nDisplacementDofsGlobal = meshPartition->nDofsGlobal();
  int nPressureDofsGlobal = meshPartitionPressure->nDofsGlobal();

  int nDisplacementDofsLocal = meshPartition->nDofsLocalWithoutGhosts();
  int nPressureDofsLocal = meshPartitionPressure->nDofsLocalWithoutGhosts();

  // get global matrix sizes
  nRowsGlobal = nDisplacementComponents*nDisplacementDofsGlobal + nPressureDofsGlobal;
  nColumnsGlobal = nRowsGlobal;

  int nRowsMatrixNonBc = 0;
  int nColumnsMatrixNonBc = 0;
  ierr = MatGetSize(matrix, &nRowsMatrixNonBc, &nColumnsMatrixNonBc); CHKERRV(ierr);


  VLOG(1) << "nRowsGlobal: " << nRowsGlobal << ", nRows actual matrix: " << nRowsMatrixNonBc;
  VLOG(1) << "nColumnsGlobal: " << nColumnsGlobal << ", nColumns actual matrix: " << nColumnsMatrixNonBc;

  // get displacement value columns
  const int *ownershipRanges;
  ierr = MatGetOwnershipRanges(matrix, &ownershipRanges); CHKERRV(ierr);

  int ownershipBegin = 0;
  int ownershipEnd = 0;
  ierr = MatGetOwnershipRange(matrix, &ownershipBegin, &ownershipEnd); CHKERRV(ierr);

  int nRowsMatrixNonBcLocal = ownershipEnd - ownershipBegin;
  VLOG(1) << "nRowsMatrixNonBcLocal: " << nRowsMatrixNonBcLocal
    << ", nDisplacementDofsLocal: " << nDisplacementDofsLocal << "/" << nDisplacementDofsGlobal << ", nPressureDofsLocal: " << nPressureDofsLocal << "/" << nPressureDofsGlobal;

  // get global matrix values

  // get values in row-major format with global indexing
  std::vector<double> matrixValuesLocal(nRowsMatrixNonBcLocal*nColumnsMatrixNonBc, 0.0);

  std::vector<int> rowIndices(nRowsMatrixNonBcLocal);
  std::vector<int> columnIndices(nColumnsMatrixNonBc);

  std::iota(rowIndices.begin(), rowIndices.end(), ownershipBegin);
  std::iota(columnIndices.begin(), columnIndices.end(), 0);

  ierr = MatGetValues(matrix, nRowsMatrixNonBcLocal, rowIndices.data(), nColumnsMatrixNonBc, columnIndices.data(), matrixValuesLocal.data());


  std::vector<double> matrixGlobalNonBc;
  if (ownRankNo == 0)
    matrixGlobalNonBc.resize(nRowsMatrixNonBc*nColumnsMatrixNonBc);

  std::vector<int> sizesOnRanks(nRanks);
  std::vector<int> offsets(nRanks);

  for (int rankNo = 0; rankNo < nRanks; rankNo++)
  {
    int nRowsOnRank = ownershipRanges[rankNo+1] - ownershipRanges[rankNo];
    sizesOnRanks[rankNo] = nRowsOnRank * nColumnsMatrixNonBc;

    // setup offsets for MPI_Gatherv
    if (rankNo == 0)
    {
      offsets[rankNo] = 0;
    }
    else
    {
      offsets[rankNo] = offsets[rankNo-1] + sizesOnRanks[rankNo-1];
    }
  }

  MPI_Gatherv(matrixValuesLocal.data(), nRowsMatrixNonBcLocal*nColumnsMatrixNonBc, MPI_DOUBLE,
              matrixGlobalNonBc.data(), sizesOnRanks.data(), offsets.data(), MPI_DOUBLE, 0, meshPartition->mpiCommunicator());

  // gather number of local dofs
  const int nComponents = nDisplacementComponents+1;
  std::vector<int> nDofsLocalRanks(nRanks*2);   // nDisplacementDofs, nPressureDofs for every rank

  std::vector<int> nDofsLocal{nDisplacementDofsLocal, nPressureDofsLocal};

  MPI_Gather(nDofsLocal.data(), 2, MPI_INT,
             nDofsLocalRanks.data(), 2, MPI_INT, 0, meshPartition->mpiCommunicator());

  VLOG(1) << "nDofsLocal: " << nDofsLocal << ", nDofsLocalRanks: " << nDofsLocalRanks;

  // gather dofNoLocalToDofNoNonBcGlobal_
  std::array<std::vector<std::vector<int>>,nComponents> dofNoLocalToDofNoNonBcGlobalRanks;
  std::array<std::vector<int>,nComponents> dofNoLocalToDofNoNonBcGlobalRanksBuffer;

  // and gather localDofNo to globalNaturalDofNo mappings
  std::array<std::vector<std::vector<global_no_t>>,nComponents> dofNosGlobalNaturalRanks;   // dofNosGlobalNaturalRanks[componentNo][rankNo][dofNoLocal];
  std::array<std::vector<global_no_t>,nComponents> dofNosGlobalNaturalRanksBuffer;

  for (int componentNo = 0; componentNo < nComponents; componentNo++)
  {
    // setup offsets for MPI_Gatherv
    std::vector<int> offsets(nRanks);
    std::vector<int> sizes(nRanks);

    int totalSize = 0;
    if (ownRankNo == 0)
    {
      offsets[0] = 0;

      if (componentNo < nDisplacementComponents)
      {
        sizes[0] = nDisplacementDofsLocal;
      }
      else
      {
        sizes[0] = nPressureDofsLocal;
      }
      totalSize = sizes[0];

      for (int rankNo = 1; rankNo < nRanks; rankNo++)
      {
        if (componentNo < nDisplacementComponents)
        {
          sizes[rankNo] = nDofsLocalRanks[2*rankNo+0];
        }
        else
        {
          sizes[rankNo] = nDofsLocalRanks[2*rankNo+1];
        }

        offsets[rankNo] = offsets[rankNo-1] + sizes[rankNo-1];
        totalSize += sizes[rankNo];
      }
    }
    else
    {
      if (componentNo < nDisplacementComponents)
      {
        sizes[ownRankNo] = nDisplacementDofsLocal;
      }
      else
      {
        sizes[ownRankNo] = nPressureDofsLocal;
      }
    }

    VLOG(1) << "component " << componentNo << ", sizes: " << sizes << ", offsets: " << offsets << ", totalSize: " << totalSize;

    // gather dofNoLocalToDofNoNonBcGlobal_
    assert(partitionedPetscVecForHyperelasticity_->dofNoLocalToDofNoNonBcGlobal()[componentNo].size() >= sizes[ownRankNo]);  // dofNoLocalToDofNoNonBcGlobal_ may contain ghost values

    dofNoLocalToDofNoNonBcGlobalRanksBuffer[componentNo].resize(totalSize);

    VLOG(1) << "MPI_Gatherv, send " << sizes[ownRankNo] << " ints, sizes: " << sizes << ", offsets: " << offsets << ", totalSize of recv buffer: " << totalSize;
    MPI_Gatherv(partitionedPetscVecForHyperelasticity_->dofNoLocalToDofNoNonBcGlobal()[componentNo].data(), sizes[ownRankNo], MPI_INT,
                dofNoLocalToDofNoNonBcGlobalRanksBuffer[componentNo].data(), sizes.data(), offsets.data(),
                MPI_INT, 0, meshPartition->mpiCommunicator());

    std::vector<global_no_t> dofNosGlobalNatural;

    if (componentNo < nDisplacementComponents)
    {
      meshPartition->getDofNosGlobalNatural(dofNosGlobalNatural);
    }
    else
    {
      meshPartitionPressure->getDofNosGlobalNatural(dofNosGlobalNatural);
    }
    assert(dofNosGlobalNatural.size() == sizes[ownRankNo]);

    dofNosGlobalNaturalRanksBuffer[componentNo].resize(totalSize);

    MPI_Gatherv(dofNosGlobalNatural.data(), sizes[ownRankNo], MPI_UNSIGNED_LONG_LONG,
                dofNosGlobalNaturalRanksBuffer[componentNo].data(), sizes.data(), offsets.data(),
                MPI_UNSIGNED_LONG_LONG, 0, meshPartition->mpiCommunicator());

    // copy the data, split by rankNo
    if (ownRankNo == 0)
    {
      dofNosGlobalNaturalRanks[componentNo].resize(nRanks);
      dofNoLocalToDofNoNonBcGlobalRanks[componentNo].resize(nRanks);

      for (int rankNo = 0; rankNo < nRanks; rankNo++)
      {
        dofNosGlobalNaturalRanks[componentNo][rankNo].resize(sizes[rankNo]);
        dofNoLocalToDofNoNonBcGlobalRanks[componentNo][rankNo].resize(sizes[rankNo]);

        for (int dofNoLocal = 0; dofNoLocal < sizes[rankNo]; dofNoLocal++)
        {
          dofNosGlobalNaturalRanks[componentNo][rankNo][dofNoLocal] = dofNosGlobalNaturalRanksBuffer[componentNo][offsets[rankNo] + dofNoLocal];
          dofNoLocalToDofNoNonBcGlobalRanks[componentNo][rankNo][dofNoLocal] = dofNoLocalToDofNoNonBcGlobalRanksBuffer[componentNo][offsets[rankNo] + dofNoLocal];
        }
      }
    }
  }

  if (ownRankNo == 0)
  {
    if (VLOG_IS_ON(1))
    {
      VLOG(1) << "dofNosGlobalNaturalRanks: " << dofNosGlobalNaturalRanks;
      VLOG(1) << "dofNoLocalToDofNoNonBcGlobalRanks: " << dofNoLocalToDofNoNonBcGlobalRanks;

      if (nRowsMatrixNonBc < 200)
      {
        std::stringstream s;
        std::stringstream s2;
        for (int j = 0; j < nRowsMatrixNonBc; j++)
        {
          for (int i = 0; i < nColumnsMatrixNonBc; i++)
          {
            s << matrixGlobalNonBc[j*nColumnsMatrixNonBc + i] << ", ";
            if (fabs(matrixGlobalNonBc[j*nColumnsMatrixNonBc + i]) < 1e-7)
              s2 << ".";
            else if (fabs(matrixGlobalNonBc[j*nColumnsMatrixNonBc + i] - 1.0) < 1e-7)
              s2 << "1";
            else if (fabs(matrixGlobalNonBc[j*nColumnsMatrixNonBc + i] - matrixGlobalNonBc[i*nColumnsMatrixNonBc + j]) < 2e-4)
              s2 << "s";
            else
              s2 << "x";
          }
          s << std::endl;
          s2 << std::endl;
        }
        VLOG(1) << "matrixGlobalNonBc: (" << nRowsMatrixNonBc << "x" << nColumnsMatrixNonBc << "):" << std::endl << s.str();
        // LOG(DEBUG) << "matrixGlobalNonBc: (" << nRowsMatrixNonBc << "x" << nColumnsMatrixNonBc << ") nonzeros:" << std::endl << s2.str();
      }
    }

    // copy matrixGlobalNonBc values to right place
    // old values in matrixGlobalNonBc are in dofNoNonBcGlobal indexing, new values in matrixGlobalNatural will be in global natural ordering
    std::vector<double> matrixGlobalNatural(nRowsGlobal*nColumnsGlobal, 0.0);

    // set to identity
    for (int i = 0; i < nRowsGlobal; i++)
      matrixGlobalNatural[i*nColumnsGlobal + i] = 1.0;

    // copy entries of old matrix to new one
    int rowNo = 0;
    for (int rankNoRow = 0; rankNoRow < nRanks; rankNoRow++)
    {
      for (int componentNoRow = 0; componentNoRow < nComponents; componentNoRow++)
      {
        int nDofsLocalRankRow = nDofsLocalRanks[2*rankNoRow + 0];
        if (componentNoRow == nDisplacementComponents)   // pressure value
        {
          nDofsLocalRankRow = nDofsLocalRanks[2*rankNoRow + 1];
        }

        for (int dofNoLocalRow = 0; dofNoLocalRow < nDofsLocalRankRow; dofNoLocalRow++, rowNo++)
        {
          int dofNoNonBcGlobalRow = dofNoLocalToDofNoNonBcGlobalRanks[componentNoRow][rankNoRow][dofNoLocalRow];

          VLOG(2) << "row [rank,comp,dof] = [" << rankNoRow << "," << componentNoRow << "," << dofNoLocalRow << "], nonBcGlobal: " << dofNoNonBcGlobalRow;

          // if dof is dirichlet bc, do not set any values (row is already set to identity matrix row)
          if (dofNoNonBcGlobalRow == -1)
            continue;

          int dofNoGlobalNaturalRow = componentNoRow*nDisplacementDofsGlobal + dofNosGlobalNaturalRanks[componentNoRow][rankNoRow][dofNoLocalRow];

          VLOG(2) << "    globalNatural: " << dofNoGlobalNaturalRow;

          int columnNo = 0;
          for (int rankNoColumn = 0; rankNoColumn < nRanks; rankNoColumn++)
          {
            for (int componentNoColumn = 0; componentNoColumn < nComponents; componentNoColumn++)
            {
              int nDofsLocalRankColumn = nDofsLocalRanks[2*rankNoColumn + 0];
              if (componentNoColumn == nDisplacementComponents)   // pressure value
              {
                nDofsLocalRankColumn = nDofsLocalRanks[2*rankNoColumn + 1];
              }

              for (int dofNoLocalColumn = 0; dofNoLocalColumn < nDofsLocalRankColumn; dofNoLocalColumn++, columnNo++)
              {
                int dofNoNonBcGlobalColumn = dofNoLocalToDofNoNonBcGlobalRanks[componentNoColumn][rankNoColumn][dofNoLocalColumn];

                VLOG(2) << "column [rank,comp,dof] = [" << rankNoColumn << "," << componentNoColumn << "," << dofNoLocalColumn << "], nonBcGlobal: " << dofNoNonBcGlobalColumn;

                if (dofNoNonBcGlobalColumn == -1)
                  continue;

                int dofNoGlobalNaturalColumn = componentNoColumn*nDisplacementDofsGlobal +dofNosGlobalNaturalRanks[componentNoColumn][rankNoColumn][dofNoLocalColumn];

                VLOG(2) << "    globalNatural: " << dofNoGlobalNaturalColumn;
                VLOG(2) << "    value[" << dofNoNonBcGlobalRow << "," << dofNoNonBcGlobalColumn << "]=" << matrixGlobalNonBc[dofNoNonBcGlobalRow*nColumnsMatrixNonBc + dofNoNonBcGlobalColumn]
                  << " at matrix[" << dofNoGlobalNaturalRow << "," << dofNoGlobalNaturalColumn << "]";

                matrixGlobalNatural[dofNoGlobalNaturalRow*nColumnsGlobal + dofNoGlobalNaturalColumn] = matrixGlobalNonBc[dofNoNonBcGlobalRow*nColumnsMatrixNonBc + dofNoNonBcGlobalColumn];
              }
            }
          }
        }
      }
    }

    if (nRowsGlobal < 200 && VLOG_IS_ON(1))
    {
      std::stringstream s2;
      for (int j = 0; j < nRowsGlobal; j++)
      {
        for (int i = 0; i < nColumnsGlobal; i++)
        {
          if (fabs(matrixGlobalNatural[j*nColumnsGlobal + i]) < 1e-7)
            s2 << ".";
          else if (fabs(matrixGlobalNatural[j*nColumnsGlobal + i] - 1.0) < 1e-7)
            s2 << "1";
          else if (fabs(matrixGlobalNatural[j*nColumnsGlobal + i] - matrixGlobalNatural[i*nColumnsGlobal + j]) < 2e-4)
            s2 << "s";
          else
            s2 << "x";
        }
        s2 << std::endl;
      }

      VLOG(1) << "final matrix " << nRowsGlobal << "x" << nColumnsGlobal << " entries, " << nDisplacementDofsGlobal << " dofs_u, " << nPressureDofsGlobal << " dofs_p, nonzeros: " << std::endl << s2.str();
    }

    // write file
    std::ofstream file;
    std::stringstream matrixName;

    if (filename.find("/") != std::string::npos)
    {
      matrixName << filename.substr(filename.rfind("/")+1);
    }
    else
    {
      matrixName << filename;
    }
    //matrixName << "r" << nRanks;

    filename += std::string(".m");
    OutputWriter::Generic::openFile(file, filename);

    // write header
    file << "% " << nRowsGlobal << "x" << nColumnsGlobal << " entries, " << nDisplacementDofsGlobal << " dofs_u, " << nPressureDofsGlobal << " dofs_p, "
      << meshPartition->nRanks() << " MPI ranks" << std::endl
      << matrixName.str() << " = ..." << std::endl << "[";

    for (int j = 0; j < nRowsGlobal; j++)
    {
      if (j != 0)
        file << "; ..." << std::endl << " ";

      for (int i = 0; i < nColumnsGlobal; i++)
      {
        if (i != 0)
          file << ", ";

        file << matrixGlobalNatural[j*nColumnsGlobal + i];
      }
    }
    file << "];" << std::endl;

    file.close();

    LOG(INFO) << "Matrix written to \"" << filename << "\".";
  }
}

//! get a submatrix of the upper left part (only displacements)
template<typename DisplacementsFunctionSpaceType, typename PressureFunctionSpaceType, int nDisplacementComponents>
Mat PartitionedPetscMatForHyperelasticity<DisplacementsFunctionSpaceType,PressureFunctionSpaceType,nDisplacementComponents>::
getSubmatrix(int rowVariableNo, int columnVariableNo)
{
  // assert that rowVariableNo and columnVariableNo have valid values
  if (nDisplacementComponents == 3)
  {
    assert(rowVariableNo >= 0 && rowVariableNo < 2);
    assert(columnVariableNo >= 0 && columnVariableNo < 2);
  }
  else if (nDisplacementComponents == 6)
  {
    assert(rowVariableNo >= 0 && rowVariableNo < 3);
    assert(columnVariableNo >= 0 && columnVariableNo < 3);
  }

  MPI_Comm mpiCommunicator = partitionedPetscVecForHyperelasticity_->meshPartition()->mpiCommunicator();

  IS indexSetRows = nullptr;
  IS indexSetColumns = nullptr;

  // no velocity components
  if (nDisplacementComponents == 3)
  {
    if (rowVariableNo == 0)
    {
      indexSetRows = partitionedPetscVecForHyperelasticity_->displacementDofsGlobal();
    }
    else if (rowVariableNo == 1)
    {
      indexSetRows = partitionedPetscVecForHyperelasticity_->pressureDofsGlobal();
    }

    if (columnVariableNo == 0)
    {
      indexSetColumns = partitionedPetscVecForHyperelasticity_->displacementDofsGlobal();
    }
    else if (columnVariableNo == 1)
    {
      indexSetColumns = partitionedPetscVecForHyperelasticity_->pressureDofsGlobal();
    }
  }
  else    // with velocity components
  {
    if (rowVariableNo == 0)
    {
      indexSetRows = partitionedPetscVecForHyperelasticity_->displacementDofsGlobal();
    }
    else if (rowVariableNo == 1)
    {
      indexSetRows = partitionedPetscVecForHyperelasticity_->velocityDofsGlobal();
    }
    else if (rowVariableNo == 2)
    {
      indexSetRows = partitionedPetscVecForHyperelasticity_->pressureDofsGlobal();
    }

    if (columnVariableNo == 0)
    {
      indexSetColumns = partitionedPetscVecForHyperelasticity_->displacementDofsGlobal();
    }
    else if (columnVariableNo == 1)
    {
      indexSetColumns = partitionedPetscVecForHyperelasticity_->velocityDofsGlobal();
    }
    else if (columnVariableNo == 2)
    {
      indexSetColumns = partitionedPetscVecForHyperelasticity_->pressureDofsGlobal();
    }
  }

  assert(indexSetRows);
  assert(indexSetColumns);

  Mat submatrix;
  PetscErrorCode ierr;
  ierr = MatGetSubMatrix(this->globalMatrix_, indexSetRows, indexSetColumns, MAT_INITIAL_MATRIX, &submatrix); CHKERRABORT(mpiCommunicator,ierr);

  return submatrix;
}
