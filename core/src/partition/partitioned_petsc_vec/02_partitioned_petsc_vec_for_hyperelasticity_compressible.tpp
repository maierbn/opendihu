#include "partition/partitioned_petsc_vec/02_partitioned_petsc_vec_for_hyperelasticity.h"

#include "utility/mpi_utility.h"

// ---- compressible case ----
//! constructor
template<typename DisplacementsFunctionSpaceType, typename PressureFunctionSpaceType, typename Term, int nComponents>
PartitionedPetscVecForHyperelasticity<DisplacementsFunctionSpaceType,PressureFunctionSpaceType,Term,nComponents,std::enable_if_t<!Term::isIncompressible,Term>>::
PartitionedPetscVecForHyperelasticity(
  std::shared_ptr<Partition::MeshPartition<DisplacementsFunctionSpaceType>> meshPartitionDisplacements,
  std::shared_ptr<Partition::MeshPartition<PressureFunctionSpaceType>> meshPartitionPressure,
  std::shared_ptr<SpatialDiscretization::DirichletBoundaryConditions<DisplacementsFunctionSpaceType,nComponents>> dirichletBoundaryConditions, std::string name) :
  PartitionedPetscVecWithDirichletBc<DisplacementsFunctionSpaceType,nComponents,nComponents>(meshPartitionDisplacements, dirichletBoundaryConditions, name, true),
  meshPartitionPressure_(meshPartitionPressure)
{
  // initialize variables for 3 or 6 displacement/velocity components
  LOG(DEBUG) << "\"" << this->name_ << "\" initialize PartitionedPetscVecForHyperelasticity";
  this->initialize(0);

  // create the Petsc Vec
  this->createVector();
}

template<typename DisplacementsFunctionSpaceType, typename PressureFunctionSpaceType, typename Term, int nComponents>
std::shared_ptr<Partition::MeshPartition<PressureFunctionSpaceType>> PartitionedPetscVecForHyperelasticity<DisplacementsFunctionSpaceType,PressureFunctionSpaceType,Term,nComponents,std::enable_if_t<!Term::isIncompressible,Term>>::
meshPartitionPressure()
{
  return meshPartitionPressure_;
}

template<typename DisplacementsFunctionSpaceType, typename PressureFunctionSpaceType, typename Term, int nComponents>
void PartitionedPetscVecForHyperelasticity<DisplacementsFunctionSpaceType,PressureFunctionSpaceType,Term,nComponents,std::enable_if_t<!Term::isIncompressible,Term>>::
dumpGlobalNatural(std::string filename)
{
  // write file
  std::ofstream file;
  std::string vectorName = filename;
  if (vectorName.find("/") != std::string::npos)
    vectorName = vectorName.substr(vectorName.rfind("/")+1);
  filename += std::string(".m");
  OutputWriter::Generic::openFile(file, filename);

  // write header
  file << "% " << this->name_ << ", " << this->nEntriesGlobal_ << " entries, " << this->meshPartition_->nRanks() << " MPI ranks" << std::endl;
  file << getString(false, vectorName);
  LOG(INFO) << "Vector \"" << filename << "\" written.";
  file.close();
}

template<typename DisplacementsFunctionSpaceType, typename PressureFunctionSpaceType, typename Term, int nComponents>
std::string PartitionedPetscVecForHyperelasticity<DisplacementsFunctionSpaceType,PressureFunctionSpaceType,Term,nComponents,std::enable_if_t<!Term::isIncompressible,Term>>::
getString(bool horizontal, std::string vectorName) const
{
  /*does not compile for composite meshes*/
#if 0
  // do not assemble a horizontal string for console in release mode, because this is only needed for debugging output
#ifdef NDEBUG
  if (horizontal)
    return std::string("");
#endif

  std::stringstream result;
  int ownRankNo = this->meshPartition_->ownRankNo();

  // the master rank collects all data and writes the file
  if (ownRankNo == 0)
  {
    // handle displacement values
    std::vector<std::vector<global_no_t>> dofNosGlobalNatural(this->meshPartition_->nRanks());
    std::vector<std::vector<double>> values(this->meshPartition_->nRanks());
    std::vector<int> nDofsOnRank(this->meshPartition_->nRanks());
    std::vector<int> nPressureDofsOnRank(this->meshPartition_->nRanks());

    VLOG(1) << "values: " << values;

    // handle own rank
    // get dof nos in global natural ordering
    int nDofsLocalWithoutGhosts = this->meshPartition_->nDofsLocalWithoutGhosts();
    nDofsOnRank[0] = nDofsLocalWithoutGhosts;

    std::vector<global_no_t> dofNosGlobalNaturalOwn;
    this->meshPartition_->getDofNosGlobalNatural(dofNosGlobalNaturalOwn);

    dofNosGlobalNatural[0].resize(nDofsLocalWithoutGhosts);
    std::copy(dofNosGlobalNaturalOwn.begin(), dofNosGlobalNaturalOwn.end(), dofNosGlobalNatural[0].begin());

    // get displacement values
    values[0].resize(nComponents*nDofsLocalWithoutGhosts);
    for (int i = 0; i < nComponents; i++)
    {
      this->getValues(i, nDofsLocalWithoutGhosts, this->meshPartition_->dofNosLocal().data(), values[0].data() + i*nDofsLocalWithoutGhosts);
    }

    VLOG(1) << "get own displacement values: " << values;

    for (int rankNoK = 0; rankNoK < this->meshPartition_->nRanks(2); rankNoK++)
    {
      for (int rankNoJ = 0; rankNoJ < this->meshPartition_->nRanks(1); rankNoJ++)
      {
        for (int rankNoI = 0; rankNoI < this->meshPartition_->nRanks(0); rankNoI++)
        {
          int rankNo = rankNoK*this->meshPartition_->nRanks(1)*this->meshPartition_->nRanks(0)
            + rankNoJ*this->meshPartition_->nRanks(0) + rankNoI;

          VLOG(1) << "rank " << rankNo;

          // determine number of dofs on rank
          nDofsOnRank[rankNo] = this->meshPartition_->nNodesLocalWithoutGhosts(2, rankNoK)
          * this->meshPartition_->nNodesLocalWithoutGhosts(1, rankNoJ)
          * this->meshPartition_->nNodesLocalWithoutGhosts(0, rankNoI)
          * DisplacementsFunctionSpaceType::nDofsPerNode();

          if (rankNo == 0)
            continue;

          VLOG(1) << ", nDofsOnRank: " << nDofsOnRank[rankNo];

          VLOG(1) << "recv from " << rankNo << " " << nDofsOnRank[rankNo] << " dofs and displacements values";

          // receive dof nos
          dofNosGlobalNatural[rankNo].resize(nDofsOnRank[rankNo]);
          MPIUtility::handleReturnValue(MPI_Recv(dofNosGlobalNatural[rankNo].data(), nDofsOnRank[rankNo], MPI_UNSIGNED_LONG_LONG,
                                                rankNo, 0, this->meshPartition_->rankSubset()->mpiCommunicator(), MPI_STATUS_IGNORE), "MPI_Recv");

          VLOG(1) << "received displacements dofs: " << dofNosGlobalNatural[rankNo];

          // receive values
          values[rankNo].resize(nComponents*nDofsOnRank[rankNo]);
          MPIUtility::handleReturnValue(MPI_Recv(values[rankNo].data(), nComponents*nDofsOnRank[rankNo], MPI_DOUBLE,
                                                rankNo, 0, this->meshPartition_->rankSubset()->mpiCommunicator(), MPI_STATUS_IGNORE), "MPI_Recv");

          VLOG(1) << "received displacements values: " << values[rankNo];
        }
      }
    }

    VLOG(1) << "values: " << values;

    // sort displacement values according to global natural dof no
    std::vector<std::pair<global_no_t,std::array<double,nComponents>>> displacementEntries;
    displacementEntries.reserve(this->nEntriesGlobal_);

    for (int rankNo = 0; rankNo < this->meshPartition_->nRanks(); rankNo++)
    {
      assert(nDofsOnRank[rankNo] == dofNosGlobalNatural[rankNo].size());
      for (int i = 0; i < nDofsOnRank[rankNo]; i++)
      {
        std::pair<global_no_t,std::array<double,nComponents>> displacementEntry;
        displacementEntry.first = dofNosGlobalNatural[rankNo][i];

        std::array<double,nComponents> valuesDof;
        for (int componentNo = 0; componentNo < nComponents; componentNo++)
        {
          valuesDof[componentNo] = values[rankNo][componentNo*nDofsOnRank[rankNo] + i];
        }
        displacementEntry.second = valuesDof;

        displacementEntries.push_back(displacementEntry);
      }
    }

    // sort list according to dof no.s
    std::sort(displacementEntries.begin(), displacementEntries.end(),
              [&](std::pair<global_no_t,std::array<double,nComponents>> a, std::pair<global_no_t,std::array<double,nComponents>> b)
    {
      return a.first < b.first;
    });

    if (VLOG_IS_ON(1))
    {
      VLOG(1) << "dofNosGlobalNatural: " << dofNosGlobalNatural;
      VLOG(1) << "values: " << values;
      VLOG(1) << "nDofsOnRank: " << nDofsOnRank;
      VLOG(1) << "displacementEntries: " << displacementEntries;
    }


    // write file
    std::string newline = "\n";
    std::string separator = ", ";
    std::array<std::string,nComponents> componentNames;

    componentNames[0] = "ux";
    componentNames[1] = "uy";
    componentNames[2] = "uz";
    if (nComponents == 6)
    {
      componentNames[3] = "vx";
      componentNames[4] = "vy";
      componentNames[5] = "vz";
    }

    if (horizontal)
    {
      newline = "";
      separator = ", ";
      result << std::endl;
    }
    else
    {
      result << vectorName << "r" << this->meshPartitionPressure_->nRanks() << " = [";
    }

    // loop over not-pressure components (u and possibly v)
    for (int componentNo = 0; componentNo < nComponents; componentNo++)
    {
      // start of component
      if (horizontal)
      {
        result << componentNames[componentNo] << " = [" << newline;  // print e.g. "ux = ["
      }

      // write displacement values
      for (int i = 0; i < displacementEntries.size(); i++)
      {
        if (i != 0)
          result << separator;
        result << displacementEntries[i].second[componentNo];
      }

      // end of component
      if (horizontal)
      {
        result << newline << "]; " << std::endl;
      }
      else
      {
        result << "]; " << std::endl;
      }
    }
  }
  else
  {
    // send global natural dof nos for displacements
    int nDofsLocalWithoutGhosts = this->meshPartition_->nDofsLocalWithoutGhosts();

    std::vector<global_no_t> dofNosGlobalNatural;
    this->meshPartition_->getDofNosGlobalNatural(dofNosGlobalNatural);

    assert(dofNosGlobalNatural.size() == nDofsLocalWithoutGhosts);

    VLOG(1) << "send to 0 " << nDofsLocalWithoutGhosts << " dofs and displacements values";

    MPIUtility::handleReturnValue(MPI_Send(dofNosGlobalNatural.data(), nDofsLocalWithoutGhosts, MPI_UNSIGNED_LONG_LONG,
                                           0, 0, this->meshPartition_->rankSubset()->mpiCommunicator()), "MPI_Send");

    VLOG(1) << "sent displacements dofs: " << dofNosGlobalNatural;

    // send displacement values
    std::vector<double> values(nComponents*nDofsLocalWithoutGhosts);

    for (int componentNo = 0; componentNo < nComponents; componentNo++)
    {
      this->getValues(componentNo, nDofsLocalWithoutGhosts, this->meshPartition_->dofNosLocal().data(), values.data() + componentNo*nDofsLocalWithoutGhosts);
    }

    MPIUtility::handleReturnValue(MPI_Send(values.data(), nComponents*nDofsLocalWithoutGhosts, MPI_DOUBLE,
                                           0, 0, this->meshPartition_->rankSubset()->mpiCommunicator()), "MPI_Send");

    VLOG(1) << "sent displacements values: " << values;
  }
  return result.str();
#else
  return std::string("");
#endif
}

template<typename DisplacementsFunctionSpaceType, typename PressureFunctionSpaceType, typename Term, int nComponents>
IS PartitionedPetscVecForHyperelasticity<DisplacementsFunctionSpaceType,PressureFunctionSpaceType,Term,nComponents,std::enable_if_t<!Term::isIncompressible,Term>>::
displacementDofsGlobal()
{
  MPI_Comm mpiCommunicator = this->meshPartition_->mpiCommunicator();

  // create index sets of rows of displacement dofs
  int nDisplacementDofsLocal = nDisplacementDofsWithoutBcLocal();
  std::vector<dof_no_t> indices(nDisplacementDofsLocal);
  std::iota(indices.begin(), indices.end(), this->nonBcDofNoGlobalBegin_);

  IS indexSet;
  PetscErrorCode ierr;
  ierr = ISCreateGeneral(mpiCommunicator, nDisplacementDofsLocal, indices.data(), PETSC_COPY_VALUES, &indexSet); CHKERRABORT(mpiCommunicator,ierr);

  return indexSet;
}

template<typename DisplacementsFunctionSpaceType, typename PressureFunctionSpaceType, typename Term, int nComponents>
IS PartitionedPetscVecForHyperelasticity<DisplacementsFunctionSpaceType,PressureFunctionSpaceType,Term,nComponents,std::enable_if_t<!Term::isIncompressible,Term>>::
velocityDofsGlobal()
{
  MPI_Comm mpiCommunicator = this->meshPartition_->mpiCommunicator();

  // create index sets of rows of displacement dofs
  int nDisplacementDofsLocal = nDisplacementDofsWithoutBcLocal();
  int nVelocityDofsLocal = nVelocityDofsWithoutBcLocal();
  std::vector<dof_no_t> indices(nVelocityDofsLocal);
  std::iota(indices.begin(), indices.end(), this->nonBcDofNoGlobalBegin_ + nDisplacementDofsLocal);

  IS indexSet;
  PetscErrorCode ierr;
  ierr = ISCreateGeneral(mpiCommunicator, nVelocityDofsLocal, indices.data(), PETSC_COPY_VALUES, &indexSet); CHKERRABORT(mpiCommunicator,ierr);

  return indexSet;
}

template<typename DisplacementsFunctionSpaceType, typename PressureFunctionSpaceType, typename Term, int nComponents>
IS PartitionedPetscVecForHyperelasticity<DisplacementsFunctionSpaceType,PressureFunctionSpaceType,Term,nComponents,std::enable_if_t<!Term::isIncompressible,Term>>::
pressureDofsGlobal()
{
  MPI_Comm mpiCommunicator = this->meshPartition_->mpiCommunicator();

  // create index sets of rows of pressure dofs
  int nPressureDofsLocal = 0;
  std::vector<dof_no_t> indices(nPressureDofsLocal);
  IS indexSet;
  PetscErrorCode ierr;
  ierr = ISCreateGeneral(mpiCommunicator, nPressureDofsLocal, indices.data(), PETSC_COPY_VALUES, &indexSet); CHKERRABORT(mpiCommunicator,ierr);

  return indexSet;
}

template<typename DisplacementsFunctionSpaceType, typename PressureFunctionSpaceType, typename Term, int nComponents>
dof_no_t PartitionedPetscVecForHyperelasticity<DisplacementsFunctionSpaceType,PressureFunctionSpaceType,Term,nComponents,std::enable_if_t<!Term::isIncompressible,Term>>::
nDisplacementDofsWithoutBcLocal()
{
  dof_no_t result = 0;
  for (int componentNo = 0; componentNo < 3; componentNo++)
  {
    result += this->nNonBcDofsWithoutGhosts_[componentNo];
  }
  return result;
}

template<typename DisplacementsFunctionSpaceType, typename PressureFunctionSpaceType, typename Term, int nComponents>
dof_no_t PartitionedPetscVecForHyperelasticity<DisplacementsFunctionSpaceType,PressureFunctionSpaceType,Term,nComponents,std::enable_if_t<!Term::isIncompressible,Term>>::
nVelocityDofsWithoutBcLocal()
{
  dof_no_t result = 0;
  for (int componentNo = 3; componentNo < nComponents; componentNo++)
  {
    result += this->nNonBcDofsWithoutGhosts_[componentNo];
  }
  return result;
}

template<typename DisplacementsFunctionSpaceType, typename PressureFunctionSpaceType, typename Term, int nComponents>
bool PartitionedPetscVecForHyperelasticity<DisplacementsFunctionSpaceType,PressureFunctionSpaceType,Term,nComponents,std::enable_if_t<!Term::isIncompressible,Term>>::
containsNanOrInf()
{
  // get all local values
  std::vector<double> values;

  for (int componentNo = 0; componentNo < nComponents; componentNo++)
  {
    values.resize(this->meshPartition_->nDofsLocalWithoutGhosts());

    // void getValues(int componentNo, PetscInt ni, const PetscInt ix[], PetscScalar y[]) const;
    this->getValues(componentNo, this->meshPartition_->nDofsLocalWithoutGhosts(), this->meshPartition_->dofNosLocal().data(), values.data());

    // loop over values and check if they are neither nan nor inf
    for (int i = 0; i < values.size(); i++)
    {
      if (!std::isfinite(values[i]) || fabs(values[i]) > 1e+75)
      {
        LOG(ERROR) << "containsNanOrInf(): value " << i << "/" << values.size() << ", component " << componentNo << "/" << nComponents << ": " << values[i];
        return true;
      }
    }
  }
  return false;
}

template<typename DisplacementsFunctionSpaceType, typename PressureFunctionSpaceType, typename Term, int nComponents>
void PartitionedPetscVecForHyperelasticity<DisplacementsFunctionSpaceType,PressureFunctionSpaceType,Term,nComponents,std::enable_if_t<!Term::isIncompressible,Term>>::
interpolateNonDofValues(std::shared_ptr<DisplacementsFunctionSpaceType> displacementsFunctionSpace, std::shared_ptr<PressureFunctionSpaceType> pressureFunctionSpace)
{
  if (!displacementsFunctionSpace->hasTriangleCorners())
    return;

  const bool useRealTriangleShapes = false;   // if real triangles should be used instead of the "deformed" quadratic 6-node-triangles

  this->startGhostManipulation();
  this->zeroGhostBuffer();

  // get all local values
  std::vector<double> valuesLocal;

  // for displacement and velocity components
  for (int componentNo = 0; componentNo < nComponents; componentNo++)
  {
    const int nDofsPerElement = 27;
    ::Mesh::face_or_edge_t edge;

    valuesLocal.resize(this->meshPartition_->nDofsLocalWithoutGhosts());

    // get local values
    // void getValues(int componentNo, PetscInt ni, const PetscInt ix[], PetscScalar y[]) const;
    this->getValues(componentNo, this->meshPartition_->nDofsLocalWithoutGhosts(), this->meshPartition_->dofNosLocal().data(), valuesLocal.data());

    std::vector<double> valuesLocalNew = valuesLocal;

    dof_no_t nDofsLocalWithoutGhosts = displacementsFunctionSpace->nDofsLocalWithoutGhosts();

    // iterate over elements
    for (element_no_t elementNoLocal = 0; elementNoLocal < displacementsFunctionSpace->nElementsLocal(); elementNoLocal++)
    {
      // if the element is a triangle at the corner
      if (this->meshPartition_->elementIsAtCorner(elementNoLocal, edge))
      {
        std::array<dof_no_t,nDofsPerElement> dofNosLocal = displacementsFunctionSpace->getElementDofNosLocal(elementNoLocal);

        if (VLOG_IS_ON(1))
          VLOG(1) << "quadratic element " << elementNoLocal << ", dofNosLocal: " << dofNosLocal;

        int maxZLevel = 3;
        if (dofNosLocal[18] >= nDofsLocalWithoutGhosts)
          maxZLevel = 2;

        // iterate over z levels
        for (int i = 0; i < maxZLevel; i++)
        {
          std::stringstream message;
          if (VLOG_IS_ON(1))
          {
            message << "el " << elementNoLocal << " " << ::Mesh::getString(edge) << " "
              << "(" << valuesLocal[dofNosLocal[i*9 + 1]] << "," << valuesLocal[dofNosLocal[i*9 + 3]]
              << "," << valuesLocal[dofNosLocal[i*9 + 4]] << "," << valuesLocal[dofNosLocal[i*9 + 5]]
              << "," << valuesLocal[dofNosLocal[i*9 + 7]] << ")";
          }

          // depending on the orientation of the triangle
          switch(edge)
          {
          // 0-1-
          case ::Mesh::face_or_edge_t::edge0Minus1Minus:
            if (useRealTriangleShapes)
            {
              valuesLocal[dofNosLocal[i*9 + 0]] = 0.5*(valuesLocal[dofNosLocal[i*9 + 2]] + valuesLocal[dofNosLocal[i*9 + 6]]);
              valuesLocalNew[dofNosLocal[i*9 + 0]] = valuesLocal[dofNosLocal[i*9 + 0]];
            }

            valuesLocalNew[dofNosLocal[i*9 + 1]] = 0.5*(valuesLocal[dofNosLocal[i*9 + 0]] + valuesLocal[dofNosLocal[i*9 + 2]]);
            valuesLocalNew[dofNosLocal[i*9 + 3]] = 0.5*(valuesLocal[dofNosLocal[i*9 + 0]] + valuesLocal[dofNosLocal[i*9 + 6]]);
            valuesLocalNew[dofNosLocal[i*9 + 4]] = 0.5*(valuesLocal[dofNosLocal[i*9 + 0]] + valuesLocal[dofNosLocal[i*9 + 8]]);
            break;

          // 0+1-
          case ::Mesh::face_or_edge_t::edge0Plus1Minus:
            if (useRealTriangleShapes)
            {
              valuesLocal[dofNosLocal[i*9 + 2]] = 0.5*(valuesLocal[dofNosLocal[i*9 + 0]] + valuesLocal[dofNosLocal[i*9 + 8]]);
              valuesLocalNew[dofNosLocal[i*9 + 2]] = valuesLocal[dofNosLocal[i*9 + 2]];
            }
            valuesLocalNew[dofNosLocal[i*9 + 1]] = 0.5*(valuesLocal[dofNosLocal[i*9 + 0]] + valuesLocal[dofNosLocal[i*9 + 2]]);
            valuesLocalNew[dofNosLocal[i*9 + 4]] = 0.5*(valuesLocal[dofNosLocal[i*9 + 2]] + valuesLocal[dofNosLocal[i*9 + 6]]);
            valuesLocalNew[dofNosLocal[i*9 + 5]] = 0.5*(valuesLocal[dofNosLocal[i*9 + 2]] + valuesLocal[dofNosLocal[i*9 + 8]]);
            break;

          // 0-1+
          case ::Mesh::face_or_edge_t::edge0Minus1Plus:
            if (useRealTriangleShapes)
            {
              valuesLocal[dofNosLocal[i*9 + 6]] = 0.5*(valuesLocal[dofNosLocal[i*9 + 0]] + valuesLocal[dofNosLocal[i*9 + 8]]);
              valuesLocalNew[dofNosLocal[i*9 + 6]] = valuesLocal[dofNosLocal[i*9 + 6]];
            }
            valuesLocalNew[dofNosLocal[i*9 + 3]] = 0.5*(valuesLocal[dofNosLocal[i*9 + 0]] + valuesLocal[dofNosLocal[i*9 + 6]]);
            valuesLocalNew[dofNosLocal[i*9 + 4]] = 0.5*(valuesLocal[dofNosLocal[i*9 + 2]] + valuesLocal[dofNosLocal[i*9 + 6]]);
            valuesLocalNew[dofNosLocal[i*9 + 7]] = 0.5*(valuesLocal[dofNosLocal[i*9 + 6]] + valuesLocal[dofNosLocal[i*9 + 8]]);
            break;

          // 0+1+
          case ::Mesh::face_or_edge_t::edge0Plus1Plus:
            if (useRealTriangleShapes)
            {
              valuesLocal[dofNosLocal[i*9 + 8]] = 0.5*(valuesLocal[dofNosLocal[i*9 + 2]] + valuesLocal[dofNosLocal[i*9 + 6]]);
              valuesLocalNew[dofNosLocal[i*9 + 8]] = valuesLocal[dofNosLocal[i*9 + 8]];
            }
            valuesLocalNew[dofNosLocal[i*9 + 4]] = 0.5*(valuesLocal[dofNosLocal[i*9 + 0]] + valuesLocal[dofNosLocal[i*9 + 8]]);
            valuesLocalNew[dofNosLocal[i*9 + 5]] = 0.5*(valuesLocal[dofNosLocal[i*9 + 2]] + valuesLocal[dofNosLocal[i*9 + 8]]);
            valuesLocalNew[dofNosLocal[i*9 + 7]] = 0.5*(valuesLocal[dofNosLocal[i*9 + 6]] + valuesLocal[dofNosLocal[i*9 + 8]]);
            break;

          default:
            break;
          }

          if (VLOG_IS_ON(1))
          {
            message << " -> "
              << "(" << valuesLocalNew[dofNosLocal[i*9 + 1]] << "," << valuesLocalNew[dofNosLocal[i*9 + 3]]
              << "," << valuesLocalNew[dofNosLocal[i*9 + 4]] << "," << valuesLocalNew[dofNosLocal[i*9 + 5]]
              << "," << valuesLocalNew[dofNosLocal[i*9 + 7]] << ")";
            VLOG(1) << message.str();
          }
        }
      }
    }
    this->setValues(componentNo, this->meshPartition_->nDofsLocalWithoutGhosts(), this->meshPartition_->dofNosLocal().data(), valuesLocalNew.data());
  }

  this->finishGhostManipulation();     // communicate and add up values in ghost buffers
}
