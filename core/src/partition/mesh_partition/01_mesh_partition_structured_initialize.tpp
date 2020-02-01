#include "partition/mesh_partition/01_mesh_partition.h"

#include <cstdlib>
#include "utility/vector_operators.h"
#include "function_space/00_function_space_base_dim.h"

namespace Partition
{

template<typename MeshType,typename BasisFunctionType>
void MeshPartition<FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>,Mesh::isStructured<MeshType>>::
initializeHasFullNumberOfNodes()
{
  // determine if the local partition is at the x+/y+/z+ border of the global domain
  for (int i = 0; i < MeshType::dim(); i++)
  {
    assert (beginElementGlobal_[i] + nElementsLocal_[i] <= (int)nElementsGlobal_[i]);

    if (beginElementGlobal_[i] + nElementsLocal_[i] == (int)nElementsGlobal_[i])
    {
      hasFullNumberOfNodes_[i] = true;
    }
  }
}

template<typename MeshType,typename BasisFunctionType>
void MeshPartition<FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>,Mesh::isStructured<MeshType>>::
initializeDofNosLocalNaturalOrdering(std::shared_ptr<FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>> functionSpace)
{
  if (dofNosLocalNaturalOrdering_.empty())
  {

    // resize the vector to hold number of localWithGhosts dofs
    dofNosLocalNaturalOrdering_.resize(nDofsLocalWithGhosts());

    const int nDofsPerNode = FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>::nDofsPerNode();
    int index = 0;
    if (MeshType::dim() == 1)
    {
      // loop over local nodes in local natural numbering
      for (node_no_t nodeX = 0; nodeX < nNodesLocalWithGhosts(0); nodeX++)
      {
        std::array<int,MeshType::dim()> coordinates({(int)nodeX});

        // loop over dofs of node
        for (int dofOnNodeIndex = 0; dofOnNodeIndex < nDofsPerNode; dofOnNodeIndex++)
        {
          dof_no_t dofNoLocal = functionSpace->getNodeNo(coordinates)*nDofsPerNode + dofOnNodeIndex;
          dofNosLocalNaturalOrdering_[index++] = dofNoLocal;
        }
      }
    }
    else if (MeshType::dim() == 2)
    {
      // loop over local nodes in local natural numbering
      for (node_no_t nodeY = 0; nodeY < nNodesLocalWithGhosts(1); nodeY++)
      {
        for (node_no_t nodeX = 0; nodeX < nNodesLocalWithGhosts(0); nodeX++)
        {
          std::array<int,MeshType::dim()> coordinates;
          coordinates[0] = (int)nodeX;
          coordinates[1] = (int)nodeY;

          // loop over dofs of node
          for (int dofOnNodeIndex = 0; dofOnNodeIndex < nDofsPerNode; dofOnNodeIndex++)
          {
            dof_no_t dofNoLocal = functionSpace->getNodeNo(coordinates)*nDofsPerNode + dofOnNodeIndex;
            dofNosLocalNaturalOrdering_[index++] = dofNoLocal;
          }
        }
      }
    }
    else if (MeshType::dim() == 3)
    {
      // loop over local nodes in local natural numbering
      for (node_no_t nodeZ = 0; nodeZ < nNodesLocalWithGhosts(2); nodeZ++)
      {
        for (node_no_t nodeY = 0; nodeY < nNodesLocalWithGhosts(1); nodeY++)
        {
          for (node_no_t nodeX = 0; nodeX < nNodesLocalWithGhosts(0); nodeX++)
          {
            std::array<int,MeshType::dim()> coordinates;
            coordinates[0] = (int)nodeX;
            coordinates[1] = (int)nodeY;
            coordinates[2] = (int)nodeZ;

            // loop over dofs of node
            for (int dofOnNodeIndex = 0; dofOnNodeIndex < nDofsPerNode; dofOnNodeIndex++)
            {
              dof_no_t dofNoLocal = functionSpace->getNodeNo(coordinates)*nDofsPerNode + dofOnNodeIndex;
              dofNosLocalNaturalOrdering_[index++] = dofNoLocal;
            }
          }
        }
      }
    }
  }
}

template<typename MeshType,typename BasisFunctionType>
void MeshPartition<FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>,Mesh::isStructured<MeshType>>::
initialize1NodeMesh()
{
  // We initialize for a mesh with 0 elements but 1 node and dof, this is needed e.g. for a single CellML problem.
  if (rankSubset_->size() > 1)
  {
    LOG(FATAL) << "Cannot run a 1-node problem on multiple (" << rankSubset_->size() << ") ranks.";
  }

  LOG(DEBUG) << "initialize mesh partition for mesh with 1 dof";

  dmElements_ = nullptr;
  beginElementGlobal_[0] = 0;
  nElementsLocal_[0] = 0;
  nElementsGlobal_[0] = 0;
  nRanks_[0] = 1;   // 1 rank
  ownRankPartitioningIndex_[0] = 0;
  localSizesOnPartitions_[0].resize(1);
  localSizesOnPartitions_[0][0] = 1;
  hasFullNumberOfNodes_[0] = true;

  onlyNodalDofLocalNos_.resize(1);
  onlyNodalDofLocalNos_[0] = 0;

  dofNosLocalNaturalOrdering_.resize(1);
  dofNosLocalNaturalOrdering_[0] = 0;
  localToGlobalPetscMappingDofs_ = NULL;
  nDofsLocalWithoutGhosts_ = FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>::nDofsPerNode();

  this->dofNosLocal_.resize(1);
  this->dofNosLocal_[0] = 0;

  PetscErrorCode ierr;
  PetscInt index = 0;
  ierr = ISLocalToGlobalMappingCreate(rankSubset_->mpiCommunicator(), 1, 1, &index, PETSC_COPY_VALUES, &localToGlobalPetscMappingDofs_); CHKERRV(ierr);
}

template<typename MeshType,typename BasisFunctionType>
void MeshPartition<FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>,Mesh::isStructured<MeshType>>::
createDmElements()
{
  dmElements_ = std::make_shared<DM>();

  PetscErrorCode ierr;
  const int ghostLayerWidth = 1;
  const int nDofsPerElement = 1;   // a multiplicity parameter to the items for which the partitioning is generated by PETSc.

  // if there is only 1 rank, there is only one subdomain that contains all elements
  if (this->nRanks() == 1)
  {
    for (int dimensionIndex = 0; dimensionIndex < MeshType::dim(); dimensionIndex++)
    {
      beginElementGlobal_[dimensionIndex] = 0;
      nElementsLocal_[dimensionIndex] = nElementsGlobal_[dimensionIndex];
      nRanks_[dimensionIndex] = 1;
      localSizesOnPartitions_[dimensionIndex].resize(1);
      localSizesOnPartitions_[dimensionIndex][0] = nElementsGlobal_[dimensionIndex];
    }
  }
  else
  {


    // create PETSc DMDA object that is a topology interface handling parallel data layout on structured grids
    if (MeshType::dim() == 1)
    {
      // create 1d decomposition
      ierr = DMDACreate1d(mpiCommunicator(), DM_BOUNDARY_NONE, nElementsGlobal_[0], nDofsPerElement, ghostLayerWidth,
                          NULL, dmElements_.get()); CHKERRV(ierr);

      // get global coordinates of local partition
      PetscInt x, m;
      ierr = DMDAGetCorners(*dmElements_, &x, NULL, NULL, &m, NULL, NULL); CHKERRV(ierr);
      beginElementGlobal_[0] = (global_no_t)x;
      nElementsLocal_[0] = (element_no_t)m;

      // get number of ranks in each coordinate direction
      std::array<PetscInt,1> nRanks;
      ierr = DMDAGetInfo(*dmElements_, NULL, NULL, NULL, NULL, &nRanks[0], NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL); CHKERRV(ierr);
      nRanks_[0] = nRanks[0];

      // get local sizes on the ranks
      const PetscInt *lxData;
      ierr = DMDAGetOwnershipRanges(*dmElements_, &lxData, NULL, NULL);

      VLOG(1) << "nRanks_[0] = " << nRanks_[0];
      VLOG(1) << "lxData: " << intptr_t(lxData);
      localSizesOnPartitions_[0].resize(nRanks_[0]);
      for (int i = 0; i < nRanks_[0]; i++)
      {
        PetscInt l = lxData[i];
        VLOG(1) << "set localSizesOnPartitions_[0][" << i<< "]=" << l;
        localSizesOnPartitions_[0][i] = l;
      }
    }
    else if (MeshType::dim() == 2)
    {
      // create 2d decomposition
      ierr = DMDACreate2d(mpiCommunicator(), DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DMDA_STENCIL_BOX,
                          nElementsGlobal_[0], nElementsGlobal_[1], PETSC_DECIDE, PETSC_DECIDE,
                          nDofsPerElement, ghostLayerWidth, NULL, NULL, dmElements_.get()); CHKERRV(ierr);

      // get global coordinates of local partition
      PetscInt x, y, m, n;
      ierr = DMDAGetCorners(*dmElements_, &x, &y, NULL, &m, &n, NULL); CHKERRV(ierr);
      beginElementGlobal_[0] = (global_no_t)x;
      beginElementGlobal_[1] = (global_no_t)y;
      nElementsLocal_[0] = (element_no_t)m;
      nElementsLocal_[1] = (element_no_t)n;

      // get number of ranks in each coordinate direction
      std::array<PetscInt,2> nRanks;
      ierr = DMDAGetInfo(*dmElements_, NULL, NULL, NULL, NULL, &nRanks[0], &nRanks[1], NULL, NULL, NULL, NULL, NULL, NULL, NULL); CHKERRV(ierr);
      nRanks_[0] = nRanks[0];
      nRanks_[1] = nRanks[1];

      // get local sizes on the ranks
      const PetscInt *lxData;
      const PetscInt *lyData;
      ierr = DMDAGetOwnershipRanges(*dmElements_, &lxData, &lyData, NULL);
      localSizesOnPartitions_[0].assign(lxData, lxData + nRanks_[0]);
      localSizesOnPartitions_[1].assign(lyData, lyData + nRanks_[1]);

      std::array<int,2> meshIsPeriodicInDimension({false,false});
      MPI_Comm cartesianCommunicator;
      MPIUtility::handleReturnValue(
        MPI_Cart_create(mpiCommunicator(), 2, nRanks_.data(), meshIsPeriodicInDimension.data(), true, &cartesianCommunicator),
      "MPI_Cart_create");
    }
    else if (MeshType::dim() == 3)
    {
      if (nElementsGlobal_[0] != 1 && nElementsGlobal_[1] != 1 && nElementsGlobal_[2] != 1)
      {

        // create 3d decomposition
        ierr = DMDACreate3d(mpiCommunicator(), DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DMDA_STENCIL_BOX,
                            nElementsGlobal_[0], nElementsGlobal_[1], nElementsGlobal_[2],
                            PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE,
                            nDofsPerElement, ghostLayerWidth, NULL, NULL, NULL, dmElements_.get()); CHKERRV(ierr);

        // get global coordinates of local partition
        PetscInt x, y, z, m, n, p;
        ierr = DMDAGetCorners(*dmElements_, &x, &y, &z, &m, &n, &p); CHKERRV(ierr);
        beginElementGlobal_[0] = (global_no_t)x;
        beginElementGlobal_[1] = (global_no_t)y;
        beginElementGlobal_[2] = (global_no_t)z;
        nElementsLocal_[0] = (element_no_t)m;
        nElementsLocal_[1] = (element_no_t)n;
        nElementsLocal_[2] = (element_no_t)p;

        // get number of ranks in each coordinate direction
        std::array<PetscInt,3> nRanks;
        ierr = DMDAGetInfo(*dmElements_, NULL, NULL, NULL, NULL, &nRanks[0], &nRanks[1], &nRanks[2], NULL, NULL, NULL, NULL, NULL, NULL); CHKERRV(ierr);
        nRanks_[0] = nRanks[0];
        nRanks_[1] = nRanks[1];
        nRanks_[2] = nRanks[2];

        // get local sizes on the ranks
        const PetscInt *lxData;
        const PetscInt *lyData;
        const PetscInt *lzData;
        ierr = DMDAGetOwnershipRanges(*dmElements_, &lxData, &lyData, &lzData);
        localSizesOnPartitions_[0].assign(lxData, lxData + nRanks_[0]);
        localSizesOnPartitions_[1].assign(lyData, lyData + nRanks_[1]);
        localSizesOnPartitions_[2].assign(lzData, lzData + nRanks_[2]);

    #if 0  // a suggestion from a HLRS course, not needed
        // create cartesian communciator using MPI_Cart_Create
        std::array<int,3> meshIsPeriodicInDimension({false,false,false});
        MPI_Comm cartesianCommunicator;
        MPIUtility::handleReturnValue(
          MPI_Cart_create(mpiCommunicator(), 2, nRanks_.data(), meshIsPeriodicInDimension.data(), true, &cartesianCommunicator),
        "MPI_Cart_create");
    #endif

      }
      else
      {
        // if the domain is rather 2D with one dimension being only 1 element thick:
        if (nElementsGlobal_[0] == 1)
        {
          // create 2d decomposition
          ierr = DMDACreate2d(mpiCommunicator(), DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DMDA_STENCIL_BOX,
                              nElementsGlobal_[1], nElementsGlobal_[2], PETSC_DECIDE, PETSC_DECIDE,
                              nDofsPerElement, ghostLayerWidth, NULL, NULL, dmElements_.get()); CHKERRV(ierr);

          // get global coordinates of local partition
          PetscInt x, y, m, n;
          ierr = DMDAGetCorners(*dmElements_, &x, &y, NULL, &m, &n, NULL); CHKERRV(ierr);
          beginElementGlobal_[0] = 0;
          beginElementGlobal_[1] = (global_no_t)x;
          beginElementGlobal_[2] = (global_no_t)y;
          nElementsLocal_[0] = 1;
          nElementsLocal_[1] = (element_no_t)m;
          nElementsLocal_[2] = (element_no_t)n;

          // get number of ranks in each coordinate direction
          std::array<PetscInt,3> nRanks;
          ierr = DMDAGetInfo(*dmElements_, NULL, NULL, NULL, NULL, &nRanks[1], &nRanks[2], NULL, NULL, NULL, NULL, NULL, NULL, NULL); CHKERRV(ierr);
          nRanks_[0] = 1;
          nRanks_[1] = nRanks[1];
          nRanks_[2] = nRanks[2];

          // get local sizes on the ranks
          const PetscInt *lxData;
          const PetscInt *lyData;
          ierr = DMDAGetOwnershipRanges(*dmElements_, &lxData, &lyData, NULL);
          localSizesOnPartitions_[0].resize(1);
          localSizesOnPartitions_[0][0] = 1;
          localSizesOnPartitions_[1].assign(lxData, lxData + nRanks_[1]);
          localSizesOnPartitions_[2].assign(lyData, lyData + nRanks_[2]);
        }
        else if (nElementsGlobal_[1] == 1)
        {
          // create 2d decomposition
          ierr = DMDACreate2d(mpiCommunicator(), DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DMDA_STENCIL_BOX,
                              nElementsGlobal_[0], nElementsGlobal_[2], PETSC_DECIDE, PETSC_DECIDE,
                              nDofsPerElement, ghostLayerWidth, NULL, NULL, dmElements_.get()); CHKERRV(ierr);

          // get global coordinates of local partition
          PetscInt x, y, m, n;
          ierr = DMDAGetCorners(*dmElements_, &x, &y, NULL, &m, &n, NULL); CHKERRV(ierr);
          beginElementGlobal_[0] = (global_no_t)x;
          beginElementGlobal_[1] = 0;
          beginElementGlobal_[2] = (global_no_t)y;
          nElementsLocal_[0] = (element_no_t)m;
          nElementsLocal_[1] = 1;
          nElementsLocal_[2] = (element_no_t)n;

          // get number of ranks in each coordinate direction
          std::array<PetscInt,3> nRanks;
          ierr = DMDAGetInfo(*dmElements_, NULL, NULL, NULL, NULL, &nRanks[0], &nRanks[2], NULL, NULL, NULL, NULL, NULL, NULL, NULL); CHKERRV(ierr);
          nRanks_[1] = 1;
          nRanks_[0] = nRanks[0];
          nRanks_[2] = nRanks[2];

          // get local sizes on the ranks
          const PetscInt *lxData;
          const PetscInt *lyData;
          ierr = DMDAGetOwnershipRanges(*dmElements_, &lxData, &lyData, NULL);
          localSizesOnPartitions_[0].assign(lxData, lxData + nRanks_[0]);
          localSizesOnPartitions_[1].resize(1);
          localSizesOnPartitions_[1][0] = 1;
          localSizesOnPartitions_[2].assign(lyData, lyData + nRanks_[2]);
        }
        else if (nElementsGlobal_[2] == 1)
        {
          // create 2d decomposition
          ierr = DMDACreate2d(mpiCommunicator(), DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DMDA_STENCIL_BOX,
                              nElementsGlobal_[0], nElementsGlobal_[1], PETSC_DECIDE, PETSC_DECIDE,
                              nDofsPerElement, ghostLayerWidth, NULL, NULL, dmElements_.get()); CHKERRV(ierr);

          // get global coordinates of local partition
          PetscInt x, y, m, n;
          ierr = DMDAGetCorners(*dmElements_, &x, &y, NULL, &m, &n, NULL); CHKERRV(ierr);
          beginElementGlobal_[0] = (global_no_t)x;
          beginElementGlobal_[1] = (global_no_t)y;
          beginElementGlobal_[2] = 0;
          nElementsLocal_[0] = (element_no_t)m;
          nElementsLocal_[1] = (element_no_t)n;
          nElementsLocal_[2] = 1;

          // get number of ranks in each coordinate direction
          std::array<PetscInt,2> nRanks;
          ierr = DMDAGetInfo(*dmElements_, NULL, NULL, NULL, NULL, &nRanks[0], &nRanks[1], NULL, NULL, NULL, NULL, NULL, NULL, NULL); CHKERRV(ierr);
          nRanks_[2] = 1;
          nRanks_[0] = nRanks[0];
          nRanks_[1] = nRanks[1];

          // get local sizes on the ranks
          const PetscInt *lxData;
          const PetscInt *lyData;
          ierr = DMDAGetOwnershipRanges(*dmElements_, &lxData, &lyData, NULL);
          localSizesOnPartitions_[0].assign(lxData, lxData + nRanks_[0]);
          localSizesOnPartitions_[1].assign(lyData, lyData + nRanks_[1]);
          localSizesOnPartitions_[2].resize(1);
          localSizesOnPartitions_[2][0] = 1;
        }
      }
    }
  }

  initializeHasFullNumberOfNodes();
  setOwnRankPartitioningIndex();

  VLOG(1) << "createDmElements determined the following parameters: "
    << "beginElementGlobal_: " << beginElementGlobal_
    << ", nElementsLocal_: " << nElementsLocal_
    << ", localSizesOnPartitions_: " << localSizesOnPartitions_
    << ", hasFullNumberOfNodes_: " << hasFullNumberOfNodes_
    << ", ownRankPartitioningIndex/nRanks: " << ownRankPartitioningIndex_ << " / " << nRanks_;
}

template<typename MeshType,typename BasisFunctionType>
void MeshPartition<FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>,Mesh::isStructured<MeshType>>::
setNDofsLocalWithoutGhosts()
{
  // compute nDofsLocalWithoutGhosts
  const int nDofsPerNode = FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>::nDofsPerNode();
  this->nDofsLocalWithoutGhosts_ = nNodesLocalWithoutGhosts() * nDofsPerNode;
}

//! determine the values of ownRankPartitioningIndex_
template<typename MeshType,typename BasisFunctionType>
void MeshPartition<FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>,Mesh::isStructured<MeshType>>::
setOwnRankPartitioningIndex()
{
  VLOG(1) << "setOwnRankPartitioningIndex";
  std::array<global_no_t,MeshType::dim()> nodeNoGlobalNatural;

  for (int i = 0; i < MeshType::dim(); i++)
  {
    nodeNoGlobalNatural[i] = beginNodeGlobalNatural(i);
  }

  ownRankPartitioningIndex_ = getPartitioningIndex(nodeNoGlobalNatural);
}

//! fill the dofLocalNo vectors
template<typename MeshType,typename BasisFunctionType>
void MeshPartition<FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>,Mesh::isStructured<MeshType>>::
createLocalDofOrderings()
{
  VLOG(1) << "--------------------";
  VLOG(1) << "createLocalDofOrderings " << MeshType::dim() << "D";

  MeshPartitionBase::createLocalDofOrderings(nDofsLocalWithGhosts());

  // fill onlyNodalDofLocalNos_
  const int nDofsPerNode = FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>::nDofsPerNode();

  // fill the vector of local dof nos, it contains for each node only the first dof (i.e. not derivatives for Hermite)
  int nNodalDofs = nDofsLocalWithGhosts() / nDofsPerNode;
  onlyNodalDofLocalNos_.resize(nNodalDofs);

  dof_no_t localDofNo = 0;
  for (node_no_t localNodeNo = 0; localNodeNo < nNodesLocalWithGhosts(); localNodeNo++)
  {
    onlyNodalDofLocalNos_[localNodeNo] = localDofNo;
    localDofNo += nDofsPerNode;
  }

  // fill ghostDofNosGlobalPetsc_
  int resultIndex = 0;
  int nGhostDofs = nDofsLocalWithGhosts() - nDofsLocalWithoutGhosts();
  ghostDofNosGlobalPetsc_.resize(nGhostDofs);
  if (MeshType::dim() == 1)
  {
    VLOG(1) << "ghost nodes range (global natural) x: [" << beginNodeGlobalNatural(0) + nNodesLocalWithoutGhosts(0) << "," << beginNodeGlobalNatural(0) + nNodesLocalWithGhosts(0) << "]";
    for (global_no_t i = beginNodeGlobalNatural(0) + nNodesLocalWithoutGhosts(0); i < beginNodeGlobalNatural(0) + nNodesLocalWithGhosts(0); i++)
    {
      for (int dofOnNodeIndex = 0; dofOnNodeIndex < nDofsPerNode; dofOnNodeIndex++)
      {
        ghostDofNosGlobalPetsc_[resultIndex++] = i*nDofsPerNode + dofOnNodeIndex;
      }
    }
  }
  else if (MeshType::dim() == 2)
  {
    VLOG(1) << "ghost nodes range (global natural) x: [" << beginNodeGlobalNatural(0) + nNodesLocalWithoutGhosts(0) << "," << beginNodeGlobalNatural(0) + nNodesLocalWithGhosts(0) << "]";
    VLOG(1) << "ghost nodes range (global natural) y: [" << beginNodeGlobalNatural(1) + nNodesLocalWithoutGhosts(1) << "," << beginNodeGlobalNatural(1) + nNodesLocalWithGhosts(1) << "]";
    for (global_no_t j = beginNodeGlobalNatural(1); j < beginNodeGlobalNatural(1) + nNodesLocalWithGhosts(1); j++)
    {
      for (global_no_t i = beginNodeGlobalNatural(0); i < beginNodeGlobalNatural(0) + nNodesLocalWithGhosts(0); i++)
      {
        // if node (i,j) is a ghost node
        bool isAtYPlus = (j >= beginNodeGlobalNatural(1) + nNodesLocalWithoutGhosts(1) && j < beginNodeGlobalNatural(1) + nNodesLocalWithGhosts(1));
        bool isAtXPlus = (i >= beginNodeGlobalNatural(0) + nNodesLocalWithoutGhosts(0) && i < beginNodeGlobalNatural(0) + nNodesLocalWithGhosts(0));

        if (isAtXPlus || isAtYPlus)
        {

          std::array<global_no_t,MeshType::dim()> coordinatesGlobal;
          coordinatesGlobal[0] = i;
          coordinatesGlobal[1] = j;

          global_no_t nodeNoGlobalPetsc = getNodeNoGlobalPetsc(coordinatesGlobal);

          // loop over dofs of this node and assign global dof nos
          for (int dofOnNodeIndex = 0; dofOnNodeIndex < nDofsPerNode; dofOnNodeIndex++)
          {
            VLOG(1) << "j=" << j << " i=" << i << " nodeNoGlobalPetsc=" << nodeNoGlobalPetsc
              << ", resultIndex=" << resultIndex;
            ghostDofNosGlobalPetsc_[resultIndex++] = nodeNoGlobalPetsc*nDofsPerNode + dofOnNodeIndex;
          }
        }
      }
    }
  }
  else if (MeshType::dim() == 3)
  {
    VLOG(1) << "ghost nodes range (global natural) x: [" << beginNodeGlobalNatural(0) + nNodesLocalWithoutGhosts(0) << "," << beginNodeGlobalNatural(0) + nNodesLocalWithGhosts(0) << "]";
    VLOG(1) << "ghost nodes range (global natural) y: [" << beginNodeGlobalNatural(1) + nNodesLocalWithoutGhosts(1) << "," << beginNodeGlobalNatural(1) + nNodesLocalWithGhosts(1) << "]";
    VLOG(1) << "ghost nodes range (global natural) z: [" << beginNodeGlobalNatural(2) + nNodesLocalWithoutGhosts(2) << "," << beginNodeGlobalNatural(2) + nNodesLocalWithGhosts(2) << "]";

    for (global_no_t k = beginNodeGlobalNatural(2); k < beginNodeGlobalNatural(2) + nNodesLocalWithGhosts(2); k++)
    {
      for (global_no_t j = beginNodeGlobalNatural(1); j < beginNodeGlobalNatural(1) + nNodesLocalWithGhosts(1); j++)
      {
        for (global_no_t i = beginNodeGlobalNatural(0); i < beginNodeGlobalNatural(0) + nNodesLocalWithGhosts(0); i++)
        {
          bool isAtZPlus = (k >= beginNodeGlobalNatural(2) + nNodesLocalWithoutGhosts(2) && k < beginNodeGlobalNatural(2) + nNodesLocalWithGhosts(2));
          bool isAtYPlus = (j >= beginNodeGlobalNatural(1) + nNodesLocalWithoutGhosts(1) && j < beginNodeGlobalNatural(1) + nNodesLocalWithGhosts(1));
          bool isAtXPlus = (i >= beginNodeGlobalNatural(0) + nNodesLocalWithoutGhosts(0) && i < beginNodeGlobalNatural(0) + nNodesLocalWithGhosts(0));

          // if node (i,j,k) is a ghost node
          if (isAtZPlus || isAtYPlus || isAtXPlus)
          {
            std::array<global_no_t,MeshType::dim()> coordinatesGlobal;
            coordinatesGlobal[0] = i;
            coordinatesGlobal[1] = j;
            coordinatesGlobal[2] = k;

            global_no_t nodeNoGlobalPetsc = getNodeNoGlobalPetsc(coordinatesGlobal);

            for (int dofOnNodeIndex = 0; dofOnNodeIndex < nDofsPerNode; dofOnNodeIndex++)
            {
              VLOG(1) << "k=" << k << " j=" << j << " i=" << i << " nodeNoGlobalPetsc=" << nodeNoGlobalPetsc
                << ", resultIndex=" << resultIndex;
              ghostDofNosGlobalPetsc_[resultIndex++] = nodeNoGlobalPetsc*nDofsPerNode + dofOnNodeIndex;
            }
          }
        }
      }
    }
  }

  VLOG(1) << "VecCreateGhost, nDofsLocalWithoutGhosts: " << nDofsLocalWithoutGhosts() << ", nDofsGlobal: " << nDofsGlobal()
    << ", nGhostDofs: " << nGhostDofs << ", ghostDofNosGlobalPetsc_: " << ghostDofNosGlobalPetsc_;

  // create localToGlobalPetscMappingDofs_
  PetscErrorCode ierr;
  Vec temporaryVector;
  ierr = VecCreateGhost(this->mpiCommunicator(), nDofsLocalWithoutGhosts(),
                        nDofsGlobal(), nGhostDofs, ghostDofNosGlobalPetsc_.data(), &temporaryVector); CHKERRV(ierr);

  // retrieve local to global mapping
  ierr = VecGetLocalToGlobalMapping(temporaryVector, &localToGlobalPetscMappingDofs_); CHKERRV(ierr);
  //ierr = VecDestroy(&temporaryVector); CHKERRV(ierr);

  // VecCreateGhost(MPI_Comm comm,PetscInt n,PetscInt N,PetscInt nghost,const PetscInt ghosts[],Vec *vv)

  VLOG(1) << "n=" << nDofsLocalWithoutGhosts() << ", N=" << nDofsGlobal() << ", nghost=" << nGhostDofs << " ghosts:" << ghostDofNosGlobalPetsc_;
  VLOG(1) << "Result: " << localToGlobalPetscMappingDofs_;
}




}  // namespace
