#include "partition/mesh_partition/01_mesh_partition_composite.h"

namespace Partition
{

//! number of elements in the current partition
template<int D, typename BasisFunctionType>
element_no_t MeshPartition<FunctionSpace::FunctionSpace<Mesh::CompositeOfDimension<D>,BasisFunctionType>,Mesh::CompositeOfDimension<D>>::
nElementsLocal() const
{
  return nElementsLocal_;
}

//! number of elements in total
template<int D, typename BasisFunctionType>
global_no_t MeshPartition<FunctionSpace::FunctionSpace<Mesh::CompositeOfDimension<D>,BasisFunctionType>,Mesh::CompositeOfDimension<D>>::
nElementsGlobal() const
{
  return nElementsGlobal_;
}

//! number of dofs in the local partition
template<int D, typename BasisFunctionType>
dof_no_t MeshPartition<FunctionSpace::FunctionSpace<Mesh::CompositeOfDimension<D>,BasisFunctionType>,Mesh::CompositeOfDimension<D>>::
nDofsLocalWithGhosts() const
{
  const int nDofsPerNode = FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<D>,BasisFunctionType>::nDofsPerNode();
  return nNodesLocalWithGhosts_ * nDofsPerNode;
}

//! number of dofs in the local partition, without ghosts
template<int D, typename BasisFunctionType>
dof_no_t MeshPartition<FunctionSpace::FunctionSpace<Mesh::CompositeOfDimension<D>,BasisFunctionType>,Mesh::CompositeOfDimension<D>>::
nDofsLocalWithoutGhosts() const
{
  const int nDofsPerNode = FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<D>,BasisFunctionType>::nDofsPerNode();
  return nNodesLocalWithoutGhosts_ * nDofsPerNode;
}

//! number of dofs in total
template<int D, typename BasisFunctionType>
global_no_t MeshPartition<FunctionSpace::FunctionSpace<Mesh::CompositeOfDimension<D>,BasisFunctionType>,Mesh::CompositeOfDimension<D>>::
nDofsGlobal() const
{
  const int nDofsPerNode = FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<D>,BasisFunctionType>::nDofsPerNode();
  return nNodesGlobal_ * nDofsPerNode;
}

//! number of dofs when summing up all global dofs of the sub meshes, this counts the shared dofs multiple times, for each mesh
template<int D, typename BasisFunctionType>
global_no_t MeshPartition<FunctionSpace::FunctionSpace<Mesh::CompositeOfDimension<D>,BasisFunctionType>,Mesh::CompositeOfDimension<D>>::
nDofsGlobalForBoundaryConditions() const
{
  global_no_t result = 0;

  for (int subMeshIndex = 0; subMeshIndex < nSubMeshes_; subMeshIndex++)
  {
    result += subFunctionSpaces_[subMeshIndex]->nDofsGlobal();
  }
  return result;
}

//! number of nodes in the local partition
template<int D, typename BasisFunctionType>
node_no_t MeshPartition<FunctionSpace::FunctionSpace<Mesh::CompositeOfDimension<D>,BasisFunctionType>,Mesh::CompositeOfDimension<D>>::
nNodesLocalWithGhosts() const
{
  return nNodesLocalWithGhosts_;
}

//! number of nodes in the local partition
template<int D, typename BasisFunctionType>
node_no_t MeshPartition<FunctionSpace::FunctionSpace<Mesh::CompositeOfDimension<D>,BasisFunctionType>,Mesh::CompositeOfDimension<D>>::
nNodesLocalWithoutGhosts() const
{
  return nNodesLocalWithoutGhosts_;
}

//! number of nodes in total
template<int D, typename BasisFunctionType>
global_no_t MeshPartition<FunctionSpace::FunctionSpace<Mesh::CompositeOfDimension<D>,BasisFunctionType>,Mesh::CompositeOfDimension<D>>::
nNodesGlobal() const
{
  return nNodesGlobal_;
}

template<int D, typename BasisFunctionType>
global_no_t MeshPartition<FunctionSpace::FunctionSpace<Mesh::CompositeOfDimension<D>,BasisFunctionType>,Mesh::CompositeOfDimension<D>>::
nNodesGlobal(int coordinateDirection) const
{
  assert (nSubMeshes_ > 0);
  return subFunctionSpaces_[0]->nNodesGlobal(coordinateDirection);
}

//! get the number of nodes in the global Petsc ordering that are in partitions prior to the own rank
template<int D, typename BasisFunctionType>
global_no_t MeshPartition<FunctionSpace::FunctionSpace<Mesh::CompositeOfDimension<D>,BasisFunctionType>,Mesh::CompositeOfDimension<D>>::
beginNodeGlobalPetsc() const
{
  return nonDuplicateNodeNoGlobalBegin_;
}

//! returns the number of submeshes
template<int D, typename BasisFunctionType>
int MeshPartition<FunctionSpace::FunctionSpace<Mesh::CompositeOfDimension<D>,BasisFunctionType>,Mesh::CompositeOfDimension<D>>::
nSubMeshes() const
{
  return nSubMeshes_;
}

//! get the local to global mapping for the current partition, for the dof numbering
template<int D, typename BasisFunctionType>
ISLocalToGlobalMapping MeshPartition<FunctionSpace::FunctionSpace<Mesh::CompositeOfDimension<D>,BasisFunctionType>,Mesh::CompositeOfDimension<D>>::
localToGlobalMappingDofs()
{
  return localToGlobalPetscMappingDofs_;
}

//! get the global natural element no for a local element no
template<int D, typename BasisFunctionType>
global_no_t MeshPartition<FunctionSpace::FunctionSpace<Mesh::CompositeOfDimension<D>,BasisFunctionType>,Mesh::CompositeOfDimension<D>>::
getElementNoGlobalNatural(element_no_t elementNoLocal) const
{
  // from a local element no in the composite numbering get the subMeshNo and the no in the submesh-based numbering
  int subMeshNo = 0;
  element_no_t elementOnMeshNoLocal = 0;
  getSubMeshNoAndElementNoLocal(elementNoLocal, subMeshNo, elementOnMeshNoLocal);

  // add the global number on previous sub meshes
  global_no_t elementNoGlobalNatural = 0;
  for (int subMeshIndex = 0; subMeshIndex < subMeshNo; subMeshIndex++)
  {
    elementNoGlobalNatural += subFunctionSpaces_[subMeshIndex]->nElementsGlobal();
  }

  // add the elementNo in the current submesh
  elementNoGlobalNatural += subFunctionSpaces_[subMeshNo]->meshPartition()->getElementNoGlobalNatural(elementOnMeshNoLocal);

  return elementNoGlobalNatural;
}

//! get the global natural element no for a local element no
template<int D, typename BasisFunctionType>
global_no_t MeshPartition<FunctionSpace::FunctionSpace<Mesh::CompositeOfDimension<D>,BasisFunctionType>,Mesh::CompositeOfDimension<D>>::
getNodeNoGlobalNatural(element_no_t elementNoLocal, int nodeIndex) const
{
  // from a local element no in the composite numbering get the subMeshNo and the no in the submesh-based numbering
  int subMeshNo = 0;
  element_no_t elementOnMeshNoLocal = 0;
  getSubMeshNoAndElementNoLocal(elementNoLocal, subMeshNo, elementOnMeshNoLocal);

  // get node no on subMesh
  node_no_t nodeNoLocal = subFunctionSpaces_[subMeshNo]->getNodeNo(elementOnMeshNoLocal, nodeIndex);

  // get global natural no on subMesh
  std::array<global_no_t,D> coordinatesGlobal = subFunctionSpaces_[subMeshNo]->meshPartition()->getCoordinatesGlobal(nodeNoLocal);
  global_no_t nodeNoGlobalNaturalSubMesh = subFunctionSpaces_[subMeshNo]->meshPartition()->getNodeNoGlobalNatural(coordinatesGlobal);

  // add the global number on previous sub meshes
  global_no_t nodeNoGlobalNatural = 0;
  for (int subMeshIndex = 0; subMeshIndex < subMeshNo; subMeshIndex++)
  {
    nodeNoGlobalNatural += subFunctionSpaces_[subMeshIndex]->nNodesGlobal();
  }

  // add the nodeNoGlobalNatural in the current submesh
  nodeNoGlobalNatural += nodeNoGlobalNaturalSubMesh;

  return nodeNoGlobalNatural;
}

//! get the node no in global petsc ordering from a local node no
template<int D, typename BasisFunctionType>
global_no_t MeshPartition<FunctionSpace::FunctionSpace<Mesh::CompositeOfDimension<D>,BasisFunctionType>,Mesh::CompositeOfDimension<D>>::
getNodeNoGlobalPetsc(node_no_t nodeNoLocal) const
{
  std::pair<int,node_no_t> result = nodeNoNonDuplicateLocalToMeshAndDuplicateLocal_[nodeNoLocal];
  int subMeshNo = result.first;
  nodeNoLocal = result.second;

  return meshAndNodeNoLocalToNodeNoNonDuplicateGlobal_[subMeshNo][nodeNoLocal];
}

//! transfer the local nos in global dof nos, using the PETSc localToGlobal mapping for the dofs
template<int D, typename BasisFunctionType>
void MeshPartition<FunctionSpace::FunctionSpace<Mesh::CompositeOfDimension<D>,BasisFunctionType>,Mesh::CompositeOfDimension<D>>::
getDofNoGlobalPetsc(const std::vector<dof_no_t> &dofNosLocal, std::vector<PetscInt> &dofNosGlobalPetsc) const
{
  dofNosGlobalPetsc.resize(dofNosLocal.size());

  // transfer the local indices to global indices
  PetscErrorCode ierr;
  ierr = ISLocalToGlobalMappingApply(localToGlobalPetscMappingDofs_, dofNosLocal.size(), dofNosLocal.data(), dofNosGlobalPetsc.data()); CHKERRV(ierr);
}

//! get the global petsc dof no for the local no, using the PETSc localToGlobal mapping for the dofs
template<int D, typename BasisFunctionType>
global_no_t MeshPartition<FunctionSpace::FunctionSpace<Mesh::CompositeOfDimension<D>,BasisFunctionType>,Mesh::CompositeOfDimension<D>>::
getDofNoGlobalPetsc(dof_no_t dofNoLocal) const
{
  PetscInt dofNoGlobal;
  PetscErrorCode ierr;
  ierr = ISLocalToGlobalMappingApply(localToGlobalPetscMappingDofs_, 1, &dofNoLocal, &dofNoGlobal); CHKERRQ(ierr);
  return (global_no_t)dofNoGlobal;
}

//! get the local element no. from the global no., set isOnLocalDomain to true if the node with global coordinates is in the local domain
template<int D, typename BasisFunctionType>
element_no_t MeshPartition<FunctionSpace::FunctionSpace<Mesh::CompositeOfDimension<D>,BasisFunctionType>,Mesh::CompositeOfDimension<D>>::
getElementNoLocal(global_no_t elementNoGlobalPetsc, bool &isOnLocalDomain) const
{
  if (elementNoGlobalBegin_ <= elementNoGlobalPetsc && elementNoGlobalPetsc < elementNoGlobalBegin_+nElementsLocal_)
  {
    isOnLocalDomain = true;
    return elementNoGlobalPetsc - elementNoGlobalBegin_;
  }
  else
  {
    isOnLocalDomain = false;
    return -1;
  }
}

//! from the submesh no and the local element no in the submesh numbering get the local element no in the composite numbering
template<int D, typename BasisFunctionType>
element_no_t MeshPartition<FunctionSpace::FunctionSpace<Mesh::CompositeOfDimension<D>,BasisFunctionType>,Mesh::CompositeOfDimension<D>>::
getElementNoLocalFromSubmesh(int subMeshNo, element_no_t elementNoLocalOnSubMesh)
{
  assert (subMeshNo >= 0 && subMeshNo < nSubMeshes_);

  element_no_t elementNoLocal = 0;
  for (int subMeshIndex = 0; subMeshIndex < subMeshNo; subMeshIndex++)
  {
    elementNoLocal += subFunctionSpaces_[subMeshIndex]->nElementsLocal();
  }
  elementNoLocal += elementNoLocalOnSubMesh;
  return elementNoLocal;
}

//! get the local node no for a global petsc node no, does not work for ghost nodes
template<int D, typename BasisFunctionType>
node_no_t MeshPartition<FunctionSpace::FunctionSpace<Mesh::CompositeOfDimension<D>,BasisFunctionType>,Mesh::CompositeOfDimension<D>>::
getNodeNoLocal(global_no_t nodeNoGlobalPetsc, bool &isLocal) const
{
  node_no_t nodeNoLocal = (node_no_t)(nodeNoGlobalPetsc - nonDuplicateNodeNoGlobalBegin_);
  isLocal = (nodeNoLocal >= 0) && (nodeNoLocal < nDofsLocalWithoutGhosts());
  return nodeNoLocal;
}

//! get the local dof no for a global petsc dof no, does not work for ghost nodes
template<int D, typename BasisFunctionType>
dof_no_t MeshPartition<FunctionSpace::FunctionSpace<Mesh::CompositeOfDimension<D>,BasisFunctionType>,Mesh::CompositeOfDimension<D>>::
getDofNoLocal(global_no_t dofNoGlobalPetsc, bool &isLocal) const
{
  const int nDofsPerNode = FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<D>,BasisFunctionType>::nDofsPerNode();
  global_no_t nodeNoGlobalPetsc = dofNoGlobalPetsc / nDofsPerNode;
  int nodalDofIndex = dofNoGlobalPetsc % nDofsPerNode;
  return getNodeNoLocal(nodeNoGlobalPetsc, isLocal) * nDofsPerNode + nodalDofIndex;
}

//! get a vector of global natural dof nos of the locally stored non-ghost dofs, needed for setParameters callback function in cellml adapter
template<int D, typename BasisFunctionType>
void MeshPartition<FunctionSpace::FunctionSpace<Mesh::CompositeOfDimension<D>,BasisFunctionType>,Mesh::CompositeOfDimension<D>>::
getDofNosGlobalNatural(std::vector<global_no_t> &dofNosGlobalNatural) const
{
  dofNosGlobalNatural.resize(this->nDofsLocalWithoutGhosts());
  int nDofsPerNode = FunctionSpaceType::nDofsPerNode();

  // loop over submeshes
  for (int subMeshNo = 0; subMeshNo < nSubMeshes_; subMeshNo++)
  {
    std::vector<global_no_t> dofNosGlobalNaturalSubMesh;
    subFunctionSpaces_[subMeshNo]->meshPartition()->getDofNosGlobalNatural(dofNosGlobalNaturalSubMesh);

    // add global natural dofs of submesh to vector
    // loop over local dofs of submesh

    // loop over local nodes of the submeshes, also including removed nodes
    for (node_no_t nodeNoLocal = 0; nodeNoLocal < subFunctionSpaces_[subMeshNo]->nNodesLocalWithoutGhosts(); nodeNoLocal++)
    {
      // if currently considered node is not removed
      if (meshAndNodeNoLocalToNodeNoNonDuplicateGlobal_[subMeshNo][nodeNoLocal] != -1)
      {
        node_no_t nodeNoNonDuplicateLocal = meshAndNodeNoLocalToNodeNoNonDuplicateLocal_[subMeshNo][nodeNoLocal];

        for (int nodalDofIndex = 0; nodalDofIndex < nDofsPerNode; nodalDofIndex++)
        {
          dof_no_t dofNoNonDuplicateLocal = nodeNoNonDuplicateLocal*nDofsPerNode + nodalDofIndex;
          dofNosGlobalNatural[dofNoNonDuplicateLocal] = dofNosGlobalNaturalSubMesh[nodeNoLocal*nDofsPerNode + nodalDofIndex];
        }
      }
    }
  }
}

//! from a vector of values of global/natural node numbers remove all that are non-local, nComponents consecutive values for each dof are assumed
template<int D, typename BasisFunctionType>
template <typename T>
void MeshPartition<FunctionSpace::FunctionSpace<Mesh::CompositeOfDimension<D>,BasisFunctionType>,Mesh::CompositeOfDimension<D>>::
extractLocalNodesWithoutGhosts(std::vector<T> &values, int nComponents) const
{
  std::vector<T> result(nNodesLocalWithoutGhosts()*nComponents);

  LOG(DEBUG) << "extractLocalNodesWithoutGhosts, initialize result with size " << result.size()
    << " (nNodesLocalWithoutGhosts: " << nNodesLocalWithoutGhosts() << ", nComponents: " << nComponents
      << "), values.size(): " << values.size();

  // loop over submeshes
  for (int subMeshNo = 0; subMeshNo < nSubMeshes_; subMeshNo++)
  {
    // loop over local nodes of the submeshes, also including removed nodes
    for (node_no_t nodeNoLocal = 0; nodeNoLocal < subFunctionSpaces_[subMeshNo]->nNodesLocalWithoutGhosts(); nodeNoLocal++)
    {
      // if currently considered node is not removed
      if (meshAndNodeNoLocalToNodeNoNonDuplicateGlobal_[subMeshNo][nodeNoLocal] != -1)
      {
        // get the global number of the composite numbering that maps to the current node
        global_no_t nodeNoNonDuplicateGlobal = meshAndNodeNoLocalToNodeNoNonDuplicateGlobal_[subMeshNo][nodeNoLocal];

        // copy the value from the input vector at the global no location to the local location in the result vector
        for (int componentNo = 0; componentNo < nComponents; componentNo++)
        {
          assert(nodeNoNonDuplicateGlobal*nComponents + componentNo < values.size());

          result[nodeNoLocal*nComponents + componentNo] = values[nodeNoNonDuplicateGlobal*nComponents + componentNo];
        }
      }
    }
  }

  // assign the result values to the values vector
  values.assign(result.begin(), result.end());
}

//! from a vector of values of global/natural dofs remove all that are non-local
template<int D, typename BasisFunctionType>
template <typename T>
void MeshPartition<FunctionSpace::FunctionSpace<Mesh::CompositeOfDimension<D>,BasisFunctionType>,Mesh::CompositeOfDimension<D>>::
extractLocalDofsWithoutGhosts(std::vector<T> &values) const
{
  const int nDofsPerNode = FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<D>,BasisFunctionType>::nDofsPerNode();

  std::vector<T> result(nNodesLocalWithoutGhosts()*nDofsPerNode);

  LOG(DEBUG) << "extractLocalDofsWithoutGhosts, initialize result with size " << result.size()
    << " (nNodesLocalWithoutGhosts: " << nNodesLocalWithoutGhosts() << ", nDofsPerNode: " << nDofsPerNode
      << "), values.size(): " << values.size();

  // loop over submeshes
  for (int subMeshNo = 0; subMeshNo < nSubMeshes_; subMeshNo++)
  {
    // loop over local nodes of the submeshes, also including removed nodes
    for (node_no_t nodeNoLocal = 0; nodeNoLocal < subFunctionSpaces_[subMeshNo]->nNodesLocalWithoutGhosts(); nodeNoLocal++)
    {
      assert(meshAndNodeNoLocalToNodeNoNonDuplicateGlobal_[subMeshNo].size() > nodeNoLocal);

      // if currently considered node is not removed
      if (meshAndNodeNoLocalToNodeNoNonDuplicateGlobal_[subMeshNo][nodeNoLocal] != -1)
      {
        // get the global number of the composite numbering that maps to the current node
        global_no_t nodeNoNonDuplicateGlobal = meshAndNodeNoLocalToNodeNoNonDuplicateGlobal_[subMeshNo][nodeNoLocal];

        // copy the value from the input vector at the global no location to the local location in the result vector

        for (int nodalDofIndex = 0; nodalDofIndex < nDofsPerNode; nodalDofIndex++)
        {
          dof_no_t dofNoLocal = nodeNoLocal*nDofsPerNode + nodalDofIndex;
          global_no_t dofNoNonDuplicateGlobal = nodeNoNonDuplicateGlobal*nDofsPerNode + nodalDofIndex;
          assert(dofNoNonDuplicateGlobal < values.size());

          result[dofNoLocal] = values[dofNoNonDuplicateGlobal];
        }
      }
    }
  }

  // assign the result values to the values vector
  values.assign(result.begin(), result.end());
}

template<int D, typename BasisFunctionType>
void MeshPartition<FunctionSpace::FunctionSpace<Mesh::CompositeOfDimension<D>,BasisFunctionType>,Mesh::CompositeOfDimension<D>>::
extractLocalDofsWithoutGhosts(std::vector<double> &vector) const
{
  this->template extractLocalDofsWithoutGhosts<double>(vector);
}

//! output to stream for debugging
template<int D, typename BasisFunctionType>
void MeshPartition<FunctionSpace::FunctionSpace<Mesh::CompositeOfDimension<D>,BasisFunctionType>,Mesh::CompositeOfDimension<D>>::
output(std::ostream &stream)
{
  stream << "CompositeMesh ";
  for (int subMeshNo = 0; subMeshNo < nSubMeshes_; subMeshNo++)
  {
    stream << "subMesh " << subMeshNo << "/" << nSubMeshes_ << ": ";
    subFunctionSpaces_[subMeshNo]->meshPartition()->output(stream);
    stream << "  ";
  }
  stream << getString();
}

//! get a vector of local dof nos, range [0,nDofsLocalWithoutGhosts] are the dofs without ghost dofs, the whole vector are the dofs with ghost dofs
//! @param onlyNodalValues: if for Hermite only get every second dof such that derivatives are not returned
template<int D, typename BasisFunctionType>
const std::vector<PetscInt> &MeshPartition<FunctionSpace::FunctionSpace<Mesh::CompositeOfDimension<D>,BasisFunctionType>,Mesh::CompositeOfDimension<D>>::
dofNosLocal(bool onlyNodalValues) const
{
  if (onlyNodalValues)
  {
    return onlyNodalDofLocalNos_;
  }
  else
  {
    return this->dofNosLocal_;
  }
}

// use getDofNoGlobalPetsc(dofNosLocal(), ...) to get dofNosGlobalPetsc

//! get the global dof nos of the ghost dofs in the local partition
template<int D, typename BasisFunctionType>
const std::vector<PetscInt> &MeshPartition<FunctionSpace::FunctionSpace<Mesh::CompositeOfDimension<D>,BasisFunctionType>,Mesh::CompositeOfDimension<D>>::
ghostDofNosGlobalPetsc() const
{
  return ghostDofNosGlobalPetsc_;
}

//! check if the given dof is owned by the own rank, then return true, if not, neighbourRankNo is set to the rank by which the dof is owned
template<int D, typename BasisFunctionType>
bool MeshPartition<FunctionSpace::FunctionSpace<Mesh::CompositeOfDimension<D>,BasisFunctionType>,Mesh::CompositeOfDimension<D>>::
isNonGhost(node_no_t nodeNoLocal, int &neighbourRankNo) const
{
  std::vector<std::pair<int,node_no_t>> subMeshesWithNodes;
  getSubMeshesWithNodes(nodeNoLocal, subMeshesWithNodes);

  VLOG(2) << "isNonGhost, nodeNoLocal: " << nodeNoLocal << ", subMeshesWithNodes (submesh no, node no in that submesh): " << subMeshesWithNodes;

  for (const std::pair<int,node_no_t> &subMeshWithNodes : subMeshesWithNodes)
  {
    int subMeshNo = subMeshWithNodes.first;
    node_no_t nodeNoLocalOnMesh = subMeshWithNodes.second;
    if (subFunctionSpaces_[subMeshNo]->meshPartition()->isNonGhost(nodeNoLocalOnMesh, neighbourRankNo))
    {
      VLOG(2) << "yes, isNonGhost on submesh " << subMeshNo << " (local node no " << nodeNoLocalOnMesh << ") on rank " << neighbourRankNo;

      return true;
    }
  }

  //return false;

  // node was no non-ghost on any submesh, this means it is a ghost node
  // the rank that owns it is the one from the submesh where the node is located
  std::pair<int,node_no_t> result = nodeNoNonDuplicateLocalToMeshAndDuplicateLocal_[nodeNoLocal];
  int subMeshNo = result.first;
  node_no_t nodeNoLocalOnMesh = result.second;

  subFunctionSpaces_[subMeshNo]->meshPartition()->isNonGhost(nodeNoLocalOnMesh, neighbourRankNo);

  VLOG(2) << "isNonGhost, nodeNoLocal: " << nodeNoLocal << ", subMeshesWithNodes: " << subMeshesWithNodes;
  VLOG(2) << "  -> is ghost, subMeshNo: " << subMeshNo << ", nodeNoLocalOnMesh: " << nodeNoLocalOnMesh << ", neighbourRankNo: " << neighbourRankNo;

  return false;
}

//! get information about neighbouring rank and boundary elements for specified face,
//! @param neighbourRankNo: the rank of the neighbouring process that shares the face, @param nElements: Size of one-layer mesh that contains boundary elements that touch the neighbouring process
template<int D, typename BasisFunctionType>
void MeshPartition<FunctionSpace::FunctionSpace<Mesh::CompositeOfDimension<D>,BasisFunctionType>,Mesh::CompositeOfDimension<D>>::
getBoundaryElements(Mesh::face_t face, int &neighbourRankNo, std::array<element_no_t,D> &nBoundaryElements, std::vector<dof_no_t> &dofNos)
{
  // This method is only needed for parallel fiber estimation
  LOG(FATAL) << "\"getBoundaryElements\" is not implemented for composite mesh.";
}

//! get the rank no of the neighbour in direction face, -1 if there is no such neighbour
template<int D, typename BasisFunctionType>
int MeshPartition<FunctionSpace::FunctionSpace<Mesh::CompositeOfDimension<D>,BasisFunctionType>,Mesh::CompositeOfDimension<D>>::
neighbourRank(Mesh::face_t face)
{
  // This method is only needed for parallel fiber estimation
  LOG(FATAL) << "\"getBoundaryElements\" is not implemented for composite mesh.";
  return 0;
}

//! from a local node no in the composite numbering get the subMeshNo and the no in the submesh-based numbering
template<int D, typename BasisFunctionType>
void MeshPartition<FunctionSpace::FunctionSpace<Mesh::CompositeOfDimension<D>,BasisFunctionType>,Mesh::CompositeOfDimension<D>>::
getSubMeshNoAndNodeNoLocal(node_no_t nodeNoLocal, int &subMeshNo, node_no_t &nodeOnMeshNoLocal) const
{
  // loop over submeshes
  node_no_t nNodesPreviousSubMeshes = 0;
  for (int subMeshIndex = 0; subMeshIndex < nSubMeshes_; subMeshIndex++)
  {
    node_no_t nNodesCurrentSubMesh = nNonDuplicateNodesWithoutGhosts_[subMeshIndex];
    VLOG(2) << "submesh " << subMeshIndex << " has " << nNodesCurrentSubMesh << " nodes, nNodesPreviousSubMeshes: " << nNodesPreviousSubMeshes;

    if (nodeNoLocal < nNodesPreviousSubMeshes + nNodesCurrentSubMesh)
    {
      VLOG(2) << "local node " << nodeNoLocal << " is on that submesh";

      subMeshNo = subMeshIndex;
      node_no_t nodeNoCompositeOnMeshNoLocal = nodeNoLocal - nNodesPreviousSubMeshes;

      nodeOnMeshNoLocal = nodeNoCompositeOnMeshNoLocal;

      for (node_no_t nodeNo = 0; nodeNo <= nodeOnMeshNoLocal; nodeNo++)
      {
        // if node is removed
        if (removedSharedNodes_[subMeshIndex].find(nodeNo) != removedSharedNodes_[subMeshIndex].end())
        {
          nodeOnMeshNoLocal++;
        }
      }

      return;
    }

    nNodesPreviousSubMeshes += nNodesCurrentSubMesh;
  }

  VLOG(2) << "node local " << nodeNoLocal << " is a ghost node";

  // if we are here, the requested node is a ghost node
  assert (nNodesPreviousSubMeshes == nNodesLocalWithoutGhosts_);

  for (int subMeshIndex = 0; subMeshIndex < nSubMeshes_; subMeshIndex++)
  {
    node_no_t nGhostNodesCurrentSubMesh = nNonDuplicateGhostNodes_[subMeshIndex];
    VLOG(2) << "submesh " << subMeshIndex << " has " << nGhostNodesCurrentSubMesh << " ghost nodes, nNodesPreviousSubMeshes: " << nNodesPreviousSubMeshes;

    if (nodeNoLocal < nNodesPreviousSubMeshes + nGhostNodesCurrentSubMesh)
    {
      VLOG(2) << "local node " << nodeNoLocal << " is on that submesh";

      subMeshNo = subMeshIndex;

      //nodeOnMeshNoLocal = nodeNoLocal - nNodesPreviousSubMeshes + nRemovedNodesNonGhost_[subMeshNo];
      nodeOnMeshNoLocal = nodeNoLocal - nNodesPreviousSubMeshes + nNonDuplicateNodesWithoutGhosts_[subMeshNo] + nRemovedNodesNonGhost_[subMeshNo];
      
      //LOG(DEBUG) << "nodeOnMeshNoLocal = " << nodeNoLocal << " - " << nNodesPreviousSubMeshes << " + " << nRemovedNodesNonGhost_[subMeshNo] << " = " << nodeOnMeshNoLocal;
      VLOG(2) << "nodeOnMeshNoLocal = " << nodeNoLocal << " - " << nNodesPreviousSubMeshes << " + " << nNonDuplicateNodesWithoutGhosts_[subMeshNo]
         << " + " << nRemovedNodesNonGhost_[subMeshNo] << " = " << nodeOnMeshNoLocal;
      VLOG(2) << "nNonDuplicateNodesWithoutGhosts_: " << nNonDuplicateNodesWithoutGhosts_ << ", nRemovedNodesNonGhost_: " << nRemovedNodesNonGhost_;
      break;
    }

    nNodesPreviousSubMeshes += nGhostNodesCurrentSubMesh;
  }
}

//! from a local element no in the composite numbering get the subMeshNo and the no in the submesh-based numbering
template<int D, typename BasisFunctionType>
void MeshPartition<FunctionSpace::FunctionSpace<Mesh::CompositeOfDimension<D>,BasisFunctionType>,Mesh::CompositeOfDimension<D>>::
getSubMeshNoAndElementNoLocal(element_no_t elementNoLocal, int &subMeshNo, element_no_t &elementOnMeshNoLocal) const
{
  element_no_t nElementsPreviousSubMeshes = 0;

  // loop over all meshes
  for (int subMeshIndex = 0; subMeshIndex < nSubMeshes_; subMeshIndex++)
  {
    element_no_t nElementsLocal = subFunctionSpaces_[subMeshIndex]->nElementsLocal();

    if (elementNoLocal >= nElementsPreviousSubMeshes && elementNoLocal < nElementsPreviousSubMeshes + nElementsLocal)
    {
      subMeshNo = subMeshIndex;
      elementOnMeshNoLocal = elementNoLocal - nElementsPreviousSubMeshes;
      return;
    }
    nElementsPreviousSubMeshes += nElementsLocal;
  }
}

//! from a local node no in the composite numbering get the subMeshNo and the no in the submesh-based numbering
template<int D, typename BasisFunctionType>
void MeshPartition<FunctionSpace::FunctionSpace<Mesh::CompositeOfDimension<D>,BasisFunctionType>,Mesh::CompositeOfDimension<D>>::
getSubMeshesWithNodes(node_no_t nodeNoLocal, std::vector<std::pair<int,node_no_t>> &subMeshesWithNodes) const
{
  int subMeshNo = 0;
  node_no_t nodeOnMeshNoLocal = 0;
  getSubMeshNoAndNodeNoLocal(nodeNoLocal, subMeshNo, nodeOnMeshNoLocal);

  VLOG(1) << "node local " << nodeNoLocal << " is node local " << nodeOnMeshNoLocal << " on sub mesh " << subMeshNo;

  subMeshesWithNodes.push_back(std::make_pair(subMeshNo, nodeOnMeshNoLocal));

  // loop over all meshes and shared nodes
  for (int subMeshIndex = 0; subMeshIndex < nSubMeshes_; subMeshIndex++)
  {
    for (std::map<node_no_t,std::pair<int,node_no_t>>::const_iterator iter = removedSharedNodes_[subMeshIndex].begin();
         iter != removedSharedNodes_[subMeshIndex].end(); iter++)
    {
      int meshNoOtherMesh = iter->second.first;
      node_no_t nodeNoOnOtherMesh = iter->second.second;

      if (meshNoOtherMesh == subMeshNo && nodeNoOnOtherMesh == nodeOnMeshNoLocal)
      {
        subMeshesWithNodes.push_back(std::make_pair(subMeshIndex, iter->first));
      }
    }
  }
}

//! from the submesh no and the local node no in the submesh numbering get the local node no in the composite numbering
template<int D, typename BasisFunctionType>
node_no_t MeshPartition<FunctionSpace::FunctionSpace<Mesh::CompositeOfDimension<D>,BasisFunctionType>,Mesh::CompositeOfDimension<D>>::
getNodeNoLocalFromSubmesh(int subMeshNo, node_no_t nodeNoDuplicateOnSubmesh, bool &nodeIsSharedAndRemovedInCurrentMesh) const
{
  assert(subMeshNo >= 0 && subMeshNo < meshAndNodeNoLocalToNodeNoNonDuplicateLocal_.size());
  assert(nodeNoDuplicateOnSubmesh >= 0 && nodeNoDuplicateOnSubmesh < meshAndNodeNoLocalToNodeNoNonDuplicateLocal_[subMeshNo].size());

  nodeIsSharedAndRemovedInCurrentMesh = isDuplicate_[subMeshNo][nodeNoDuplicateOnSubmesh];

  // std::vector<std::map<node_no_t,std::pair<int,node_no_t>>> removedSharedNodes_;   //< removedSharedNodes_[meshNo][nodeNo] = <sameAsInMeshNo,nodeNoOfThatMesh> nodes that are shared between function spaces, they appear only once in the second function space and are removed there (not included in the composite mapping)

  // if the node on the mesh subMeshNo is a duplicate node and therefore its entry in meshAndNodeNoLocalToNodeNoNonDuplicateLocal_ is -1
  if (nodeIsSharedAndRemovedInCurrentMesh)
  {
    // get the same shared node in the other mesh that has its entry in meshAndNodeNoLocalToNodeNoNonDuplicateLocal_
    std::pair<int,node_no_t> sharedNode = removedSharedNodes_[subMeshNo].at(nodeNoDuplicateOnSubmesh);

    // use this entry
    return meshAndNodeNoLocalToNodeNoNonDuplicateLocal_[sharedNode.first][sharedNode.second];
  }

  // node is a normal, non-shared node
  return meshAndNodeNoLocalToNodeNoNonDuplicateLocal_[subMeshNo][nodeNoDuplicateOnSubmesh];
}

}  // namespace

