#include "partition/mesh_partition/01_mesh_partition_composite.h"

namespace Partition
{

template<int D, typename BasisFunctionType>
MeshPartition<FunctionSpace::FunctionSpace<Mesh::CompositeOfDimension<D>,BasisFunctionType>,Mesh::CompositeOfDimension<D>>::
MeshPartition(const std::vector<std::shared_ptr<FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<D>,BasisFunctionType>>> &subFunctionSpaces) :
  subFunctionSpaces_(subFunctionSpaces)
{
  // initialize number of local and global elements
  this->nElementsLocal_ = 0;
  this->nElementsGlobal_ = 0;

  // iterate over submeshes
  for(const std::shared_ptr<FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<D>,BasisFunctionType>> &subFunctionSpace : subFunctionSpaces_)
  {
    // count number of elements
    this->nElementsLocal_ += subFunctionSpace->nElementsLocal();
    this->nElementsGlobal_ += subFunctionSpace->nElementsGlobal();
  }

  // determine nodes that are the same on multiple meshes

  // get node positions of all meshes
  std::vector<std::vector<Vec3>> nodePositions(subFunctionSpaces_.size());
  std::vector<std::vector<std::pair<Vec3,dof_no_t>>> nodePositionsDofs(subFunctionSpaces_.size());  // the node positions with dofs for every submesh

  // iterate over submeshes and save all node positions
  int i = 0;
  for(const std::shared_ptr<FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<D>,BasisFunctionType>> &subFunctionSpace : subFunctionSpaces_)

  {
    // get geometry field
    subFunctionSpace->geometry()->getValuesWithGhosts(nodePositions[i]);

    for (dof_no_t dofNoLocal = 0; dofNoLocal < subFunctionSpace->nDofsLocalWithGhosts(); dofNoLocal++)
    {
      nodePositionsDofs[i].push_back(std::make_pair(nodePositions[i][dofNoLocal], dofNoLocal));
    }

    // sort according to x coordinate of node positions
    std::sort(nodePositionsDofs[i].begin(), nodePositionsDofs[i].end(), [](const std::pair<Vec3,dof_no_t> &a, const std::pair<Vec3,dof_no_t> &b)
    {
      return a.first[0] < b.first[0];
    });
    i++;
  }


  sharedDofs_.resize(subFunctionSpaces_.size());
  // std::vector<std::map<dof_no_t,std::pair<int,dof_no_t>>>

  // iterate over submeshes
  const double nodePositionEqualTolerance = 1e-3;
  i = 0;
  for(const typename std::vector<std::shared_ptr<FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<D>,BasisFunctionType>>>::const_iterator iter = subFunctionSpaces_.cbegin();
      iter != subFunctionSpaces_.cend(); iter++, i++)
  {
    // find shared dofs in all next meshes

    // loop over node positions of this mesh
    for (dof_no_t dofNoLocal = 0; dofNoLocal < nodePositions[i].size(); dofNoLocal++)
    {
      Vec3 position = nodePositions[i][dofNoLocal];

      // iterate over further submeshes
      int indexOtherMesh = 0;
      for(int indexOtherMesh = i; indexOtherMesh < subFunctionSpaces_.size(); indexOtherMesh++)
      {
        // find node that has closest x coordinate

        // get node with last x coordinate that is lower by more than tolerance
        int k = nodePositionsDofs[indexOtherMesh].size() / 2;
        int kPrevious = -1;
        int lower = 0;
        int upper = nodePositionsDofs[indexOtherMesh].size();

        while (k != kPrevious)
        {
          Vec3 currentNodePosition = nodePositionsDofs[indexOtherMesh][k].first;
          if (currentNodePosition[0] < position[0]-nodePositionEqualTolerance)
          {
            lower = k;
          }
          else
          {
            upper = k;
          }
          k = (upper + lower) / 2;
        }

        // check all node positions
        for (;;)
        {
          Vec3 nodePositionOtherMesh = nodePositionsDofs[indexOtherMesh][k].first;
          dof_no_t dofNoLocalOtherMesh = nodePositionsDofs[indexOtherMesh][k].second;

          if (nodePositionOtherMesh[0] > position[0]+nodePositionEqualTolerance)
          {
            break;
          }

          double distance = MathUtility::distance<3>(position, nodePositionOtherMesh);

          // if the other mesh node is at the same position as the first node
          if (distance <= nodePositionEqualTolerance)
          {
            dof_no_t dofNoOtherMesh = k;
            if (sharedDofs_[indexOtherMesh].find(dofNoOtherMesh) != sharedDofs_[indexOtherMesh].end())
            {
              sharedDofs_[indexOtherMesh][dofNoOtherMesh] = std::make_pair(dofNoLocalOtherMesh, indexOtherMesh);
            }
            break;
          }
        }
      }
    }
  }

  // count number of shared local dofs
  nDofsSharedLocal_ = 0;
  nGhostDofsSharedLocal_ = 0;
  for (int meshIndex = 0; meshIndex < sharedDofs_.size(); meshIndex++)
  {
    std::map<dof_no_t,std::pair<int,dof_no_t>> &sharedDofsInMesh = sharedDofs_[meshIndex];

    for (std::map<dof_no_t,std::pair<int,dof_no_t>>::iterator iter2 = sharedDofsInMesh.begin(); iter2 != sharedDofsInMesh.end(); iter2++)
    {
      if (iter2->second < subFunctionSpaces_[iter2->first].nDofsLocalWithoutGhosts())
      {
        nDofsSharedLocal_++;
      }
      else
      {
        nGhostDofsSharedLocal_++;
      }
    }
  }
}

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
  dof_no_t nDofsLocalWithGhosts = 0;

  // iterate over submeshes
  for(const std::shared_ptr<FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<D>,BasisFunctionType>> &subFunctionSpace : subFunctionSpaces_)
  {
    nDofsLocalWithGhosts += subFunctionSpace->nDofsLocalWithGhosts();
  }

  nDofsLocalWithGhosts -= nDofsSharedLocal_;
  nDofsLocalWithGhosts -= nGhostDofsSharedLocal_;

  return nDofsLocalWithGhosts;
}

//! number of dofs in the local partition, without ghosts
template<int D, typename BasisFunctionType>
dof_no_t MeshPartition<FunctionSpace::FunctionSpace<Mesh::CompositeOfDimension<D>,BasisFunctionType>,Mesh::CompositeOfDimension<D>>::
nDofsLocalWithoutGhosts() const
{
  dof_no_t nDofsLocalWithoutGhosts = 0;

  // iterate over submeshes
  for(const std::shared_ptr<FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<D>,BasisFunctionType>> &subFunctionSpace : subFunctionSpaces_)
  {
    nDofsLocalWithoutGhosts += subFunctionSpace->nDofsLocalWithoutGhosts();
  }
  nDofsLocalWithoutGhosts -= nDofsSharedLocal_;

  return nDofsLocalWithoutGhosts;
}

/*
//! number of dofs in total
template<int D, typename BasisFunctionType>
global_no_t MeshPartition<FunctionSpace::FunctionSpace<Mesh::CompositeOfDimension<D>,BasisFunctionType>,Mesh::CompositeOfDimension<D>>::
nDofsGlobal() const
{
  int nDofsGlobal = 0;

  // iterate over submeshes
  for(const std::shared_ptr<FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<D>,BasisFunctionType>> &subFunctionSpace : subFunctionSpaces_)
  {
    nDofsGlobal += subFunctionSpace->nDofsGlobal();
  }
  return nDofsGlobal;
}

//! number of nodes in the local partition
template<int D, typename BasisFunctionType>
node_no_t MeshPartition<FunctionSpace::FunctionSpace<Mesh::CompositeOfDimension<D>,BasisFunctionType>,Mesh::CompositeOfDimension<D>>::
nNodesLocalWithGhosts() const
{
  int nNodesLocalWithGhosts = 0;

  // iterate over submeshes
  for(const std::shared_ptr<FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<D>,BasisFunctionType>> &subFunctionSpace : subFunctionSpaces_)
  {
    nNodesLocalWithGhosts += subFunctionSpace->nNodesLocalWithGhosts();
  }
  return nNodesLocalWithGhosts;
}

//! number of nodes in the local partition
template<int D, typename BasisFunctionType>
node_no_t MeshPartition<FunctionSpace::FunctionSpace<Mesh::CompositeOfDimension<D>,BasisFunctionType>,Mesh::CompositeOfDimension<D>>::
nNodesLocalWithoutGhosts() const
{
  int nNodesLocalWithoutGhosts = 0;

  // iterate over submeshes
  for(const std::shared_ptr<FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<D>,BasisFunctionType>> &subFunctionSpace : subFunctionSpaces_)
  {
    nNodesLocalWithoutGhosts += subFunctionSpace->nNodesLocalWithoutGhosts();
  }
  return nNodesLocalWithoutGhosts;
}

//! number of nodes in total
template<int D, typename BasisFunctionType>
global_no_t MeshPartition<FunctionSpace::FunctionSpace<Mesh::CompositeOfDimension<D>,BasisFunctionType>,Mesh::CompositeOfDimension<D>>::
nNodesGlobal() const
{
  int nNodesGlobal = 0;

  // iterate over submeshes
  for(const std::shared_ptr<FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<D>,BasisFunctionType>> &subFunctionSpace : subFunctionSpaces_)
  {
    nNodesGlobal += subFunctionSpace->nNodesGlobal();
  }
  return nNodesGlobal;
}

//! get the number of nodes in the global Petsc ordering that are in partitions prior to the own rank
template<int D, typename BasisFunctionType>
global_no_t MeshPartition<FunctionSpace::FunctionSpace<Mesh::CompositeOfDimension<D>,BasisFunctionType>,Mesh::CompositeOfDimension<D>>::
beginNodeGlobalPetsc() const
{
  global_no_t beginNodeGlobalPetsc = 0;

  // iterate over submeshes
  for(const std::shared_ptr<FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<D>,BasisFunctionType>> &subFunctionSpace : subFunctionSpaces_)
  {
    beginNodeGlobalPetsc += subFunctionSpace->beginNodeGlobalPetsc();
  }


  return beginNodeGlobalPetsc;

}

//! number of nodes in total
template<int D, typename BasisFunctionType>
global_no_t MeshPartition<FunctionSpace::FunctionSpace<Mesh::CompositeOfDimension<D>,BasisFunctionType>,Mesh::CompositeOfDimension<D>>::
nNodesGlobal(int coordinateDirection) const;

//! get if there are nodes on both borders in the given coordinate direction
//! this is the case if the partition touches the right/top/back border
//! Consider the partition specified by partitionIndex or the current partition if partitionIndex == -1.
template<int D, typename BasisFunctionType>
bool MeshPartition<FunctionSpace::FunctionSpace<Mesh::CompositeOfDimension<D>,BasisFunctionType>,Mesh::CompositeOfDimension<D>>::
hasFullNumberOfNodes(int coordinateDirection, int partitionIndex = -1) const;

//! get a vector with the local sizes on every rank
template<int D, typename BasisFunctionType>
const std::vector<element_no_t> &MeshPartition<FunctionSpace::FunctionSpace<Mesh::CompositeOfDimension<D>,BasisFunctionType>,Mesh::CompositeOfDimension<D>>::
localSizesOnRanks(int coordinateDirection) const;

//! get the local to global mapping for the current partition, for the dof numbering
template<int D, typename BasisFunctionType>
ISLocalToGlobalMapping MeshPartition<FunctionSpace::FunctionSpace<Mesh::CompositeOfDimension<D>,BasisFunctionType>,Mesh::CompositeOfDimension<D>>::
localToGlobalMappingDofs();

//! get the global natural element no for a local element no
template<int D, typename BasisFunctionType>
global_no_t MeshPartition<FunctionSpace::FunctionSpace<Mesh::CompositeOfDimension<D>,BasisFunctionType>,Mesh::CompositeOfDimension<D>>::
getElementNoGlobalNatural(element_no_t elementNoLocal) const;

//! get the global natural node no for the global coordinates of this node, this can be combined with getCoordinatesGlobal
template<int D, typename BasisFunctionType>
global_no_t MeshPartition<FunctionSpace::FunctionSpace<Mesh::CompositeOfDimension<D>,BasisFunctionType>,Mesh::CompositeOfDimension<D>>::
getNodeNoGlobalNatural(std::array<global_no_t,MeshType::dim()> coordinatesGlobal) const;

//! get the node no in global petsc ordering from a local node no
template<int D, typename BasisFunctionType>
global_no_t MeshPartition<FunctionSpace::FunctionSpace<Mesh::CompositeOfDimension<D>,BasisFunctionType>,Mesh::CompositeOfDimension<D>>::
getNodeNoGlobalPetsc(node_no_t nodeNoLocal) const;

//! transfer the local nos in global dof nos, using the PETSc localToGlobal mapping for the dofs
template<int D, typename BasisFunctionType>
void getDofNoGlobalPetsc(const std::vector<dof_no_t> &dofNosLocal, std::vector<PetscInt> &dofNosGlobalPetsc) const;

//! get the global petsc dof no for the local no, using the PETSc localToGlobal mapping for the dofs
template<int D, typename BasisFunctionType>
global_no_t MeshPartition<FunctionSpace::FunctionSpace<Mesh::CompositeOfDimension<D>,BasisFunctionType>,Mesh::CompositeOfDimension<D>>::
getDofNoGlobalPetsc(dof_no_t dofNoLocal) const;

//! get the global node coordinates (x,y,z) of the node given by its local node no. This also works for ghost nodes.
template<int D, typename BasisFunctionType>
std::array<global_no_t,MeshType::dim()> MeshPartition<FunctionSpace::FunctionSpace<Mesh::CompositeOfDimension<D>,BasisFunctionType>,Mesh::CompositeOfDimension<D>>::
getCoordinatesGlobal(node_no_t nodeNoLocal) const;

//! get the local coordinates for a local node no, also for ghost nodes. With this method and functionSpace->getNodeNo(coordinatesLocal) it is possible to implement a global-to-local mapping.
template<int D, typename BasisFunctionType>
std::array<int,MeshType::dim()> MeshPartition<FunctionSpace::FunctionSpace<Mesh::CompositeOfDimension<D>,BasisFunctionType>,Mesh::CompositeOfDimension<D>>::
getCoordinatesLocal(node_no_t nodeNoLocal) const;

//! from global natural coordinates compute the local coordinates, set isOnLocalDomain to true if the node with global coordinates is in the local domain
template<int D, typename BasisFunctionType>
std::array<int,MeshType::dim()> MeshPartition<FunctionSpace::FunctionSpace<Mesh::CompositeOfDimension<D>,BasisFunctionType>,Mesh::CompositeOfDimension<D>>::
getCoordinatesLocal(std::array<global_no_t,MeshType::dim()> coordinatesGlobal, bool &isOnLocalDomain) const;

//! get the local coordinates for a local element no
template<int D, typename BasisFunctionType>
std::array<int,MeshType::dim()> MeshPartition<FunctionSpace::FunctionSpace<Mesh::CompositeOfDimension<D>,BasisFunctionType>,Mesh::CompositeOfDimension<D>>::
getElementCoordinatesLocal(element_no_t elementNoLocal) const;

//! get the local element no. from coordinates
template<int D, typename BasisFunctionType>
element_no_t MeshPartition<FunctionSpace::FunctionSpace<Mesh::CompositeOfDimension<D>,BasisFunctionType>,Mesh::CompositeOfDimension<D>>::
getElementNoLocal(std::array<int,MeshType::dim()> elementCoordinates) const;

//! get the local element no. from the global no., set isOnLocalDomain to true if the node with global coordinates is in the local domain
template<int D, typename BasisFunctionType>
element_no_t MeshPartition<FunctionSpace::FunctionSpace<Mesh::CompositeOfDimension<D>,BasisFunctionType>,Mesh::CompositeOfDimension<D>>::
getElementNoLocal(global_no_t elementNoGlobalPetsc, bool &isOnLocalDomain) const;

//! get the local node no for a global petsc node no, does not work for ghost nodes
template<int D, typename BasisFunctionType>
node_no_t MeshPartition<FunctionSpace::FunctionSpace<Mesh::CompositeOfDimension<D>,BasisFunctionType>,Mesh::CompositeOfDimension<D>>::
getNodeNoLocal(global_no_t nodeNoGlobalPetsc, bool &isLocal) const;

//! get the local dof no for a global petsc dof no, does not work for ghost nodes
template<int D, typename BasisFunctionType>
dof_no_t MeshPartition<FunctionSpace::FunctionSpace<Mesh::CompositeOfDimension<D>,BasisFunctionType>,Mesh::CompositeOfDimension<D>>::
getDofNoLocal(global_no_t dofNoGlobalPetsc, bool &isLocal) const;

//! from a vector of values of global/natural node numbers remove all that are non-local, nComponents consecutive values for each dof are assumed
template<int D, typename BasisFunctionType>
template <typename T>
void MeshPartition<FunctionSpace::FunctionSpace<Mesh::CompositeOfDimension<D>,BasisFunctionType>,Mesh::CompositeOfDimension<D>>::
extractLocalNodesWithoutGhosts(std::vector<T> &vector, int nComponents=1) const;

//! from a vector of values of global/natural dofs remove all that are non-local
template<int D, typename BasisFunctionType>
template <typename T>
void MeshPartition<FunctionSpace::FunctionSpace<Mesh::CompositeOfDimension<D>,BasisFunctionType>,Mesh::CompositeOfDimension<D>>::
extractLocalDofsWithoutGhosts(std::vector<T> &values) const;

//! from a vector of values of global/natural dofs remove all that are non-local
template<int D, typename BasisFunctionType>
void MeshPartition<FunctionSpace::FunctionSpace<Mesh::CompositeOfDimension<D>,BasisFunctionType>,Mesh::CompositeOfDimension<D>>::
extractLocalDofsWithoutGhosts(std::vector<double> &values) const;

//! get the partition index in a given coordinate direction from the rankNo
template<int D, typename BasisFunctionType>
int MeshPartition<FunctionSpace::FunctionSpace<Mesh::CompositeOfDimension<D>,BasisFunctionType>,Mesh::CompositeOfDimension<D>>::
convertRankNoToPartitionIndex(int coordinateDirection, int rankNo);

//! output to stream for debugging
template<int D, typename BasisFunctionType>
void MeshPartition<FunctionSpace::FunctionSpace<Mesh::CompositeOfDimension<D>,BasisFunctionType>,Mesh::CompositeOfDimension<D>>::
output(std::ostream &stream);

//! get a vector of local dof nos, range [0,nDofsLocalWithoutGhosts] are the dofs without ghost dofs, the whole vector are the dofs with ghost dofs
//! @param onlyNodalValues: if for Hermite only get every second dof such that derivatives are not returned
template<int D, typename BasisFunctionType>
const std::vector<PetscInt> &MeshPartition<FunctionSpace::FunctionSpace<Mesh::CompositeOfDimension<D>,BasisFunctionType>,Mesh::CompositeOfDimension<D>>::
dofNosLocal(bool onlyNodalValues=false) const;

// use getDofNoGlobalPetsc(dofNosLocal(), ...) to get dofNosGlobalPetsc

//! get a vector of global natural dof nos of the locally stored non-ghost dofs, needed for setParameters callback function in cellml adapter
template<int D, typename BasisFunctionType>
void MeshPartition<FunctionSpace::FunctionSpace<Mesh::CompositeOfDimension<D>,BasisFunctionType>,Mesh::CompositeOfDimension<D>>::
getDofNosGlobalNatural(std::vector<global_no_t> &dofNosGlobalNatural) const;

//! get the global dof nos of the ghost dofs in the local partition
template<int D, typename BasisFunctionType>
const std::vector<PetscInt> &MeshPartition<FunctionSpace::FunctionSpace<Mesh::CompositeOfDimension<D>,BasisFunctionType>,Mesh::CompositeOfDimension<D>>::
ghostDofNosGlobalPetsc() const;

//! Initialize the vector dofNosLocalNaturalOrdering_, this needs the functionSpace and has to be called before dofNosLocalNaturalOrdering() can be used.
//! If the vector is already initialized by a previous call to this method, it has no effect.
template<int D, typename BasisFunctionType>
void MeshPartition<FunctionSpace::FunctionSpace<Mesh::CompositeOfDimension<D>,BasisFunctionType>,Mesh::CompositeOfDimension<D>>::
initializeDofNosLocalNaturalOrdering(std::shared_ptr<FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>> functionSpace);

//! Get a vector of local dof nos in local natural ordering, initializeDofNosLocalNaturalOrdering has to be called beforehand.
template<int D, typename BasisFunctionType>
const std::vector<dof_no_t> &MeshPartition<FunctionSpace::FunctionSpace<Mesh::CompositeOfDimension<D>,BasisFunctionType>,Mesh::CompositeOfDimension<D>>::
dofNosLocalNaturalOrdering() const;

//! check if the given dof is owned by the own rank, then return true, if not, neighbourRankNo is set to the rank by which the dof is owned
template<int D, typename BasisFunctionType>
bool MeshPartition<FunctionSpace::FunctionSpace<Mesh::CompositeOfDimension<D>,BasisFunctionType>,Mesh::CompositeOfDimension<D>>::
isNonGhost(node_no_t nodeNoLocal, int &neighbourRankNo) const;

//! get information about neighbouring rank and boundary elements for specified face,
//! @param neighbourRankNo: the rank of the neighbouring process that shares the face, @param nElements: Size of one-layer mesh that contains boundary elements that touch the neighbouring process
template<int D, typename BasisFunctionType>
void MeshPartition<FunctionSpace::FunctionSpace<Mesh::CompositeOfDimension<D>,BasisFunctionType>,Mesh::CompositeOfDimension<D>>::
getBoundaryElements(Mesh::face_t face, int &neighbourRankNo, std::array<element_no_t,MeshType::dim()> &nBoundaryElements, std::vector<dof_no_t> &dofNos);

//! get the rank no of the neighbour in direction face, -1 if there is no such neighbour
template<int D, typename BasisFunctionType>
int MeshPartition<FunctionSpace::FunctionSpace<Mesh::CompositeOfDimension<D>,BasisFunctionType>,Mesh::CompositeOfDimension<D>>::
neighbourRank(Mesh::face_t face);

//! get the partitioning index in the coordinate direction, i.e. the no. of this rank in this direction, the total number of ranks in each direction can be retrieved by nRanks
template<int D, typename BasisFunctionType>
int MeshPartition<FunctionSpace::FunctionSpace<Mesh::CompositeOfDimension<D>,BasisFunctionType>,Mesh::CompositeOfDimension<D>>::
ownRankPartitioningIndex(int coordinateDirection);

//! refine the partitioning by multiplying the number of elements by refinementFactor
template<int D, typename BasisFunctionType>
void MeshPartition<FunctionSpace::FunctionSpace<Mesh::CompositeOfDimension<D>,BasisFunctionType>,Mesh::CompositeOfDimension<D>>::
refine(std::array<int,MeshType::dim()> refinementFactor);
*/

}  // namespace

