#include "partition/mesh_partition/01_mesh_partition.h"

#include <petscsys.h>

#include "mesh/mesh.h"
#include "utility/mpi_utility.h"

//! output mesh partition
template<typename FunctionSpaceType>
std::ostream &operator<<(std::ostream &stream, std::shared_ptr<Partition::MeshPartition<FunctionSpaceType>> meshPartition)
{
  meshPartition->output(stream);
  return stream;
}


//! output mesh partition
template<typename FunctionSpaceType>
std::ostream &operator<<(std::ostream &stream, Partition::MeshPartition<FunctionSpaceType> meshPartition)
{
  meshPartition.output(stream);
  return stream;
}

namespace Partition
{

template<typename MeshType,typename BasisFunctionType>
void MeshPartition<FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>,Mesh::isStructured<MeshType>>::
output(std::ostream &stream)
{
  stream << "MeshPartition<structured>, nElementsGlobal: ";
  for (int i = 0; i < MeshType::dim(); i++)
    stream << nElementsGlobal_[i] << ",";

  stream << " nElementsLocal: ";
  for (int i = 0; i < MeshType::dim(); i++)
    stream << nElementsLocal_[i] << ",";

  stream << " nRanks: ";
  for (int i = 0; i < MeshType::dim(); i++)
    stream << nRanks_[i] << ",";

  stream << " beginElementGlobal: ";
  for (int i = 0; i < MeshType::dim(); i++)
    stream << beginElementGlobal_[i] << ",";

  stream << " hasFullNumberOfNodes: " ;
  for (int i = 0; i < MeshType::dim(); i++)
    stream << hasFullNumberOfNodes_[i] << ",";

  stream << " localSizesOnPartitions: ";
  for (int i = 0; i < MeshType::dim(); i++)
  {
    stream << "(";
    for (int j = 0; j < localSizesOnPartitions_[i].size(); j++)
      stream << localSizesOnPartitions_[i][j] << ",";
    stream << ")";
  }

  stream << " nNodesGlobal: ";
  for (int i = 0; i < MeshType::dim(); i++)
    stream << nNodesGlobal(i) << ",";
  stream << "total " << nNodesGlobal() << ",";

  stream << " nNodesLocalWithGhosts: ";
  for (int i = 0; i < MeshType::dim(); i++)
    stream << nNodesLocalWithGhosts(i) << ",";
  stream << "total " << nNodesLocalWithGhosts() << ",";

  stream << " nNodesLocalWithoutGhosts: ";
  for (int i = 0; i < MeshType::dim(); i++)
    stream << nNodesLocalWithoutGhosts(i) << ",";
  stream << "total " << nNodesLocalWithoutGhosts()
    << ", dofNosLocal: [";

  int dofNosLocalEnd = std::min(100, (int)this->dofNosLocal_.size());
  if (VLOG_IS_ON(2))
  {
    dofNosLocalEnd = this->dofNosLocal_.size();
  }
  for (int i = 0; i < dofNosLocalEnd; i++)
  {
    stream << this->dofNosLocal_[i] << " ";
  }
  if (dofNosLocalEnd < this->dofNosLocal_.size())
    stream << " ... ( " << this->dofNosLocal_.size() << " local dof nos)";
  stream << "], ghostDofNosGlobalPetsc: [";

  for (int i = 0; i < std::min(100,(int)ghostDofNosGlobalPetsc_.size()); i++)
    stream << ghostDofNosGlobalPetsc_[i] << " ";
  stream << "], localToGlobalPetscMappingDofs_: ";
  if (localToGlobalPetscMappingDofs_)
  {
    stream << localToGlobalPetscMappingDofs_;
  }
  else
  {
    stream << "(none)";
  }
}

} // namespace
