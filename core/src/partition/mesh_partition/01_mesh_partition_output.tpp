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
