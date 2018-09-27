#include "partition/01_mesh_partition.h"

#include <petscsys.h>

#include "mesh/mesh.h"
#include "utility/mpi_utility.h"


// output local to global mapping
std::ostream &operator<<(std::ostream &stream, std::shared_ptr<ISLocalToGlobalMapping> localToGlobalMapping)
{
  stream << *localToGlobalMapping;
  return stream;
}

// output local to global mapping
std::ostream &operator<<(std::ostream &stream, ISLocalToGlobalMapping localToGlobalMapping)
{
  PetscErrorCode ierr;
  PetscInt nproc;
  PetscInt *procs;
  PetscInt *numprocs;
  PetscInt **indices;
  ierr = ISLocalToGlobalMappingGetInfo(localToGlobalMapping, &nproc, &procs, &numprocs, &indices); CHKERRABORT(MPI_COMM_WORLD, ierr);
  stream << "[localToGlobalMapping, " << nproc << " ranks: (";
  
  // output ranks
  if (nproc > 0)
  {
    stream << procs[0];
  }
  for (int i = 1; i < nproc; i++)
  {
    stream << ", " << procs[i];
  }
  stream << "), indices of nodes (in local numbering) shared with neighbors (sorted by global numbering) : (";
  
  // output boundary elements
  if (nproc > 0)
  {
    stream << "0=";
    for (int j = 0; j < numprocs[0]; j++) 
    {
      stream << indices[0][j] << " ";
    }
  }
  for (int i = 1; i < nproc; i++)
  {
    stream << "; " << i << "=";
    for (int j = 0; j < numprocs[i]; j++) 
    {
      stream << indices[i][j] << " ";
    }
  }
  stream << "), ";
  
  PetscInt nLocalIndices;
  ierr = ISLocalToGlobalMappingGetSize(localToGlobalMapping, &nLocalIndices); CHKERRABORT(MPI_COMM_WORLD, ierr);
  stream << nLocalIndices << " local to Petsc: [";
  
  PetscInt const *localIndices;
  ierr = ISLocalToGlobalMappingGetIndices(localToGlobalMapping, &localIndices); CHKERRABORT(MPI_COMM_WORLD, ierr);
  
  if (nLocalIndices > 0)
  {
    stream << "l0=p" << localIndices[0];
  }

  int indicesEnd = std::min(100, nLocalIndices);
  if (VLOG_IS_ON(1))
  {
    indicesEnd = nLocalIndices;
  }

  for (int i = 1; i < indicesEnd; i++)
  {
    stream << ",l" << i << "=p" << localIndices[i];
  }
  if (indicesEnd < nLocalIndices)
  {
    stream << " ... (" << nLocalIndices << " local indices)";
  }
  stream << "]]";
  
  return stream;
}
