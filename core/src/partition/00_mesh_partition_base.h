#pragma once

#include <memory>
#include <petscdmda.h>

#include "control/types.h"
#include "partition/rank_subset.h"

namespace Partition
{

/** base class for mesh partition */
class MeshPartitionBase
{
public:
 
  //! constructor, store the rankSubset
  MeshPartitionBase(std::shared_ptr<RankSubset> rankSubset);
  
  //! virtual destructor
  virtual ~MeshPartitionBase();
 
  //! number of ranks
  int nRanks();
  
  //! number of entries in the current partition (this usually refers to the elements)
  virtual element_no_t nElementsLocal() = 0;
  
  //! number of nodes in total
  virtual global_no_t nElementsGlobal() = 0;
  
  //! return reference to a vector containing all local dofs, i.e. a vector with {0,1,2,...,nLocalDofsWithGhosts-1}
  std::vector<PetscInt> &localDofNos();
  
  //! remove all dofs from the vector that are not handled in the local partition
  virtual void extractLocalDofsWithoutGhosts(std::vector<double> &values) = 0;
  
  //! get the MPI communicator that is needed for the work portion
  MPI_Comm mpiCommunicator();
  
protected:
   
  //! initialize the localNodeNos_ vector to values {0,1,...,nLocalDofsWithGhosts-1}
  void initializeLocalDofsVector(node_no_t nLocalDofsWithGhosts);
 
  //TODO localDofNos wg Hermite
  std::vector<PetscInt> localDofNos_;     ///< list of local dofs for a field variable with 1 component (1 dof per node). This is {0,1,...,nLocalDofsWithGhosts()-1}. It is needed for calls to Petsc functions that access all local data, e.g. within fieldVariable->getValues
  std::shared_ptr<RankSubset> rankSubset_;  ///< the set of ranks that compute something where this partition is a part of, also holds the MPI communciator
};

}  // namespace
