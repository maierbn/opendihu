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
  int nRanks() const;

  //! own rank id in the current communcator
  int ownRankNo();
  
  //! number of entries in the current partition (this usually refers to the elements)
  virtual element_no_t nElementsLocal() const = 0;
  
  //! number of nodes in total
  virtual global_no_t nElementsGlobal() const = 0;
  
  //! remove all dofs from the vector that are not handled in the local partition
  virtual void extractLocalDofsWithoutGhosts(std::vector<double> &values) const = 0;
  
  //! get the MPI communicator that is needed for the work portion
  MPI_Comm mpiCommunicator() const;
  
  //! fill the dofNosLocal vector
  void createLocalDofOrderings(dof_no_t nDofsLocal);
  
  //! get a vector of local dof nos, range [0,nDofsLocalWithoutGhosts] are the dofs without ghost dofs, the whole vector are the dofs with ghost dofs (only for structured mesh)e
  const std::vector<PetscInt> &dofNosLocal() const;
  
  //! number of dofs in total
  virtual global_no_t nDofsGlobal() const = 0;

  //! get a vector of local dof nos including ghost dofs, in the natural ordering
  virtual const std::vector<dof_no_t> &dofNosLocalNaturalOrdering() const;

  //! get a vector of global natural dof nos of the locally stored non-ghost dofs, needed for setParameters callback function in cellml adapter
  virtual void getDofNosGlobalNatural(std::vector<global_no_t> &dofNosGlobalNatural) const = 0;

  //! get a PETSc IS (index set) with the same information as dofNosLocal_
  const IS &dofNosLocalIS() const;
  
protected:
   
  std::shared_ptr<RankSubset> rankSubset_;  ///< the set of ranks that compute something where this partition is a part of, also holds the MPI communciator
  
  std::vector<dof_no_t> dofNosLocal_;   ///< vector of all local nos of non-ghost dofs followed by the ghost dofs
  IS dofNosLocalIS_;   ///< index set (IS) with the indices of the local dof nos (including ghosts)
};

}  // namespace
