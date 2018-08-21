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
  
  //! number of entries in the current partition (this usually refers to the elements)
  virtual element_no_t nElementsLocal() const = 0;
  
  //! number of nodes in total
  virtual global_no_t nElementsGlobal() const = 0;
  
  //! remove all dofs from the vector that are not handled in the local partition
  virtual void extractLocalDofsWithoutGhosts(std::vector<double> &values) const = 0;
  
  //! get the MPI communicator that is needed for the work portion
  MPI_Comm mpiCommunicator() const;
  
protected:
   
  std::shared_ptr<RankSubset> rankSubset_;  ///< the set of ranks that compute something where this partition is a part of, also holds the MPI communciator
};

}  // namespace
