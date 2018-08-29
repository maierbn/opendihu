#pragma once

#include <vector>

#include "control/types.h"

namespace Partition
{

/** This is a list of ranks that perform a task, e.g. compute a partition.
 */
class RankSubset
{
public:
  
  //! constructor that constructs a rank subset with a single rank
  RankSubset(int singleRank);
  
  //! constructor that constructs a whole vector of ranks
  RankSubset(std::vector<int> &ranks);
 
  //! constructor that constructs a rank subset with all ranks (MPICommWorld)
  RankSubset();
 
  //! number of ranks in the current rank list
  element_no_t size() const;

  //! get the own rank id of the mpi Communicator
  element_no_t ownRankNo();

  //! first entry of the rank list
  std::vector<int>::const_iterator begin();
  
  //! one after last  entry of the rank list
  std::vector<int>::const_iterator end();
  
  //! get the MPI communicator that contains all ranks of this subset
  MPI_Comm mpiCommunicator() const;
  
protected:
 
  std::vector<int> rankNo_;  ///< the list of ranks
  int ownRankNo_;             ///< own rank id of this rankSubset
  MPI_Comm mpiCommunicator_;    ///< the MPI communicator that contains only the ranks of this rank subset
  
};

//! output rank subset
std::ostream &operator<<(std::ostream &stream, RankSubset rankSubset);

}  // namespace
