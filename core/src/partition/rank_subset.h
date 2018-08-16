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
  RankSubset(std::vector<int> ranks);
 
  //! constructor that constructs a rank subset with all ranks (MPICommWorld)
  RankSubset();
 
  //! number of ranks in the current rank list
  element_no_t size();

  //! first entry of the rank list
  std::vector<int>::const_iterator begin();
  
  //! one after last  entry of the rank list
  std::vector<int>::const_iterator end();
  
  //! get the MPI communicator that contains all ranks of this subset
  MPI_Comm mpiCommunicator();
  
protected:
 
  std::vector<int> rankNo_;  ///< the list of ranks
  MPI_Comm mpiCommunicator_;    ///< the MPI communicator that contains only the ranks of this rank subset
  
};

//! output rank subset
std::ostream &operator<<(std::ostream &stream, RankSubset rankSubset);

}  // namespace
