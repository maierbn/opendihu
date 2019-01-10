#pragma once

#include <set>

#include "control/types.h"

namespace Partition
{

/** This is a list of ranks that perform a task, e.g. compute a partition.
 */
class RankSubset
{
public:

  //! constructor that constructs a rank subset with a single rank
  RankSubset(int singleRank, MPI_Comm mpiCommunicator);

  //! constructor that constructs a rank subset with a single rank, from MPICommWorld
  RankSubset(int singleRank);
  
  //! constructor that constructs a whole set of ranks
  template<typename Iter>
  RankSubset(Iter ranksBegin, Iter ranksEnd);
 
  //! constructor that constructs a rank subset with all ranks (MPICommWorld)
  RankSubset();


  //! number of ranks in the current rank list
  element_no_t size() const;

  //! get the own rank id of the mpi Communicator
  element_no_t ownRankNo();

  //! check if the own rank from MPICommWorld is contained in the current rankSubset
  bool ownRankIsContained() const;

  //! first entry of the rank list
  std::set<int>::const_iterator begin();
  
  //! one after last  entry of the rank list
  std::set<int>::const_iterator end();
  
  //! get the MPI communicator that contains all ranks of this subset
  MPI_Comm mpiCommunicator() const;
  
protected:
 
  std::set<int> rankNo_;  ///< the list of ranks
  int ownRankNo_;             ///< own rank id of this rankSubset
  MPI_Comm mpiCommunicator_;    ///< the MPI communicator that contains only the ranks of this rank subset
  
};

//! output rank subset
std::ostream &operator<<(std::ostream &stream, RankSubset rankSubset);

}  // namespace

#include "partition/rank_subset.tpp"
