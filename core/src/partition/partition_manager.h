#pragma once

#include <Python.h>  // has to be the first included header
#include <memory>

#include "partition/partition.h"

namespace Partition
{

class Manager
{
public:
  
  //! constructor
  Manager(PyObject *settings);
 
  //! create new partitioning over all available processes, respective the rank subset that was set by the last call to setRankSubsetForNextCreatedMesh
  template<typename BasisOnMesh>
  std::shared_ptr<MeshPartition<BasisOnMesh>> createPartitioning(global_no_t globalSize);

  //! store a rank subset that will be used for the next partitioning that will be created
  void setRankSubsetForNextCreatedMesh(std::shared_ptr<RankSubset> nextRankSubset);
  
  //! get the number of MPI ranks 
  int nRanks();
  
  //! get the own MPI rank no
  int rankNo();
  
private:
 
  PyObject *specificSettings_;  ///< the settings object for the partition manager
 
  int nRanksCommWorld_;  ///< number of MPI ranks, processes
  int rankNoCommWorld_;   ///< own process no.
  
  std::shared_ptr<RankSubset> nextRankSubset_;
};

};    // namespace
