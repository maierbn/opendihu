#pragma once

#include <Python.h>  // has to be the first included header
#include <map>
#include <vector>

#include "partition/rank_subset.h"

namespace Control
{

/** Helper class to send double values from any rank to any rank.
 *  Furthermore, "dofNos", i.e. arbitrary integers can be send the same way as the values on initialize.
 */
class ValueCommunicator
{
public:

  //! Initialize the communication, this needs to be called once, before the first communication.
  //! This allows to send indices in remoteDofNosAtRanks to the ranks, they will appear there in dofNos.
  //! @param remoteDofNosLocalAtRanks[in] The rank nos to send to and the local dofs at these ranks, where the values should be stored remotely
  //! @param dofNosLocal[out] The local dofs where the received values should be stored.
  //! @param rankSubset[in] the ranks that are involved in communication
  void initialize(const std::map<int,std::vector<int>> remoteDofNosAtRanks, std::vector<int> &dofNos, std::shared_ptr<Partition::RankSubset> rankSubset);

  //! communicate the values to the ranks
  void communicate(const std::map<int,std::vector<double>> valuesToSendToRanks, std::vector<double> &receivedValues);

protected:

  std::shared_ptr<Partition::RankSubset> rankSubset_;    //< the involved ranks

  std::vector<std::pair<int,int>> nDofsToReceiveFromRanks_;   //< (foreignRank,nNodes), number of nodes requested by and to be send to foreignRank
  int nValuesToReceive_;               //< total number of values to receive on own rank
};

}  // namespace Control
