#include "slot_connection/global_connections_by_slot_name.h"

#include "utility/python_utility.h"
#include "slot_connection/slot_connector_data_helper.h"

//! add all parsed slot connections to the slotsConnection_ object of a splitting scheme
template<typename SlotConnectorDataType1, typename SlotConnectorDataType2>
void GlobalConnectionsBySlotName::
addConnections(
  std::shared_ptr<std::tuple<std::shared_ptr<SlotConnectorDataType1>,std::shared_ptr<SlotConnectorDataType2>>> slotConnectorData,
  std::shared_ptr<SlotsConnection> slotsConnection
)
{
  if (!slotConnectorData)
  {
    LOG(ERROR) << "In GlobalConnectionsBySlotName::addConnections, slotConnectorData is not set.";
    return;
  }

  if (!slotsConnection)
  {
    LOG(ERROR) << "In GlobalConnectionsBySlotName::addConnections, slotsConnection is not set.";
    return;
  }

  std::vector<std::string> slotNamesTerm1;
  std::vector<std::string> slotNamesTerm2;
  SlotConnectorDataHelper<SlotConnectorDataType1>::getSlotNames(std::get<0>(*slotConnectorData), slotNamesTerm1);
  SlotConnectorDataHelper<SlotConnectorDataType2>::getSlotNames(std::get<1>(*slotConnectorData), slotNamesTerm2);

  // loop over global connections
  for (std::pair<std::string,std::string> connectedSlotNames : connections_)
  {
    // direction Term1 -> Term2

    // determine slot no.s for current slotConnectorData
    std::vector<std::string>::iterator iterFrom = std::find(slotNamesTerm1.begin(), slotNamesTerm1.end(), connectedSlotNames.first);
    std::vector<std::string>::iterator iterTo   = std::find(slotNamesTerm2.begin(), slotNamesTerm2.end(), connectedSlotNames.second);

    LOG(DEBUG) << "connectedSlotNames " << connectedSlotNames << " slotNamesTerm1: " << slotNamesTerm1 << ", slotNamesTerm2: " << slotNamesTerm2;
    LOG(DEBUG) << "  Term1->Term2, Term1 found: " << (iterFrom != slotNamesTerm1.end()) << ", Term2 found: " << (iterTo != slotNamesTerm2.end());

    // if the "from" and "to" slot names match to slot names in the current slotConnectorData
    if (iterFrom != slotNamesTerm1.end() && iterTo != slotNamesTerm2.end())
    {
      // determine slot nos
      int slotNoFrom = std::distance(slotNamesTerm1.begin(), iterFrom);
      int slotNoTo = std::distance(slotNamesTerm2.begin(), iterTo);

      LOG(DEBUG) << "Term1->Term2, add connectionn for slots " << slotNoFrom << " -> " << slotNoTo;

      // add output connection
      slotsConnection->addConnectionTerm1ToTerm2(slotNoFrom, slotNoTo);
    }

    // direction Term2 -> Term1

    // determine slot no.s for current slotConnectorData
    iterFrom = std::find(slotNamesTerm2.begin(), slotNamesTerm2.end(), connectedSlotNames.first);
    iterTo   = std::find(slotNamesTerm1.begin(), slotNamesTerm1.end(), connectedSlotNames.second);

    LOG(DEBUG) << "  Term2->Term1, Term2 found: " << (iterFrom != slotNamesTerm2.end()) << ", Term1 found: " << (iterTo != slotNamesTerm1.end());

    // if the "from" and "to" slot names match to slot names in the current slotConnectorData
    if (iterFrom != slotNamesTerm2.end() && iterTo != slotNamesTerm1.end())
    {
      // determine slot nos
      int slotNoFrom = std::distance(slotNamesTerm2.begin(), iterFrom);
      int slotNoTo = std::distance(slotNamesTerm1.begin(), iterTo);

      LOG(DEBUG) << "Term2->Term1, add connectionn for slots " << slotNoFrom << " -> " << slotNoTo;

      // add output connection
      slotsConnection->addConnectionTerm2ToTerm1(slotNoFrom, slotNoTo);
    }
  }

  // check for slots that have the same slot name in both Term1 and Term2, these will also be connected
  for (int slotNoTerm1 = 0; slotNoTerm1 < slotNamesTerm1.size(); slotNoTerm1++)
  {
    std::string slotNameTerm1 = slotNamesTerm1[slotNoTerm1];

    if (slotNameTerm1 == "")
      continue;

    std::vector<std::string>::iterator iter = std::find(slotNamesTerm2.begin(), slotNamesTerm2.end(), slotNameTerm1);

    if (iter != slotNamesTerm2.end())
    {
      int slotNoTerm2 = std::distance(slotNamesTerm2.begin(), iter);

      LOG(DEBUG) << "add two-way connections because of matching slot name " << slotNameTerm1;

      // if connection does not yet exist
      if (std::find(connectionsFromSameSlotNames_.begin(), connectionsFromSameSlotNames_.end(), slotNameTerm1) == connectionsFromSameSlotNames_.end())
      {
        // add both connections
        slotsConnection->addConnectionTerm1ToTerm2(slotNoTerm1, slotNoTerm2);
        slotsConnection->addConnectionTerm2ToTerm1(slotNoTerm2, slotNoTerm1);

        connectionsFromSameSlotNames_.push_back(slotNameTerm1);
      }
    }
  }
}
