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

    LOG(DEBUG) << "connectedSlotNames " << connectedSlotNames;

    // if the "from" and "to" slot names match to slot names in the current slotConnectorData
    if (iterFrom != slotNamesTerm1.end() && iterTo != slotNamesTerm2.end())
    {
      // determine slot nos
      int slotNoFrom = std::distance(slotNamesTerm1.begin(), iterFrom);
      int slotNoTo = std::distance(slotNamesTerm2.begin(), iterTo);

      LOG(DEBUG) << "Term1->Term2, slots " << slotNoFrom << " -> " << slotNoTo;

      // add output connection
      slotsConnection->addConnectionTerm1ToTerm2(slotNoFrom, slotNoTo);
    }

    // direction Term2 -> Term1

    // determine slot no.s for current slotConnectorData
    iterFrom = std::find(slotNamesTerm2.begin(), slotNamesTerm2.end(), connectedSlotNames.first);
    iterTo   = std::find(slotNamesTerm1.begin(), slotNamesTerm1.end(), connectedSlotNames.second);

    // if the "from" and "to" slot names match to slot names in the current slotConnectorData
    if (iterFrom != slotNamesTerm2.end() && iterTo != slotNamesTerm1.end())
    {
      // determine slot nos
      int slotNoFrom = std::distance(slotNamesTerm2.begin(), iterFrom);
      int slotNoTo = std::distance(slotNamesTerm1.begin(), iterTo);

      LOG(DEBUG) << "Term2->Term1, slots " << slotNoFrom << " -> " << slotNoTo;

      // add output connection
      slotsConnection->addConnectionTerm2ToTerm1(slotNoFrom, slotNoTo);
    }
  }
}
