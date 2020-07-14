#include "output_connector_data_transfer/global_connections_by_slot_name.h"

#include "utility/python_utility.h"
#include "output_connector_data_transfer/output_connector_data_helper.h"

//! add all parsed slot connections to the outputConnection_ object of a splitting scheme
template<typename OutputConnectorDataType1, typename OutputConnectorDataType2>
void GlobalConnectionsBySlotName::
addConnections(
  std::shared_ptr<std::tuple<std::shared_ptr<OutputConnectorDataType1>,std::shared_ptr<OutputConnectorDataType2>>> outputConnectorData,
  std::shared_ptr<OutputConnection> outputConnection
)
{
  if (!outputConnectorData)
    return;

  std::vector<std::string> slotNamesTerm1;
  std::vector<std::string> slotNamesTerm2;
  OutputConnectorDataHelper<OutputConnectorDataType1>::getSlotNames(std::get<0>(*outputConnectorData), slotNamesTerm1);
  OutputConnectorDataHelper<OutputConnectorDataType2>::getSlotNames(std::get<1>(*outputConnectorData), slotNamesTerm2);

  // loop over global connections
  for (std::pair<std::string,std::string> connectedSlotNames : connections_)
  {
    // direction Term1 -> Term2

    // determine slot no.s for current outputConnectorData
    std::vector<std::string>::iterator iterFrom = std::find(slotNamesTerm1.begin(), slotNamesTerm1.end(), connectedSlotNames.first);
    std::vector<std::string>::iterator iterTo   = std::find(slotNamesTerm2.begin(), slotNamesTerm2.end(), connectedSlotNames.second);

    // if the "from" and "to" slot names match to slot names in the current outputConnectorData
    if (iterFrom != slotNamesTerm1.end() && iterTo != slotNamesTerm2.end())
    {
      // determine slot nos
      int slotNoFrom = std::distance(slotNamesTerm1.begin(), iterFrom);
      int slotNoTo = std::distance(slotNamesTerm2.begin(), iterTo);

      // add output connection
      outputConnection->addConnectionTerm1ToTerm2(slotNoFrom, slotNoTo);
    }

    // direction Term2 -> Term1

    // determine slot no.s for current outputConnectorData
    iterFrom = std::find(slotNamesTerm2.begin(), slotNamesTerm2.end(), connectedSlotNames.first);
    iterTo   = std::find(slotNamesTerm1.begin(), slotNamesTerm1.end(), connectedSlotNames.second);

    // if the "from" and "to" slot names match to slot names in the current outputConnectorData
    if (iterFrom != slotNamesTerm2.end() && iterTo != slotNamesTerm1.end())
    {
      // determine slot nos
      int slotNoFrom = std::distance(slotNamesTerm2.begin(), iterFrom);
      int slotNoTo = std::distance(slotNamesTerm1.begin(), iterTo);

      // add output connection
      outputConnection->addConnectionTerm2ToTerm1(slotNoFrom, slotNoTo);
    }
  }
}
