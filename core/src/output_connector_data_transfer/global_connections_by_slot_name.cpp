#include "output_connector_data_transfer/global_connections_by_slot_name.h"

#include "utility/python_utility.h"

GlobalConnectionsBySlotName::
GlobalConnectionsBySlotName(PythonConfig settings)
{
  // parse the connectedSlots
  if (settings.hasKey("connectedSlots"))
  {
    settings.getOptionVector<std::pair<std::string,std::string>>("connectedSlots", connections_);
  }
}

std::string GlobalConnectionsBySlotName::getDescriptionForDiagram()
{
  std::stringstream s;

  if (!connections_.empty())
  {
    s << "The following data slot connection were given by the setting \"connectionSlots\":\n";
    for (std::vector<std::pair<std::string,std::string>>::iterator iter = connections_.begin(); iter != connections_.end(); iter++)
    {
      // check if reverse connection is also contained
      bool reverseConnectionAlsoContained = false;
      bool connectionAlreadyDisplayed = false;

      for (std::vector<std::pair<std::string,std::string>>::iterator iter2 = connections_.begin(); iter2 != connections_.end(); iter2++)
      {
        if (iter->first == iter2->second && iter->second == iter2->first)
        {
          reverseConnectionAlsoContained = true;
          if (iter2 < iter)
            connectionAlreadyDisplayed = true;

          break;
        }
      }

      if (reverseConnectionAlsoContained)
      {
        s << "  " << (*iter).first << "¤ <─> ¤" << (*iter).second << "\n";
      }
      else if (!connectionAlreadyDisplayed)
      {
        s << "  " << (*iter).first << "¤ ─> ¤" << (*iter).second << "\n";
      }
    }
  }
  return s.str();
}
