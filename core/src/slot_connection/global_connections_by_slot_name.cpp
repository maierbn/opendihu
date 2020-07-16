#include "slot_connection/global_connections_by_slot_name.h"

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
    s << "The following data slot connection were given by the setting \"connectedSlots\":\n";
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
      if (connectionAlreadyDisplayed)
        continue;

      std::string fromSlotName = (*iter).first;
      std::string toSlotName = (*iter).second;

      if (fromSlotName.length() < 6)
      {
        fromSlotName = std::string(6-fromSlotName.length(), ' ') + fromSlotName;
      }
      if (toSlotName.length() < 6)
      {
        toSlotName += std::string(6-toSlotName.length(), ' ');
      }

      if (reverseConnectionAlsoContained)
      {
        s << "  " << fromSlotName << " ¤ <─> ¤ " << toSlotName << "\n";
      }
      else
      {
        s << "  " << fromSlotName << " ¤ ──> ¤ " << toSlotName << "\n";
      }
    }
    s << "\n";
  }

  if (!connectionsFromSameSlotNames_.empty())
  {
    s << "The following data slots were connected because the names appeared in both terms of a coupling or splitting scheme:\n";

    for (std::string slotName : connectionsFromSameSlotNames_)
    {
      std::string fromSlotName = slotName;
      std::string toSlotName = slotName;

      if (fromSlotName.length() < 6)
      {
        fromSlotName = std::string(6-fromSlotName.length(), ' ') + fromSlotName;
      }
      if (toSlotName.length() < 6)
      {
        toSlotName += std::string(6-toSlotName.length(), ' ');
      }

      s << "  " << fromSlotName << " ¤ <─> ¤ " << toSlotName << "\n";
    }
    s << "\n";
  }

  return s.str();
}
