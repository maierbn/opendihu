#include "output_connector_data_transfer/output_connection.h"

#include "control/diagnostic_tool/solver_structure_visualizer.h"

OutputConnection::OutputConnection(PythonConfig settings):
  fieldVariableNamesInitialized_(false), slotInformationInitialized_(false), settings_(settings)
{
  // parse values from settings
  std::vector<int> indicesTerm1To2;
  std::vector<int> indicesTerm2To1;

  settings.getOptionVector("connectedSlotsTerm1To2", indicesTerm1To2);
  settings.getOptionVector("connectedSlotsTerm2To1", indicesTerm2To1);

  // set indices in connectorTerm1To2_ and connectorTerm2To1_
  for (std::vector<int>::iterator iter = indicesTerm1To2.begin(); iter != indicesTerm1To2.end(); iter++)
  {
    Connector connector;
    connector.index = *iter;
    connector.avoidCopyIfPossible = false;
    connectorTerm1To2_.push_back(connector);
  }

  for (std::vector<int>::iterator iter = indicesTerm2To1.begin(); iter != indicesTerm2To1.end(); iter++)
  {
    Connector connector;
    connector.index = *iter;
    connector.avoidCopyIfPossible = false;
    connectorTerm2To1_.push_back(connector);
  }

  // if field variable gets mapped in both directions, set avoidCopyIfPossible to true
  for (int i = 0; i < connectorTerm1To2_.size(); i++)
  {
    int mappedIndex = connectorTerm1To2_[i].index;

    // check if other direction is mapped the same way
    if (connectorTerm2To1_.size() > mappedIndex)
    {
      if (connectorTerm2To1_[mappedIndex].index == i)
      {
        connectorTerm1To2_[i].avoidCopyIfPossible = true;
        connectorTerm2To1_[i].avoidCopyIfPossible = true;
      }
    }
  }

}

void OutputConnection::setTransferDirection(bool term1To2)
{
  transferDirectionTerm1To2_ = term1To2;
}


//! get the connectors from term 1 to term 2
const std::vector<OutputConnection::Connector> &OutputConnection::connectorTerm1To2() const
{
  return connectorTerm1To2_;
}

//! get the connectors from term 2 to term 1
const std::vector<OutputConnection::Connector> &OutputConnection::connectorTerm2To1() const
{
  return connectorTerm2To1_;
}

std::string OutputConnection::getDebugInformation() const
{
  std::stringstream result;

  if (transferDirectionTerm1To2_)
  {
    result << "Transfer data from Term1 to Term2." << std::endl;
  }
  else
  {
    result << "Transfer data from Term2 to Term1." << std::endl;
  }

  // output slots of Term 1
  int nSlotsTerm1 = fieldVariableNamesTerm1Vector1_.size() + fieldVariableNamesTerm1Vector2_.size();
  result << "Term 1 has the following " << nSlotsTerm1 << " slot"
    << (nSlotsTerm1 == 1? "" : "s") << ":" << std::endl;

  if (nSlotsTerm1 == 0)
    result << "  (none)" << std::endl;

  int slotNo = 0;
  for (int i = 0; i < fieldVariableNamesTerm1Vector1_.size(); i++, slotNo++)
  {
    result << "  slot " << slotNo << ": " << fieldVariableNamesTerm1Vector1_[i] << " (stored in variable1[" << i << "])" << std::endl;
  }
  for (int i = 0; i < fieldVariableNamesTerm1Vector2_.size(); i++, slotNo++)
  {
    result << "  slot " << slotNo << ": " << fieldVariableNamesTerm1Vector2_[i] << " (stored in variable2[" << i << "])" << std::endl;
  }

  // output slots of Term 2
  int nSlotsTerm2 = fieldVariableNamesTerm2Vector1_.size() + fieldVariableNamesTerm2Vector2_.size();
  result << "Term 2 has the following " << nSlotsTerm2 << " slot"
    << (nSlotsTerm2 == 1? "" : "s") << ":" << std::endl;

  if (nSlotsTerm2 == 0)
    result << "  (none)" << std::endl;

  slotNo = 0;
  for (int i = 0; i < fieldVariableNamesTerm2Vector1_.size(); i++, slotNo++)
  {
    result << "  slot " << slotNo << ": " << fieldVariableNamesTerm2Vector1_[i] << " (stored in variable1[" << i << "])" << std::endl;
  }
  for (int i = 0; i < fieldVariableNamesTerm2Vector2_.size(); i++, slotNo++)
  {
    result << "  slot " << slotNo << ": " << fieldVariableNamesTerm2Vector2_[i] << " (stored in variable2[" << i << "])" << std::endl;
  }

  // output connections of Term 1
  result << "Term 1 has the following " << connectorTerm1To2_.size() << " connection"
    << (connectorTerm1To2_.size() == 1? "" : "s") << ":" << std::endl;

  if (connectorTerm1To2_.size() == 0)
  result << "  (none)" << std::endl;

  for (int i = 0; i < connectorTerm1To2_.size(); i++)
  {
    result << "  Term1.slot " << i << " (";

    int vectorNo = std::min(1,i/nFieldVariablesTerm1Vector1_);
    int vectorIndex = i - vectorNo*nFieldVariablesTerm1Vector1_;

    result << "Term1.variable" << vectorNo+1 << "[" << vectorIndex << "], ";

    if (vectorNo == 0)
    {
      result << fieldVariableNamesTerm1Vector1_[vectorIndex];
    }
    else
    {
      if (vectorIndex < fieldVariableNamesTerm1Vector2_.size())
      {
        result << fieldVariableNamesTerm1Vector2_[vectorIndex];
      }
      else
      {
        result << "(unknown field variable)";
      }
    }

    result << ") -> Term2.slot " << connectorTerm1To2_[i].index << " (";

    if (connectorTerm1To2_[i].index == -1)
    {
      result << "None)";
    }
    else
    {
      vectorNo = std::min(1,connectorTerm1To2_[i].index/nFieldVariablesTerm2Vector1_);
      vectorIndex = connectorTerm1To2_[i].index - vectorNo*nFieldVariablesTerm2Vector1_;

      result << "Term2.variable" << vectorNo+1 << "[" << vectorIndex << "], ";

      if (vectorNo == 0)
      {
        result << fieldVariableNamesTerm2Vector1_[vectorIndex];
      }
      else
      {
        if (vectorIndex < fieldVariableNamesTerm2Vector2_.size())
        {
          result << fieldVariableNamesTerm2Vector2_[vectorIndex];
        }
        else
        {
          result << "(unknown field variable)";
        }
      }

      if (connectorTerm1To2_[i].avoidCopyIfPossible)
      {
        result << ") avoidCopyIfPossible";
      }
      else
      {
        result << ") data copy";
      }
    }
    result << std::endl;
  }

  // output connections of Term 2
  result << "Term 2 has the following " << connectorTerm2To1_.size() << " connection"
    << (connectorTerm2To1_.size() == 1? "" : "s") << ":" << std::endl;

  if (connectorTerm2To1_.size() == 0)
    result << "  (none)" << std::endl;

  for (int i = 0; i < connectorTerm2To1_.size(); i++)
  {
    result << "  Term2.slot " << i << " (";

    int vectorNo = std::min(1,i/nFieldVariablesTerm2Vector1_);
    int vectorIndex = i - vectorNo*nFieldVariablesTerm2Vector1_;

    result << "Term2.variable" << vectorNo+1 << "[" << vectorIndex << "], ";

    if (vectorNo == 0)
    {
      result << fieldVariableNamesTerm2Vector1_[vectorIndex];
    }
    else
    {
      if (vectorIndex < fieldVariableNamesTerm2Vector2_.size())
      {
        result << fieldVariableNamesTerm2Vector2_[vectorIndex];
      }
      else
      {
        result << "(unknown field variable)";
      }
    }

    result << ") -> Term1.slot " << connectorTerm2To1_[i].index << " (";

    if (connectorTerm2To1_[i].index == -1)
    {
      result << "None)";
    }
    else
    {
      vectorNo = std::min(1,connectorTerm2To1_[i].index/nFieldVariablesTerm1Vector1_);
      vectorIndex = connectorTerm2To1_[i].index - vectorNo*nFieldVariablesTerm1Vector1_;

      result << "Term1.variable" << vectorNo+1 << "[" << vectorIndex << "], ";

      if (vectorNo == 0)
      {
        result << fieldVariableNamesTerm1Vector1_[vectorIndex];
      }
      else
      {
        if (vectorIndex < fieldVariableNamesTerm1Vector2_.size())
        {
          result << fieldVariableNamesTerm1Vector2_[vectorIndex];
        }
        else
        {
          result << "(unknown field variable)";
        }
      }

      if (connectorTerm2To1_[i].avoidCopyIfPossible)
      {
        result << ") avoidCopyIfPossible";
      }
      else
      {
        result << ") data copy";
      }
    }
    result << std::endl;
  }

  return result.str();
}

bool OutputConnection::getSlotInformation(int fromVectorNo, int fromVectorIndex,
                                          int &toVectorNo, int &toVectorIndex, bool &avoidCopyIfPossible, bool disableWarnings) const
{
#ifdef NDEBUG
  if (slotInformationInitialized_)
  {
    const Result &result = slotInformation_[transferDirectionTerm1To2_][fromVectorNo][fromVectorIndex];
    toVectorNo = result.toVectorNo;
    toVectorIndex = result.toVectorIndex;
    avoidCopyIfPossible = result.avoidCopyIfPossible;
    return result.successful;
  }
#endif

  // fromVectorNo and toVectorNo are 0 or 1

#ifndef NDEBUG
  VLOG(1) << "getSlotInformation(" << fromVectorNo << "," << fromVectorIndex << ")" << getDebugInformation();
#endif

  if (transferDirectionTerm1To2_)
  {
    // if we map from term 1 to term 2

    // transform from index
    int fromIndex = fromVectorNo*nFieldVariablesTerm1Vector1_ + fromVectorIndex;

    // check for wrong indices and output errors
    if (fromIndex >= connectorTerm1To2_.size())
    {
      if (!disableWarnings)
        LOG(WARNING) << "(0) Unconnected slot " << fromIndex << " of Term1 in " << settings_ << ": " << getDebugInformation() << ", fromIndex ("
          << fromIndex << "=" << fromVectorNo << "*" << nFieldVariablesTerm1Vector1_ << "+" << fromVectorIndex << ") >= " << connectorTerm1To2_.size()
          << "\nThere are only " << connectorTerm1To2_.size() << " slots to connect to, but a connection for slot no " << fromIndex << " is needed."
          << "\nMaybe not enough entries are given for \"connectedSlotsTerm1To2\" or \"connectedSlotsTerm2To1\" "
          << "or unneccessary slots have been created (e.g. intermediates in cellml that are not used further)."
          << DihuContext::solverStructureVisualizer()->getDiagram();
      return false;
    }

    if (fromVectorNo == 0 && fromVectorIndex >= nFieldVariablesTerm1Vector1_)
    {
      if (!disableWarnings)
        LOG(WARNING) << "(1) Unconnected slot " << fromIndex << " of Term1 in " << settings_ << ": " << getDebugInformation()
          << ", from vector 0 (variable1), index " << fromVectorIndex
          << ", but this vector has only " << nFieldVariablesTerm1Vector1_ << " entries." << std::endl
          << "\nMaybe wrong numbers are given for \"connectedSlotsTerm1To2\" or \"connectedSlotsTerm2To1\" "
          << "or unneccessary slots have been created (e.g. intermediates in cellml that are not used further)."
          << DihuContext::solverStructureVisualizer()->getDiagram();
      return false;
    }
    else if (fromVectorNo == 1 && fromVectorIndex >= nFieldVariablesTerm1Vector2_)
    {
      if (!disableWarnings)
        LOG(WARNING) << "(2) Unconnected slot " << fromIndex << " of Term1 in " << settings_ << ": " << getDebugInformation()
          << ", from vector 1 (variable2), index " << fromVectorIndex
          << ", but this vector has only " << nFieldVariablesTerm1Vector2_ << " entries." << std::endl
          << "\nMaybe wrong numbers are given for \"connectedSlotsTerm1To2\" or \"connectedSlotsTerm2To1\" "
          << "or unneccessary slots have been created (e.g. intermediates in cellml that are not used further)."
          << DihuContext::solverStructureVisualizer()->getDiagram();
      return false;
    }

    int toIndex = connectorTerm1To2_[fromIndex].index;

    // if toIndex is set to None in the python config, this means we do not want to have a connection in this slot and it is no warning
    if (toIndex == -1)
      return false;

    toVectorNo = std::min(1,int(toIndex/nFieldVariablesTerm2Vector1_));
    toVectorIndex = toIndex - toVectorNo*nFieldVariablesTerm2Vector1_;

    avoidCopyIfPossible = connectorTerm1To2_[fromIndex].avoidCopyIfPossible;


    if (toVectorNo == 0 && toVectorIndex >= nFieldVariablesTerm2Vector1_)
    {
      if (!disableWarnings)
        LOG(WARNING) << "(3) Unconnected slot " << toIndex << " of Term2 in " << settings_ << ": " << getDebugInformation()
          << ", from vector " << fromVectorNo << " (variable" << fromVectorNo+1 << "), index " << fromVectorIndex
          << " to vector 0 (variable1), index " << toVectorIndex
          << ", but this vector has only " << nFieldVariablesTerm2Vector1_ << " entries." << std::endl
          << "\nMaybe wrong numbers are given for \"connectedSlotsTerm1To2\" or \"connectedSlotsTerm2To1\" "
          << "or unneccessary slots have been created (e.g. intermediates in cellml that are not used further)."
          << DihuContext::solverStructureVisualizer()->getDiagram();
      return false;
    }
    else if (toVectorNo == 1 && toVectorIndex >= nFieldVariablesTerm2Vector2_)
    {
      if (!disableWarnings)
        LOG(WARNING) << "(4) Unconnected slot " << toIndex << " of Term2 in " << settings_ << ": " << getDebugInformation()
          << ", from vector " << fromVectorNo << " (variable" << fromVectorNo+1 << "), index " << fromVectorIndex
          << " to vector 1 (variable2), index " << toVectorIndex
          << ", but this vector has only " << nFieldVariablesTerm2Vector2_ << " entries." << std::endl
          << "\nMaybe wrong numbers are given for \"connectedSlotsTerm1To2\" or \"connectedSlotsTerm2To1\" "
          << "or unneccessary slots have been created (e.g. intermediates in cellml that are not used further)."
          << DihuContext::solverStructureVisualizer()->getDiagram();
      return false;
    }
  }
  else
  {
    // if we map from term 2 to term 1

    // transform from index
    int fromIndex = fromVectorNo*nFieldVariablesTerm2Vector1_ + fromVectorIndex;

    LOG(DEBUG) << "fromVectorNo=" << fromVectorNo << ", fromVectorIndex=" << fromVectorIndex << ", fromIndex=" << fromIndex;

    // check for wrong indices and output errors
    if (fromIndex >= connectorTerm2To1_.size())
    {
      if (!disableWarnings)
        LOG(WARNING) << "(5) Unconnected slot " << fromIndex << " of Term2 in " << settings_ << ": " << getDebugInformation() << ", fromIndex (" << fromIndex << ") >= " << connectorTerm2To1_.size()
          << "\nThere are only " << connectorTerm2To1_.size() << " slots to connect from, but a connection from slot no " << fromIndex << " is needed."
          << "\nMaybe not enough entries are given for \"connectedSlotsTerm1To2\" or \"connectedSlotsTerm2To1\" "
          << "or unneccessary slots have been created (e.g. intermediates in cellml that are not used further)."
          << DihuContext::solverStructureVisualizer()->getDiagram();
      return false;
    }

    if (fromVectorNo == 0 && fromVectorIndex >= nFieldVariablesTerm2Vector1_)
    {
      if (!disableWarnings)
        LOG(WARNING) << "(6) Unconnected slot " << fromIndex << " of Term2 in " << settings_ << ": " << getDebugInformation() << ", from vector 0 (variable1), index " << fromVectorIndex
          << ", but this vector has only " << nFieldVariablesTerm2Vector1_ << " entries." << std::endl
          << "\nMaybe wrong numbers are given for \"connectedSlotsTerm1To2\" or \"connectedSlotsTerm2To1\" "
          << "or unneccessary slots have been created (e.g. intermediates in cellml that are not used further)."
          << DihuContext::solverStructureVisualizer()->getDiagram();
      return false;
    }
    else if (fromVectorNo == 1 && fromVectorIndex >= nFieldVariablesTerm2Vector2_)
    {
      if (!disableWarnings)
        LOG(WARNING) << "(7) Unconnected slot " << fromIndex << " of Term2 in " << settings_ << ": " << getDebugInformation() << ", from vector 1 (variable2), index " << fromVectorIndex
          << ", but this vector has only " << nFieldVariablesTerm2Vector2_ << " entries." << std::endl
          << "\nMaybe wrong numbers are given for \"connectedSlotsTerm1To2\" or \"connectedSlotsTerm2To1\" "
          << "or unneccessary slots have been created (e.g. intermediates in cellml that are not used further)."
          << DihuContext::solverStructureVisualizer()->getDiagram();
      return false;
    }

    int toIndex = connectorTerm2To1_[fromIndex].index;

    // if toIndex is set to None in the python config, this means we do not want to have a connection in this slot and it is no warning
    if (toIndex == -1)
      return false;

    toVectorNo = std::min(1,int(toIndex/nFieldVariablesTerm1Vector1_));
    toVectorIndex = toIndex - toVectorNo*nFieldVariablesTerm1Vector1_;

    LOG(DEBUG) << "toVectorNo=" << toVectorNo << ", toVectorIndex=" << toVectorIndex << ", toIndex=" << toIndex;

    avoidCopyIfPossible = connectorTerm2To1_[fromIndex].avoidCopyIfPossible;

    if (toVectorNo == 0 && toVectorIndex >= nFieldVariablesTerm1Vector1_)
    {
      if (!disableWarnings)
        LOG(WARNING) << "(8) Unconnected slot " << toIndex << " of Term1 in " << settings_ << ": " << getDebugInformation()
          << ", from vector " << fromVectorNo << " (variable" << fromVectorNo+1 << "), index " << fromVectorIndex
          << " to vector 0 (variable1), index " << toVectorIndex
          << ", but this vector has only " << nFieldVariablesTerm1Vector1_ << " entries." << std::endl
          << "\nMaybe wrong numbers are given for \"connectedSlotsTerm1To2\" or \"connectedSlotsTerm2To1\" "
          << "or unneccessary slots have been created (e.g. intermediates in cellml that are not used further)."
          << DihuContext::solverStructureVisualizer()->getDiagram();
      return false;
    }
    else if (toVectorNo == 1 && toVectorIndex >= nFieldVariablesTerm1Vector2_)
    {
      if (!disableWarnings)
        LOG(WARNING) << "(9) Unconnected slot " << toIndex << " of Term1 in " << settings_ << ": " << getDebugInformation()
          << ", from vector " << fromVectorNo << " (variable" << fromVectorNo+1 << "), index " << fromVectorIndex
          << " to vector 1 (variable2), index " << toVectorIndex
          << ", but this vector has only " << nFieldVariablesTerm1Vector2_ << " entries." << std::endl
          << "\nMaybe wrong numbers are given for \"connectedSlotsTerm1To2\" or \"connectedSlotsTerm2To1\" "
          << "or unneccessary slots have been created (e.g. intermediates in cellml that are not used further)."
          << DihuContext::solverStructureVisualizer()->getDiagram();
      return false;
    }
  }

  // completed successfully
  return true;
}
