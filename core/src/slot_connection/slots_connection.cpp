#include "slot_connection/slots_connection.h"

#include "control/diagnostic_tool/solver_structure_visualizer.h"

SlotsConnection::SlotsConnection(PythonConfig settings):
  fieldVariableNamesInitialized_(false), transferDirectionTerm1To2_(true), slotInformationInitialized_(false),
  settings_(settings), subSlotsConnection1_(nullptr), subSlotsConnection2_(nullptr), subSlotsConnection3_(nullptr), subSlotsConnection4_(nullptr)
{
  // parse values from settings
  std::vector<int> indicesTerm1To2;
  std::vector<int> indicesTerm2To1;

  settings.getOptionVector("connectedSlotsTerm1To2", indicesTerm1To2);
  settings.getOptionVector("connectedSlotsTerm2To1", indicesTerm2To1);

  // set indices in connectorForVisualizerTerm1To2_ and connectorForVisualizerTerm1To2_
  for (std::vector<int>::iterator iter = indicesTerm1To2.begin(); iter != indicesTerm1To2.end(); iter++)
  {
    Connector connector;
    connector.index = *iter;
    connector.avoidCopyIfPossible = false;
    connectorForVisualizerTerm1To2_.push_back(connector);
  }

  for (std::vector<int>::iterator iter = indicesTerm2To1.begin(); iter != indicesTerm2To1.end(); iter++)
  {
    Connector connector;
    connector.index = *iter;
    connector.avoidCopyIfPossible = false;
    connectorForVisualizerTerm2To1_.push_back(connector);
  }

  updateAvoidCopyIfPossible();
}

//! copy constructor
SlotsConnection::SlotsConnection(const SlotsConnection &rhs) :
  fieldVariableNamesInitialized_(false), transferDirectionTerm1To2_(true), slotInformationInitialized_(false)
{
  connectorForVisualizerTerm1To2_ = rhs.connectorForVisualizerTerm1To2_;
  connectorForVisualizerTerm2To1_ = rhs.connectorForVisualizerTerm2To1_;
  transferDirectionTerm1To2_ = rhs.transferDirectionTerm1To2_;
}

void SlotsConnection::addConnectionTerm1ToTerm2(int slotNoFrom, int slotNoTo)
{
  // check if connection is already present
  if (connectorForVisualizerTerm1To2_.size() > slotNoFrom)
  {
    if (connectorForVisualizerTerm1To2_[slotNoFrom].index == slotNoTo)
    {
      LOG(DEBUG) << "addConnectionTerm1ToTerm2(" << slotNoFrom << "," << slotNoTo << "), is already contained";

      // connection is already contained, return
      return;
    }
  }

  // resize number of entries in connectorForVisualizerTerm1To2_
  if (connectorForVisualizerTerm1To2_.size() <= slotNoFrom)
  {
    Connector unconnected;
    unconnected.index = -1;
    unconnected.avoidCopyIfPossible = false;
    connectorForVisualizerTerm1To2_.resize(slotNoFrom+1, unconnected);
  }

  // set new connection
  connectorForVisualizerTerm1To2_[slotNoFrom].index = slotNoTo;

  LOG(DEBUG) << "set new connection Term1->Term2 " << slotNoFrom << ":" << slotNoTo;

  // set the avoidCopyIfPossible value
  updateAvoidCopyIfPossible();

  // look-up table has to initialized again
  slotInformationInitialized_ = false;
  fieldVariableNamesInitialized_ = false;

  std::stringstream s;
  for (int i = 0; i < connectorForVisualizerTerm1To2_.size(); i++)
  {
    s << "  " << i << "->" << connectorForVisualizerTerm1To2_[i].index << " " << connectorForVisualizerTerm1To2_[i].avoidCopyIfPossible << "\n";
  }
  LOG(DEBUG) << s.str();
}

void SlotsConnection::addConnectionTerm2ToTerm1(int slotNoFrom, int slotNoTo)
{
  // check if connection is already present
  if (connectorForVisualizerTerm2To1_.size() > slotNoFrom)
  {
    if (connectorForVisualizerTerm2To1_[slotNoFrom].index == slotNoTo)
    {
      LOG(DEBUG) << "addConnectionTerm2ToTerm1(" << slotNoFrom << "," << slotNoTo << "), is already contained";

      // connection is already contained, return
      return;
    }
  }

  // resize number of entries in connectorForVisualizerTerm2To1_
  if (connectorForVisualizerTerm2To1_.size() <= slotNoFrom)
  {
    Connector unconnected;
    unconnected.index = -1;
    unconnected.avoidCopyIfPossible = false;
    connectorForVisualizerTerm2To1_.resize(slotNoFrom+1, unconnected);
  }

  // set new connection
  connectorForVisualizerTerm2To1_[slotNoFrom].index = slotNoTo;

  LOG(DEBUG) << "set new connection Term2->Term1 " << slotNoFrom << ":" << slotNoTo;

  // set the avoidCopyIfPossible value
  updateAvoidCopyIfPossible();

  // look-up table has to initialized again
  slotInformationInitialized_ = false;
  fieldVariableNamesInitialized_ = false;

  std::stringstream s;
  for (int i = 0; i < connectorForVisualizerTerm2To1_.size(); i++)
  {
    s << "  " << i << "->" << connectorForVisualizerTerm2To1_[i].index << " " << connectorForVisualizerTerm2To1_[i].avoidCopyIfPossible << "\n";
  }
  LOG(DEBUG) << s.str();
}

void SlotsConnection::setTransferDirection(bool term1To2)
{
  transferDirectionTerm1To2_ = term1To2;

  if (subSlotsConnection1_)
    subSlotsConnection1_->setTransferDirection(term1To2);
  if (subSlotsConnection2_)
    subSlotsConnection2_->setTransferDirection(term1To2);
  if (subSlotsConnection3_)
    subSlotsConnection3_->setTransferDirection(term1To2);
  if (subSlotsConnection4_)
    subSlotsConnection4_->setTransferDirection(term1To2);
}

void SlotsConnection::updateAvoidCopyIfPossible()
{
  // if field variable gets mapped in both directions, set avoidCopyIfPossible to true
  for (int i = 0; i < connectorForVisualizerTerm1To2_.size(); i++)
  {
    int mappedIndex = connectorForVisualizerTerm1To2_[i].index;

    // check if other direction is mapped the same way
    if (connectorForVisualizerTerm2To1_.size() > mappedIndex)
    {
      if (connectorForVisualizerTerm2To1_[mappedIndex].index == i)
      {
        connectorForVisualizerTerm1To2_[i].avoidCopyIfPossible = true;
        connectorForVisualizerTerm2To1_[i].avoidCopyIfPossible = true;
      }
    }
  }
}

//! get the connectors from term 1 to term 2
const std::vector<SlotsConnection::Connector> &SlotsConnection::connectorForVisualizerTerm1To2() const
{
  return connectorForVisualizerTerm1To2_;
}

//! get the connectors from term 2 to term 1
const std::vector<SlotsConnection::Connector> &SlotsConnection::connectorForVisualizerTerm2To1() const
{
  return connectorForVisualizerTerm2To1_;
}

std::string SlotsConnection::getDebugInformation() const
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
      if (vectorIndex >= 0 && vectorIndex < fieldVariableNamesTerm1Vector1_.size())
      {
        result << fieldVariableNamesTerm1Vector1_[vectorIndex];
      }
      else
      {
        result << "(out of range)";
      }
    }
    else
    {
      if (vectorIndex >= 0 && vectorIndex < fieldVariableNamesTerm1Vector2_.size())
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
        if (vectorIndex >= 0 && vectorIndex < fieldVariableNamesTerm2Vector1_.size())
        {
          result << fieldVariableNamesTerm2Vector1_[vectorIndex];
        }
        else
        {
          result << "(out of range)";
        }
      }
      else
      {
        if (vectorIndex >= 0 && vectorIndex < fieldVariableNamesTerm2Vector2_.size())
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
      if (vectorIndex >= 0 && vectorIndex < fieldVariableNamesTerm2Vector1_.size())
      {
        result << fieldVariableNamesTerm2Vector1_[vectorIndex];
      }
      else
      {
        result << "(out of range)";
      }
    }
    else
    {
      if (vectorIndex >= 0 && vectorIndex < fieldVariableNamesTerm2Vector2_.size())
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
        if (vectorIndex >= 0 && vectorIndex < fieldVariableNamesTerm1Vector1_.size())
        {
          result << fieldVariableNamesTerm1Vector1_[vectorIndex];
        }
        else
        {
          result << "(out of range)";
        }
      }
      else
      {
        if (vectorIndex >= 0 && vectorIndex < fieldVariableNamesTerm1Vector2_.size())
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

bool SlotsConnection::getSlotInformation(int fromVectorNo, int fromVectorIndex,
                                          int &toVectorNo, int &toVectorIndex, bool &avoidCopyIfPossible, bool disableWarnings) const
{

  // in release target, use precomputed look-up table for slot connection indices
#ifdef NDEBUG
  if (slotInformationInitialized_)
  {
    if (fromVectorIndex >= slotInformation_[transferDirectionTerm1To2_][fromVectorNo].size())
    {
      return false;
    }

    const Result &result = slotInformation_[transferDirectionTerm1To2_][fromVectorNo][fromVectorIndex];
    toVectorNo = result.toVectorNo;
    toVectorIndex = result.toVectorIndex;
    avoidCopyIfPossible = result.avoidCopyIfPossible;
    return result.successful;
  }
#endif

  disableWarnings = true;   // do not show warnings, they would also appear if SlotsConnectionDataType is a tuple, this is the case for MapDofs

  // fromVectorNo and toVectorNo are 0 or 1

#ifndef NDEBUG
  VLOG(1) << "getSlotInformation(" << fromVectorNo << "," << fromVectorIndex << ")" << getDebugInformation();
  LOG(DEBUG) << "getSlotInformation(" << fromVectorNo << "," << fromVectorIndex << "), " << (transferDirectionTerm1To2_? "1->2" : "2->1");
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
        LOG(DEBUG) << "(0) Unconnected slot " << fromIndex << " of Term1 in " << settings_ << ", fromIndex ("
          << fromIndex << "=" << fromVectorNo << "*" << nFieldVariablesTerm1Vector1_ << "+" << fromVectorIndex << ") >= " << connectorTerm1To2_.size()
          << "\nThere are only " << connectorTerm1To2_.size() << " slots to connect to, but a connection for slot no " << fromIndex << " is needed."
          << "\nMaybe not enough entries are given for \"connectedSlotsTerm1To2\" or \"connectedSlotsTerm2To1\" "
          << "or unneccessary slots have been created (e.g. algebraics in cellml that are not used further)."
          << getDebugInformation()
          << DihuContext::solverStructureVisualizer()->getDiagram();
      return false;
    }

    if (fromVectorNo == 0 && fromVectorIndex >= nFieldVariablesTerm1Vector1_)
    {
      if (!disableWarnings)
        LOG(WARNING) << "(1) Unconnected slot " << fromIndex << " of Term1 in " << settings_ 
          << ", from vector 0 (variable1), index " << fromVectorIndex
          << ", but this vector has only " << nFieldVariablesTerm1Vector1_ << " entries." << std::endl
          << "\nMaybe wrong numbers are given for \"connectedSlotsTerm1To2\" or \"connectedSlotsTerm2To1\" "
          << "or unneccessary slots have been created (e.g. algebraics in cellml that are not used further)."
          << getDebugInformation()
          << DihuContext::solverStructureVisualizer()->getDiagram();
      return false;
    }
    else if (fromVectorNo == 1 && fromVectorIndex >= nFieldVariablesTerm1Vector2_)
    {
      if (!disableWarnings)
        LOG(WARNING) << "(2) Unconnected slot " << fromIndex << " of Term1 in " << settings_ 
          << ", from vector 1 (variable2), index " << fromVectorIndex
          << ", but this vector has only " << nFieldVariablesTerm1Vector2_ << " entries." << std::endl
          << "\nMaybe wrong numbers are given for \"connectedSlotsTerm1To2\" or \"connectedSlotsTerm2To1\" "
          << "or unneccessary slots have been created (e.g. algebraics in cellml that are not used further)."
          << getDebugInformation()
          << DihuContext::solverStructureVisualizer()->getDiagram();
      return false;
    }

    int toIndex = connectorTerm1To2_[fromIndex].index;

    // if toIndex is set to None in the python config, this means we do not want to have a connection in this slot and it is no warning
    if (toIndex == -1)
    {
      return false;
    }

    if (nFieldVariablesTerm2Vector1_ == 0)
      toVectorNo = 1;
    else
      toVectorNo = std::min(1,int(toIndex/nFieldVariablesTerm2Vector1_));
    toVectorIndex = toIndex - toVectorNo*nFieldVariablesTerm2Vector1_;

    avoidCopyIfPossible = connectorTerm1To2_[fromIndex].avoidCopyIfPossible;


    if (toVectorNo == 0 && toVectorIndex >= nFieldVariablesTerm2Vector1_)
    {
      if (!disableWarnings)
        LOG(WARNING) << "(3) Unconnected slot " << toIndex << " of Term2 in " << settings_ 
          << ", from vector " << fromVectorNo << " (variable" << fromVectorNo+1 << "), index " << fromVectorIndex
          << " to vector 0 (variable1), index " << toVectorIndex
          << ", but this vector has only " << nFieldVariablesTerm2Vector1_ << " entries." << std::endl
          << "\nMaybe wrong numbers are given for \"connectedSlotsTerm1To2\" or \"connectedSlotsTerm2To1\" "
          << "or unneccessary slots have been created (e.g. algebraics in cellml that are not used further)."
          << getDebugInformation()
          << DihuContext::solverStructureVisualizer()->getDiagram();
        return false;
    }
    else if (toVectorNo == 1 && toVectorIndex >= nFieldVariablesTerm2Vector2_)
    {
      if (!disableWarnings)
        LOG(WARNING) << "(4) Unconnected slot " << toIndex << " of Term2 in " << settings_ 
          << ", from vector " << fromVectorNo << " (variable" << fromVectorNo+1 << "), index " << fromVectorIndex
          << " to vector 1 (variable2), index " << toVectorIndex
          << ", but this vector has only " << nFieldVariablesTerm2Vector2_ << " entries." << std::endl
          << "\nMaybe wrong numbers are given for \"connectedSlotsTerm1To2\" or \"connectedSlotsTerm2To1\" "
          << "or unneccessary slots have been created (e.g. algebraics in cellml that are not used further)."
          << getDebugInformation()
          << DihuContext::solverStructureVisualizer()->getDiagram();
      return false;
    }

    LOG(DEBUG) << "  yes, transfer fromVectorNo=" << fromVectorNo << ", fromVectorIndex=" << fromVectorIndex << "(fromIndex " << fromIndex << ")"
      << " to toVectorNo=" << toVectorNo << ", toVectorIndex=" << toVectorIndex << ", toIndex=" << toIndex;
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
        LOG(DEBUG) << "(5) Unconnected slot " << fromIndex << " of Term2 in " << settings_ << ", fromIndex (" << fromIndex << ") >= " << connectorTerm2To1_.size()
          << "\nThere are only " << connectorTerm2To1_.size() << " slots to connect from, but a connection from slot no " << fromIndex << " is needed."
          << "\nMaybe not enough entries are given for \"connectedSlotsTerm1To2\" or \"connectedSlotsTerm2To1\" "
          << "or unneccessary slots have been created (e.g. algebraics in cellml that are not used further)."
          << getDebugInformation()
          << DihuContext::solverStructureVisualizer()->getDiagram();
      return false;
    }

    if (fromVectorNo == 0 && fromVectorIndex >= nFieldVariablesTerm2Vector1_)
    {
      if (!disableWarnings)
        LOG(WARNING) << "(6) Unconnected slot " << fromIndex << " of Term2 in " << settings_ << ", from vector 0 (variable1), index " << fromVectorIndex
          << ", but this vector has only " << nFieldVariablesTerm2Vector1_ << " entries." << std::endl
          << "\nMaybe wrong numbers are given for \"connectedSlotsTerm1To2\" or \"connectedSlotsTerm2To1\" "
          << "or unneccessary slots have been created (e.g. algebraics in cellml that are not used further)."
          << getDebugInformation()
          << DihuContext::solverStructureVisualizer()->getDiagram();
      return false;
    }
    else if (fromVectorNo == 1 && fromVectorIndex >= nFieldVariablesTerm2Vector2_)
    {
      if (!disableWarnings)
        LOG(WARNING) << "(7) Unconnected slot " << fromIndex << " of Term2 in " << settings_ << ", from vector 1 (variable2), index " << fromVectorIndex
          << ", but this vector has only " << nFieldVariablesTerm2Vector2_ << " entries." << std::endl
          << "\nMaybe wrong numbers are given for \"connectedSlotsTerm1To2\" or \"connectedSlotsTerm2To1\" "
          << "or unneccessary slots have been created (e.g. algebraics in cellml that are not used further)."
          << getDebugInformation()
          << DihuContext::solverStructureVisualizer()->getDiagram();
      return false;
    }

    int toIndex = connectorTerm2To1_[fromIndex].index;

    // if toIndex is set to None in the python config, this means we do not want to have a connection in this slot and it is no warning
    if (toIndex == -1)
    {
      return false;
    }

    if (nFieldVariablesTerm1Vector1_ == 0)
      toVectorNo = 1;
    else
      toVectorNo = std::min(1,int(toIndex/nFieldVariablesTerm1Vector1_));
    toVectorIndex = toIndex - toVectorNo*nFieldVariablesTerm1Vector1_;

    LOG(DEBUG) << "toVectorNo=" << toVectorNo << ", toVectorIndex=" << toVectorIndex << ", toIndex=" << toIndex;

    avoidCopyIfPossible = connectorTerm2To1_[fromIndex].avoidCopyIfPossible;

    if (toVectorNo == 0 && toVectorIndex >= nFieldVariablesTerm1Vector1_)
    {
      if (!disableWarnings)
        LOG(WARNING) << "(8) Unconnected slot " << toIndex << " of Term1 in " << settings_
          << ", from vector " << fromVectorNo << " (variable" << fromVectorNo+1 << "), index " << fromVectorIndex
          << " to vector 0 (variable1), index " << toVectorIndex
          << ", but this vector has only " << nFieldVariablesTerm1Vector1_ << " entries." << std::endl
          << "\nMaybe wrong numbers are given for \"connectedSlotsTerm1To2\" or \"connectedSlotsTerm2To1\" "
          << "or unneccessary slots have been created (e.g. algebraics in cellml that are not used further)."
          << getDebugInformation()
          << DihuContext::solverStructureVisualizer()->getDiagram();
      return false;
    }
    else if (toVectorNo == 1 && toVectorIndex >= nFieldVariablesTerm1Vector2_)
    {
      if (!disableWarnings)
        LOG(WARNING) << "(9) Unconnected slot " << toIndex << " of Term1 in " << settings_
          << ", from vector " << fromVectorNo << " (variable" << fromVectorNo+1 << "), index " << fromVectorIndex
          << " to vector 1 (variable2), index " << toVectorIndex
          << ", but this vector has only " << nFieldVariablesTerm1Vector2_ << " entries." << std::endl
          << "\nMaybe wrong numbers are given for \"connectedSlotsTerm1To2\" or \"connectedSlotsTerm2To1\" "
          << "or unneccessary slots have been created (e.g. algebraics in cellml that are not used further)."
          << getDebugInformation()
          << DihuContext::solverStructureVisualizer()->getDiagram();
      return false;
    }

    LOG(DEBUG) << "  yes, transfer transfer fromVectorNo=" << fromVectorNo << ", fromVectorIndex=" << fromVectorIndex << "(fromIndex " << fromIndex << ") "
      << "to toVectorNo=" << toVectorNo << ", toVectorIndex=" << toVectorIndex << ", toIndex=" << toIndex;
  }

  // completed successfully
  return true;
}

std::shared_ptr<SlotsConnection> &SlotsConnection::subSlotsConnection1()
{
  return subSlotsConnection1_;
}

std::shared_ptr<SlotsConnection> &SlotsConnection::subSlotsConnection2()
{
  return subSlotsConnection2_;
}

std::shared_ptr<SlotsConnection> &SlotsConnection::subSlotsConnection3()
{
  return subSlotsConnection3_;
}

std::shared_ptr<SlotsConnection> &SlotsConnection::subSlotsConnection4()
{
  return subSlotsConnection4_;
}
