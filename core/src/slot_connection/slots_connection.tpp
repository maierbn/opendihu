#include "slot_connection/slots_connection.h"

template<typename FunctionSpaceType1, int nComponents1a, int nComponents1b, typename FunctionSpaceType2, int nComponents2a, int nComponents2b>
void SlotsConnection::
initialize(const Data::SlotConnectorData<FunctionSpaceType1,nComponents1a,nComponents1b> &transferableSolutionData1,
           const Data::SlotConnectorData<FunctionSpaceType2,nComponents2a,nComponents2b> &transferableSolutionData2,
           int offsetSlotNoData1, int offsetSlotNoData2)
{
  LOG(DEBUG) << "SlotsConnection::initialize, " << (transferDirectionTerm1To2_? "1->2" : "2->1") << ", offsets: " << offsetSlotNoData1 << "," << offsetSlotNoData2
    << ", initialized: " << fieldVariableNamesInitialized_;

  if (transferDirectionTerm1To2_)
  {
    nFieldVariablesTerm1Vector1_ = transferableSolutionData1.variable1.size();
    nFieldVariablesTerm1Vector2_ = transferableSolutionData1.variable2.size();
    nFieldVariablesTerm2Vector1_ = transferableSolutionData2.variable1.size();
    nFieldVariablesTerm2Vector2_ = transferableSolutionData2.variable2.size();

    offsetSlotNoData1_ = offsetSlotNoData1;
    offsetSlotNoData2_ = offsetSlotNoData2;
  }
  else
  {
    nFieldVariablesTerm1Vector1_ = transferableSolutionData2.variable1.size();
    nFieldVariablesTerm1Vector2_ = transferableSolutionData2.variable2.size();
    nFieldVariablesTerm2Vector1_ = transferableSolutionData1.variable1.size();
    nFieldVariablesTerm2Vector2_ = transferableSolutionData1.variable2.size();

    offsetSlotNoData1_ = offsetSlotNoData2;
    offsetSlotNoData2_ = offsetSlotNoData1;
  }

  if (!fieldVariableNamesInitialized_)
  {
    fieldVariableNamesInitialized_ = true;

    // clear all previous data
    fieldVariableNamesTerm1Vector1_.clear();
    fieldVariableNamesTerm1Vector2_.clear();
    fieldVariableNamesTerm2Vector1_.clear();
    fieldVariableNamesTerm2Vector2_.clear();
    connectorTerm1To2_.clear();
    connectorTerm2To1_.clear();
    for (int transferDirection = 0; transferDirection < 2; transferDirection++)
    {
      for (int fromVectorNo = 0; fromVectorNo < 2; fromVectorNo++)
      {
        slotInformation_[transferDirection][fromVectorNo].clear();
      }
    }

    // collect field variable names for debugging
    for (const Data::ComponentOfFieldVariable<FunctionSpaceType1,nComponents1a> &entry : transferableSolutionData1.variable1)
    {
      std::stringstream name;
      if (entry.values)
      {
        name << entry.values->name() << "." << entry.componentNo;
      }
      else
      {
        name << "(null)." << entry.componentNo;
      }
      fieldVariableNamesTerm1Vector1_.push_back(name.str());
    }

    for (const Data::ComponentOfFieldVariable<FunctionSpaceType1,nComponents1b> &entry : transferableSolutionData1.variable2)
    {
      std::stringstream name;
      if (entry.values)
      {
        name << entry.values->name() << "." << entry.componentNo;
      }
      else
      {
        name << "(null)." << entry.componentNo;
      }
      fieldVariableNamesTerm1Vector2_.push_back(name.str());
    }
    for (const Data::ComponentOfFieldVariable<FunctionSpaceType2,nComponents2a> &entry : transferableSolutionData2.variable1)
    {
      std::stringstream name;
      if (entry.values)
      {
        name << entry.values->name() << "." << entry.componentNo;
      }
      else
      {
        name << "(null)." << entry.componentNo;
      }
      fieldVariableNamesTerm2Vector1_.push_back(name.str());
    }

    for (const Data::ComponentOfFieldVariable<FunctionSpaceType2,nComponents2b> &entry : transferableSolutionData2.variable2)
    {
      std::stringstream name;
      if (entry.values)
      {
        name << entry.values->name() << "." << entry.componentNo;
      }
      else
      {
        name << "(null)." << entry.componentNo;
      }
      fieldVariableNamesTerm2Vector2_.push_back(name.str());
    }

    // consider offsets and compute connectorTerm1To2_ and connectorTerm2To1_
    // for example, offsetSlotNoData1_=2, offsetSlotNoData2_=1
    //
    // if connectorForVisualizerTerm1To2_ is:
    // slots data1  data2
    //       [ 0 +-->0 ] <-- handled in a previous SlotsConnection object
    //       [ 1-+   1
    //         2---->2
    //         3--+  3
    //         4  +->4
    // In this example, the SlotsConnection object consideres only (global) slots 2,3,4 which become 0,1,2 for data1 and for data2 1-4 which become 0-3.
    //
    // this leads to connectorTerm1To2_:
    // slots data1  data2
    //         0---+ 0
    //         1-+ +>1
    //         2 |   2
    //           +-->3

    // compute connectorTerm1To2_ from connectorForVisualizerTerm1To2_
    for (int globalSlotNo1 = offsetSlotNoData1_; globalSlotNo1 < connectorForVisualizerTerm1To2_.size(); globalSlotNo1++)
    {
      Connector connector = connectorForVisualizerTerm1To2_[globalSlotNo1];
      if (connector.index == -1)
      {
        connectorTerm1To2_.push_back(connector);
      }
      else if (connector.index >= offsetSlotNoData2_)
      {
        connector.index -= offsetSlotNoData2_;
        connectorTerm1To2_.push_back(connector);
      }
    }

    // compute connectorTerm2To1_ from connectorForVisualizerTerm2To1_
    for (int globalSlotNo2 = offsetSlotNoData2_; globalSlotNo2 < connectorForVisualizerTerm2To1_.size(); globalSlotNo2++)
    {
      Connector connector = connectorForVisualizerTerm2To1_[globalSlotNo2];
      if (connector.index == -1)
      {
        connectorTerm2To1_.push_back(connector);
      }
      else if (connector.index >= offsetSlotNoData1_)
      {
        connector.index -= offsetSlotNoData1_;
        connectorTerm2To1_.push_back(connector);
      }
    }

    LOG(DEBUG) << "offsetSlotNoData1_: " << offsetSlotNoData1_ << ", offsetSlotNoData2_: " << offsetSlotNoData2_;

    LOG(DEBUG) << "connectorTerm1To2_: ";
    int i = 0;
    for (std::vector<Connector>::const_iterator iter = connectorTerm1To2_.begin(); iter != connectorTerm1To2_.end(); iter++, i++)
    {
      LOG(DEBUG) << "  " << i << ". " << iter->index << " (avoidCopyIfPossible=" << iter->avoidCopyIfPossible << ")";
    }
    LOG(DEBUG) << "connectorTerm2To1_: ";
    i = 0;
    for (std::vector<Connector>::const_iterator iter = connectorTerm2To1_.begin(); iter != connectorTerm2To1_.end(); iter++, i++)
    {
      LOG(DEBUG) << "  " << i << ". " << iter->index << " (avoidCopyIfPossible=" << iter->avoidCopyIfPossible << ")";
    }

    LOG(DEBUG) << "connectorForVisualizerTerm1To2_: ";
    i = 0;
    for (std::vector<Connector>::const_iterator iter = connectorForVisualizerTerm1To2_.begin(); iter != connectorForVisualizerTerm1To2_.end(); iter++, i++)
    {
      LOG(DEBUG) << "  " << i << ". " << iter->index << " (avoidCopyIfPossible=" << iter->avoidCopyIfPossible << ")";
    }
    LOG(DEBUG) << "connectorForVisualizerTerm2To1_: ";
    i = 0;
    for (std::vector<Connector>::const_iterator iter = connectorForVisualizerTerm2To1_.begin(); iter != connectorForVisualizerTerm2To1_.end(); iter++, i++)
    {
      LOG(DEBUG) << "  " << i << ". " << iter->index << " (avoidCopyIfPossible=" << iter->avoidCopyIfPossible << ")";
    }

    initializeSlotInformation(transferableSolutionData1, transferableSolutionData2);
  }
}

template<typename FunctionSpaceType1, int nComponents1a, int nComponents1b, typename FunctionSpaceType2, int nComponents2a, int nComponents2b>
void SlotsConnection::
initializeSlotInformation(const Data::SlotConnectorData<FunctionSpaceType1,nComponents1a,nComponents1b> &transferableSolutionData1,
                          const Data::SlotConnectorData<FunctionSpaceType2,nComponents2a,nComponents2b> &transferableSolutionData2)
{
  if (slotInformationInitialized_)
    return;

  // save previous value of the variable transferDirectionTerm1To2_
  bool previousTransferDirectionTerm1To2 = transferDirectionTerm1To2_;

  // variable1 from 1 to 2
  // ----------------------
  transferDirectionTerm1To2_ = true;

  // determine if for this variable there are copy and non-copy entries. If so, non-copy is not possible and all have to be set to copy.
  bool copyRequired = false;

  // iterate over the first vector of variables
  for (int i = 0; i < transferableSolutionData1.variable1.size(); i++)
  {
    VLOG(1) << "i=" << i;

    // call getSlotInformation, which outputs the warnings
    int fromVectorNo = 0;
    int fromVectorIndex = i;
    int toVectorNo = 0;
    int toVectorIndex = 0;
    bool avoidCopyIfPossible = true;
    bool slotIsConnected = getSlotInformation(fromVectorNo, fromVectorIndex, toVectorNo, toVectorIndex, avoidCopyIfPossible);

    if (!slotIsConnected)
      continue;

    if (!avoidCopyIfPossible)
    {
      copyRequired = true;
    }

    // check if other entries of the variable map to a different variable
    if (avoidCopyIfPossible)
    {
      // loop over other connections of this variable
      for (int i2 = i+1; i2 < transferableSolutionData1.variable1.size(); i2++)
      {
        VLOG(1) << "i2=" << i2;

        int fromVectorNo = 0;
        int fromVectorIndex = i2;
        int toVectorNo2 = 0;
        int toVectorIndex2 = 0;
        bool avoidCopyIfPossible = true;
        bool slotIsConnected = getSlotInformation(fromVectorNo, fromVectorIndex, toVectorNo2, toVectorIndex2, avoidCopyIfPossible);

        if (!slotIsConnected)
          continue;

        if (toVectorNo2 != toVectorNo || !avoidCopyIfPossible)
        {
          LOG(DEBUG) << "slot " << i << " which is on variable1 maps to variable " << toVectorNo+1 << ", but slot "
            << i2 << ", which is also on variable1 maps to variable " << toVectorNo2+1 << ". avoidCopyIfPossible = " << avoidCopyIfPossible;
          LOG(DEBUG) << "Therefore copyRequired.";
          copyRequired = true;
        }
      }
    }
  }

  if (copyRequired)
  {
    LOG(DEBUG) << "transferDirectionTerm1To2_=" << transferDirectionTerm1To2_
      << ", copyRequired was set, set all connections for variable1 to copy";

    // set in all entries of the first field variable avoidCopyIfPossible=false
    for (int i = 0; i < transferableSolutionData1.variable1.size(); i++)
    {
      VLOG(1) << " i=" << i;

      int fromVectorNo = 0;
      int fromVectorIndex = i;
      int toVectorNo = 0;
      int toVectorIndex = 0;
      bool avoidCopyIfPossible = true;
      bool slotIsConnected = getSlotInformation(fromVectorNo, fromVectorIndex, toVectorNo, toVectorIndex, avoidCopyIfPossible);

      if (!slotIsConnected)
        continue;

      int fromIndex = fromVectorNo*nFieldVariablesTerm1Vector1_ + fromVectorIndex;
      connectorTerm1To2_[fromIndex].avoidCopyIfPossible = false;
      if (!connectorForVisualizerTerm1To2_.empty())
        connectorForVisualizerTerm1To2_[fromIndex+offsetSlotNoData1_].avoidCopyIfPossible = false;
    }
  }

  // variable2 from 1 to 2
  // ----------------------
  transferDirectionTerm1To2_ = true;

  // determine if for this variable there are copy and non-copy entries. If so, non-copy is not possible and all have to be set to copy.
  copyRequired = false;

  // iterate over the first vector of variables
  for (int i = 0; i < transferableSolutionData1.variable2.size(); i++)
  {
    VLOG(1) << "i=" << i;

    // call getSlotInformation, which outputs the warnings
    int fromVectorNo = 1;
    int fromVectorIndex = i;
    int toVectorNo = 0;
    int toVectorIndex = 0;
    bool avoidCopyIfPossible = true;
    bool slotIsConnected = getSlotInformation(fromVectorNo, fromVectorIndex, toVectorNo, toVectorIndex, avoidCopyIfPossible);

    if (!slotIsConnected)
      continue;

    if (!avoidCopyIfPossible)
    {
      copyRequired = true;
    }

    // check if other entries of the variable map to a different variable
    if (avoidCopyIfPossible)
    {
      // loop over other connections of this variable
      for (int i2 = i+1; i2 < transferableSolutionData1.variable2.size(); i2++)
      {
        VLOG(1) << "i2=" << i2;

        int fromVectorNo = 1;
        int fromVectorIndex = i2;
        int toVectorNo2 = 0;
        int toVectorIndex2 = 0;
        bool avoidCopyIfPossible = true;
        bool slotIsConnected = getSlotInformation(fromVectorNo, fromVectorIndex, toVectorNo2, toVectorIndex2, avoidCopyIfPossible);

        if (!slotIsConnected)
          continue;

        if (toVectorNo2 != toVectorNo || !avoidCopyIfPossible)
        {
          LOG(DEBUG) << "slot " << i << " which is on variable2 maps to variable " << toVectorNo+1 << ", but slot "
            << i2 << ", which is also on variable2 maps to variable " << toVectorNo2+1 << ". avoidCopyIfPossible = " << avoidCopyIfPossible;
          LOG(DEBUG) << "Therefore copyRequired.";
          copyRequired = true;
        }
      }
    }
  }

  if (copyRequired)
  {
    LOG(DEBUG) << "transferDirectionTerm1To2_=" << transferDirectionTerm1To2_
      << ", copyRequired was set, set all connections for variable2 to copy";

    // set in all entries of the first field variable avoidCopyIfPossible=false
    for (int i = 0; i < transferableSolutionData1.variable2.size(); i++)
    {
      VLOG(1) << " i=" << i;

      int fromVectorNo = 1;
      int fromVectorIndex = i;
      int toVectorNo = 0;
      int toVectorIndex = 0;
      bool avoidCopyIfPossible = true;
      bool slotIsConnected = getSlotInformation(fromVectorNo, fromVectorIndex, toVectorNo, toVectorIndex, avoidCopyIfPossible);

      if (!slotIsConnected)
        continue;

      int fromIndex = fromVectorNo*nFieldVariablesTerm1Vector1_ + fromVectorIndex;
      connectorTerm1To2_[fromIndex].avoidCopyIfPossible = false;
      if (!connectorForVisualizerTerm1To2_.empty())
        connectorForVisualizerTerm1To2_[fromIndex+offsetSlotNoData1_].avoidCopyIfPossible = false;
    }
  }


  // variable1 from 2 to 1
  // ----------------------
  transferDirectionTerm1To2_ = false;

  // determine if for this variable there are copy and non-copy entries. If so, non-copy is not possible and all have to be set to copy.
  copyRequired = false;

  // iterate over the first vector of variables
  for (int i = 0; i < transferableSolutionData2.variable1.size(); i++)
  {
    VLOG(1) << "i=" << i;

    // call getSlotInformation, which outputs the warnings
    int fromVectorNo = 0;
    int fromVectorIndex = i;
    int toVectorNo = 0;
    int toVectorIndex = 0;
    bool avoidCopyIfPossible = true;
    bool slotIsConnected = getSlotInformation(fromVectorNo, fromVectorIndex, toVectorNo, toVectorIndex, avoidCopyIfPossible);

    if (!slotIsConnected)
      continue;

    if (!avoidCopyIfPossible)
    {
      copyRequired = true;
    }

    // check if other entries of the variable map to a different variable
    if (avoidCopyIfPossible)
    {
      // loop over other connections of this variable
      for (int i2 = i+1; i2 < transferableSolutionData2.variable1.size(); i2++)
      {
        VLOG(1) << "i2=" << i2;

        int fromVectorNo = 0;
        int fromVectorIndex = i2;
        int toVectorNo2 = 0;
        int toVectorIndex2 = 0;
        bool avoidCopyIfPossible = true;
        bool slotIsConnected = getSlotInformation(fromVectorNo, fromVectorIndex, toVectorNo2, toVectorIndex2, avoidCopyIfPossible);

        if (!slotIsConnected)
          continue;

        if (toVectorNo2 != toVectorNo || !avoidCopyIfPossible)
        {
          LOG(DEBUG) << "slot " << i << " which is on variable1 maps to variable " << toVectorNo+1 << ", but slot "
            << i2 << ", which is also on variable1 maps to variable " << toVectorNo2+1 << ". avoidCopyIfPossible = " << avoidCopyIfPossible;
          LOG(DEBUG) << "Therefore copyRequired.";
          copyRequired = true;
        }
      }
    }
  }

  if (copyRequired)
  {
    LOG(DEBUG) << "transferDirectionTerm1To2_=" << transferDirectionTerm1To2_
      << ", copyRequired was set, set all connections for variable1 to copy";

    // set in all entries of the first field variable avoidCopyIfPossible=false
    for (int i = 0; i < transferableSolutionData2.variable1.size(); i++)
    {
      VLOG(1) << " i=" << i;

      int fromVectorNo = 0;
      int fromVectorIndex = i;
      int toVectorNo = 0;
      int toVectorIndex = 0;
      bool avoidCopyIfPossible = true;
      bool slotIsConnected = getSlotInformation(fromVectorNo, fromVectorIndex, toVectorNo, toVectorIndex, avoidCopyIfPossible);

      if (!slotIsConnected)
        continue;

      int fromIndex = fromVectorNo*nFieldVariablesTerm2Vector1_ + fromVectorIndex;
      connectorTerm2To1_[fromIndex].avoidCopyIfPossible = false;
      if (!connectorForVisualizerTerm2To1_.empty())
        connectorForVisualizerTerm2To1_[fromIndex+offsetSlotNoData2_].avoidCopyIfPossible = false;
    }
  }


  // variable2 from 2 to 1
  // ----------------------
  transferDirectionTerm1To2_ = false;

  // determine if for this variable there are copy and non-copy entries. If so, non-copy is not possible and all have to be set to copy.
  copyRequired = false;

  // iterate over the first vector of variables
  for (int i = 0; i < transferableSolutionData2.variable2.size(); i++)
  {
    VLOG(1) << "i=" << i;

    // call getSlotInformation, which outputs the warnings
    int fromVectorNo = 1;
    int fromVectorIndex = i;
    int toVectorNo = 0;
    int toVectorIndex = 0;
    bool avoidCopyIfPossible = true;
    bool slotIsConnected = getSlotInformation(fromVectorNo, fromVectorIndex, toVectorNo, toVectorIndex, avoidCopyIfPossible);

    if (!slotIsConnected)
      continue;

    if (!avoidCopyIfPossible)
    {
      copyRequired = true;
    }

    // check if other entries of the variable map to a different variable
    if (avoidCopyIfPossible)
    {
      // loop over other connections of this variable
      for (int i2 = i+1; i2 < transferableSolutionData2.variable2.size(); i2++)
      {
        VLOG(1) << "i2=" << i2;

        int fromVectorNo = 1;
        int fromVectorIndex = i2;
        int toVectorNo2 = 0;
        int toVectorIndex2 = 0;
        bool avoidCopyIfPossible = true;
        bool slotIsConnected = getSlotInformation(fromVectorNo, fromVectorIndex, toVectorNo2, toVectorIndex2, avoidCopyIfPossible);

        if (!slotIsConnected)
          continue;

        if (toVectorNo2 != toVectorNo || !avoidCopyIfPossible)
        {
          LOG(DEBUG) << "slot " << i << " which is on variable2 maps to variable " << toVectorNo+1 << ", but slot "
            << i2 << ", which is also on variable2 maps to variable " << toVectorNo2+1 << ". avoidCopyIfPossible = " << avoidCopyIfPossible;
          LOG(DEBUG) << "Therefore copyRequired.";
          copyRequired = true;
        }
      }
    }
  }

  if (copyRequired)
  {
    LOG(DEBUG) << "transferDirectionTerm1To2_=" << transferDirectionTerm1To2_
      << ", copyRequired was set, set all connections for variable1 to copy";

    // set in all entries of the first field variable avoidCopyIfPossible=false
    for (int i = 0; i < transferableSolutionData2.variable2.size(); i++)
    {
      VLOG(1) << " i=" << i;

      int fromVectorNo = 1;
      int fromVectorIndex = i;
      int toVectorNo = 0;
      int toVectorIndex = 0;
      bool avoidCopyIfPossible = true;
      bool slotIsConnected = getSlotInformation(fromVectorNo, fromVectorIndex, toVectorNo, toVectorIndex, avoidCopyIfPossible);

      if (!slotIsConnected)
        continue;

      int fromIndex = fromVectorNo*nFieldVariablesTerm2Vector1_ + fromVectorIndex;
      connectorTerm2To1_[fromIndex].avoidCopyIfPossible = false;
      if (!connectorForVisualizerTerm2To1_.empty())
        connectorForVisualizerTerm2To1_[fromIndex+offsetSlotNoData2_].avoidCopyIfPossible = false;
    }
  }

  // fill look-up table slotInformation_[transferDirectionTerm1To2][fromVectorNo][fromVectorIndex]
  for (int transferDirectionTerm1To2 = 0; transferDirectionTerm1To2 < 2; transferDirectionTerm1To2++)
  {
    if (transferDirectionTerm1To2 == 0)
      LOG(DEBUG) << "initialize lookup-table for terms 2->1";
    else
      LOG(DEBUG) << "initialize lookup-table for terms 1->2";

    for (int fromVectorNo = 0; fromVectorNo < 2; fromVectorNo++)
    {
      // iterate to the next fromVectorIndex until there were 3 unconnected slots, then it is assumed that there are no more
      int nUnsuccessful = 0;
      for (int fromVectorIndex = 0; ; fromVectorIndex++)
      {
        transferDirectionTerm1To2_ = transferDirectionTerm1To2;
        Result result;
        result.successful = getSlotInformation(fromVectorNo, fromVectorIndex, result.toVectorNo, result.toVectorIndex, result.avoidCopyIfPossible, true);

        slotInformation_[transferDirectionTerm1To2][fromVectorNo].push_back(result);

        if (!result.successful)
          nUnsuccessful++;

        if (nUnsuccessful == 3)
          break;
      }
    }
  }

#ifndef NDEBUG
  LOG(DEBUG) << "lookup-table initialization done:";
  for (int transferDirectionTerm1To2 = 0; transferDirectionTerm1To2 < 2; transferDirectionTerm1To2++)
  {
    for (int fromVectorNo = 0; fromVectorNo < 2; fromVectorNo++)
    {
      for (int i = 0; i < slotInformation_[transferDirectionTerm1To2][fromVectorNo].size(); i++)
      {
        LOG(DEBUG) << "  slotInformation_[" << transferDirectionTerm1To2 << "][" << fromVectorNo << "][" << i << "] = " << slotInformation_[transferDirectionTerm1To2][fromVectorNo][i].successful;
      }
    }
  }
#endif

  // restore previous value of the variable transferDirectionTerm1To2_
  transferDirectionTerm1To2_ = previousTransferDirectionTerm1To2;
  slotInformationInitialized_ = true;
}
