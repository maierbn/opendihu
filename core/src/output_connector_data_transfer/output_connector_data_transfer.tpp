#include "output_connector_data_transfer/output_connector_data_transfer.h"

#include <vector>
#include <tuple>
#include "easylogging++.h"

/** Transfer between two field variables with given component number
 */
template<typename FunctionSpaceType1, int nComponents1a, int nComponents1b, typename FunctionSpaceType2, int nComponents2a, int nComponents2b>
void SolutionVectorMapping<
  Data::OutputConnectorData<FunctionSpaceType1,nComponents1a,nComponents1b>,
  Data::OutputConnectorData<FunctionSpaceType2,nComponents2a,nComponents2b>
>::transfer(const std::shared_ptr<Data::OutputConnectorData<FunctionSpaceType1,nComponents1a,nComponents1b>> transferableSolutionData1,
            std::shared_ptr<Data::OutputConnectorData<FunctionSpaceType2,nComponents2a,nComponents2b>> transferableSolutionData2,
            OutputConnection &outputConnection)
{
  LOG(DEBUG) << "transfer standard";

  // initialize output connection object
  outputConnection.initialize(*transferableSolutionData1, *transferableSolutionData2);

  typedef FieldVariable::FieldVariable<FunctionSpaceType1,nComponents1a> FieldVariable1;

  // for the first vector of variables (the "states" in case of CellMLAdapter)
  for (int i = 0; i < transferableSolutionData1->variable1.size(); i++)
  {
    int fromVectorNo = 0;
    int fromVectorIndex = i;
    int toVectorNo = 0;
    int toVectorIndex = 0;
    bool avoidCopyIfPossible = true;
    bool slotIsConnected = outputConnection.getSlotInformation(fromVectorNo, fromVectorIndex, toVectorNo, toVectorIndex, avoidCopyIfPossible);

    if (!slotIsConnected)
      continue;

    std::shared_ptr<FieldVariable1> fieldVariable1 = transferableSolutionData1->variable1[fromVectorIndex].values;
    int componentNo1 = transferableSolutionData1->variable1[fromVectorIndex].componentNo;

    LOG(DEBUG) << "map slot from variable1, index " << fromVectorIndex
      << " to variable" << toVectorNo+1 << ", index " << toVectorIndex;

    if (componentNo1 < 0)
    {
      LOG(DEBUG) << "do not map this slot";
      continue;
    }

    if (!fieldVariable1)
    {
      LOG(FATAL) << "FieldVariable1 is null!";
    }

    if (toVectorNo == 0)
    {
      typedef FieldVariable::FieldVariable<FunctionSpaceType2,nComponents2a> FieldVariable2;
      std::shared_ptr<FieldVariable2> &fieldVariable2 = transferableSolutionData2->variable1[toVectorIndex].values;
      int componentNo2 = transferableSolutionData2->variable1[toVectorIndex].componentNo;
      assert(fieldVariable2);

      LOG(DEBUG) << "  " << fieldVariable1->name() << "[" << componentNo1 << "] -> "
        << fieldVariable2->name() << "[" << componentNo2 << "], avoidCopyIfPossible: " << avoidCopyIfPossible;

      // perform the mapping
      DihuContext::meshManager()->template prepareMapping<FieldVariable1,FieldVariable2>(fieldVariable1, fieldVariable2);
      DihuContext::meshManager()->template map<FieldVariable1,FieldVariable2>(fieldVariable1, componentNo1, fieldVariable2, componentNo2, avoidCopyIfPossible);
      DihuContext::meshManager()->template finalizeMapping<FieldVariable1,FieldVariable2>(fieldVariable1, fieldVariable2);
    }
    else
    {
      typedef FieldVariable::FieldVariable<FunctionSpaceType2,nComponents2b> FieldVariable2;
      std::shared_ptr<FieldVariable2> &fieldVariable2 = transferableSolutionData2->variable2[toVectorIndex].values;
      int componentNo2 = transferableSolutionData2->variable2[toVectorIndex].componentNo;
      assert(fieldVariable2);


      LOG(DEBUG) << "  " << fieldVariable1->name() << "[" << componentNo1 << "] -> "
        << fieldVariable2->name() << "[" << componentNo2 << "], avoidCopyIfPossible: " << avoidCopyIfPossible;

      // perform the mapping
      DihuContext::meshManager()->template prepareMapping<FieldVariable1,FieldVariable2>(fieldVariable1, fieldVariable2);
      DihuContext::meshManager()->template map<FieldVariable1,FieldVariable2>(fieldVariable1, componentNo1, fieldVariable2, componentNo2, avoidCopyIfPossible);
      DihuContext::meshManager()->template finalizeMapping<FieldVariable1,FieldVariable2>(fieldVariable1, fieldVariable2);
    }
  }

  // for the second vector of variables (the "intermediates" in case of CellMLAdapter)
  for (int i = 0; i < transferableSolutionData1->variable2.size(); i++)
  {
    int fromVectorNo = 1;
    int fromVectorIndex = i;
    int toVectorNo = 0;
    int toVectorIndex = 0;
    bool avoidCopyIfPossible = true;
    bool slotIsConnected = outputConnection.getSlotInformation(fromVectorNo, fromVectorIndex, toVectorNo, toVectorIndex, avoidCopyIfPossible);

    if (!slotIsConnected)
      continue;

    typedef FieldVariable::FieldVariable<FunctionSpaceType1,nComponents1b> FieldVariable1;

    std::shared_ptr<FieldVariable1> fieldVariable1 = transferableSolutionData1->variable2[fromVectorIndex].values;
    int componentNo1 = transferableSolutionData1->variable2[fromVectorIndex].componentNo;

    LOG(DEBUG) << "map slot from variable2, index " << fromVectorIndex
      << " to variable" << toVectorNo+1 << ", index " << toVectorIndex;

    if (componentNo1 < 0)
    {
      LOG(DEBUG) << "do not map this slot";
      continue;
    }

    if (!fieldVariable1)
    {
      LOG(FATAL) << "FieldVariable1 is null!";
    }

    if (toVectorNo == 0)
    {
      typedef FieldVariable::FieldVariable<FunctionSpaceType2,nComponents2a> FieldVariable2;
      std::shared_ptr<FieldVariable2> &fieldVariable2 = transferableSolutionData2->variable1[toVectorIndex].values;
      int componentNo2 = transferableSolutionData2->variable1[toVectorIndex].componentNo;
      assert(fieldVariable2);

      LOG(DEBUG) << "  " << fieldVariable1->name() << "[" << componentNo1 << "] -> "
        << fieldVariable2->name() << "[" << componentNo2 << "], avoidCopyIfPossible: " << avoidCopyIfPossible;

      // perform the mapping
      DihuContext::meshManager()->template prepareMapping<FieldVariable1,FieldVariable2>(fieldVariable1, fieldVariable2);
      DihuContext::meshManager()->template map<FieldVariable1,FieldVariable2>(fieldVariable1, componentNo1, fieldVariable2, componentNo2, avoidCopyIfPossible);
      DihuContext::meshManager()->template finalizeMapping<FieldVariable1,FieldVariable2>(fieldVariable1, fieldVariable2);
    }
    else
    {
      typedef FieldVariable::FieldVariable<FunctionSpaceType2,nComponents2b> FieldVariable2;
      std::shared_ptr<FieldVariable2> &fieldVariable2 = transferableSolutionData2->variable2[toVectorIndex].values;
      int componentNo2 = transferableSolutionData2->variable2[toVectorIndex].componentNo;
      assert(fieldVariable2);

      LOG(DEBUG) << "  " << fieldVariable1->name() << "[" << componentNo1 << "] -> "
        << fieldVariable2->name() << "[" << componentNo2 << "], avoidCopyIfPossible: " << avoidCopyIfPossible;

      // perform the mapping
      DihuContext::meshManager()->template prepareMapping<FieldVariable1,FieldVariable2>(fieldVariable1, fieldVariable2);
      DihuContext::meshManager()->template map<FieldVariable1,FieldVariable2>(fieldVariable1, componentNo1, fieldVariable2, componentNo2, avoidCopyIfPossible);
      DihuContext::meshManager()->template finalizeMapping<FieldVariable1,FieldVariable2>(fieldVariable1, fieldVariable2);
    }
  }

  // transfer geometry field if it was set in transferableSolutionData1
  if (transferableSolutionData1->geometryField && !transferableSolutionData2->variable1.empty())
  {
    LOG(DEBUG) << "transfer geometry field";

    // get source field variable, this is the same for all fibers
    typedef FieldVariable::FieldVariable<FunctionSpaceType1,3> FieldVariableSource;
    typedef FieldVariable::FieldVariable<FunctionSpaceType2,3> FieldVariableTarget;

    std::shared_ptr<FieldVariableSource> geometryFieldSource = transferableSolutionData1->geometryField;
    std::shared_ptr<FieldVariableTarget> geometryFieldTarget = std::make_shared<FieldVariableTarget>(transferableSolutionData2->variable1[0].values->functionSpace()->geometryField());


    // perform the mapping
    DihuContext::meshManager()->template prepareMapping<FieldVariableSource,FieldVariableTarget>(geometryFieldSource, geometryFieldTarget);

    // map the whole geometry field (all components), do not avoid copy
    DihuContext::meshManager()->template map<FieldVariableSource,FieldVariableTarget>(geometryFieldSource, -1, geometryFieldTarget, -1, false);
    DihuContext::meshManager()->template finalizeMapping<FieldVariableSource,FieldVariableTarget>(geometryFieldSource, geometryFieldTarget);
  }

  VLOG(1) << "at the end of output_connector_data_transfer_cellml.";
  VLOG(1) << "transferableSolutionData1: " << *transferableSolutionData1;
  VLOG(1) << "transferableSolutionData2: " << *transferableSolutionData2;
}
