#include "control/map_dofs/map_dofs.h"

namespace Control
{

//! advance simulation by the given time span
template<typename FunctionSpaceType, typename NestedSolverType>
void MapDofs<FunctionSpaceType,NestedSolverType>::
advanceTimeSpan()
{
  LOG(DEBUG) << "MapDofs::advanceTimeSpan, " << mappingsBeforeComputation_.size() << " mappings before computation, "
    << mappingsAfterComputation_.size() << " mappings after computation.";
  LOG(DEBUG) << "MapDofs::performMappings beforeComputation";

  // perform mapping from settings "beforeComputation"
  performMappings(mappingsBeforeComputation_);

  // compute the simulation in the current time span with the nested solver
  nestedSolver_.advanceTimeSpan();

  LOG(DEBUG) << "MapDofs::performMappings afterComputation";
  // perform mapping from settings "afterComputation"
  performMappings(mappingsAfterComputation_);
}

//! run solution process, this first calls initialize() and then run()
template<typename FunctionSpaceType, typename NestedSolverType>
void MapDofs<FunctionSpaceType,NestedSolverType>::
run()
{
  // initialize the solver
  initialize();

  // advance one timestep
  advanceTimeSpan();
}

//! reset state of this object, such that a new initialize() is necessary ("uninitialize")
template<typename FunctionSpaceType, typename NestedSolverType>
void MapDofs<FunctionSpaceType,NestedSolverType>::
reset()
{
  nestedSolver_.reset();
}

template<typename FunctionSpaceType, typename NestedSolverType>
void MapDofs<FunctionSpaceType,NestedSolverType>::
setTimeSpan(double startTime, double endTime)
{
  nestedSolver_.setTimeSpan(startTime, endTime);
}

template<typename FunctionSpaceType, typename NestedSolverType>
void MapDofs<FunctionSpaceType,NestedSolverType>::
performMappings(std::vector<DofsMappingType> &mappings)
{
  int ownRankNo = data_.functionSpace()->meshPartition()->rankSubset()->ownRankNo();

  for (DofsMappingType &mapping : mappings)
  {
    LOG(DEBUG) << "MapDofs::perform mapping slots " << mapping.outputConnectorSlotNoFrom << " -> " << mapping.outputConnectorSlotNoTo << ", "
      << (mapping.mode == DofsMappingType::modeCopyLocal? " modeCopyLocal": (mapping.mode == DofsMappingType::modeCopyLocalIfPositive? "modeCopyLocalIfPositive": " modeCommunicate"));

    static std::map<int,std::vector<double>> valuesToSendToRanks;
    static std::vector<double> receivedValues;

    receivedValues.clear();

    // get values that will be transferred
    // loop over the ranks and the dofs that have to be send to them
    for (std::pair<int,std::vector<dof_no_t>> rankNoAndDofNosLocal : mapping.dofNosLocalOfValuesToSendToRanks)
    {
      // get rank and dof nos
      int rankNo = rankNoAndDofNosLocal.first;
      std::vector<dof_no_t> &dofNosLocal = rankNoAndDofNosLocal.second;

      // get values at those dof nos
      valuesToSendToRanks[rankNo].clear();
      slotGetValues(mapping.outputConnectorSlotNoFrom, mapping.outputConnectorArrayIndexFrom, dofNosLocal, valuesToSendToRanks[rankNo]);
    }

    LOG(DEBUG) << "get values " << valuesToSendToRanks;

    // if communication is involved
    if (mapping.mode == DofsMappingType::modeCommunicate)
    {
      mapping.valueCommunicator.communicate(valuesToSendToRanks, receivedValues);

      // set values
      // set receivedValues at dofs receivedValueDofNosLocal
      slotSetValues(mapping.outputConnectorSlotNoTo, mapping.outputConnectorArrayIndexTo, mapping.receivedValueDofNosLocal, receivedValues, INSERT_VALUES);
    }
    else  if (mapping.mode == DofsMappingType::modeCopyLocal)
    {
      // collect all values
      for (double value : valuesToSendToRanks[ownRankNo])
      {
        receivedValues.push_back(value);
      }

      // set values
      // set receivedValues at dofs receivedValueDofNosLocal
      slotSetValues(mapping.outputConnectorSlotNoTo, mapping.outputConnectorArrayIndexTo, mapping.receivedValueDofNosLocal, receivedValues, INSERT_VALUES);
    }
    else if (mapping.mode == DofsMappingType::modeCopyLocalIfPositive)
    {
      // store dof nos for positive values in maskedDofNosLocal
      mapping.maskedDofNosLocal.clear();

      // loop over the values to be set
      int i = 0;
      for (std::vector<double>::iterator iter = valuesToSendToRanks[ownRankNo].begin(); iter != valuesToSendToRanks[ownRankNo].end(); iter++, i++)
      {
        // only consider value if it is positive
        double value = *iter;

        if (value > 0)
        {
          mapping.maskedDofNosLocal.push_back(mapping.receivedValueDofNosLocal[i]);
          receivedValues.push_back(value);
        }
      }

      // set values, only where dofs are != -1
      slotSetValues(mapping.outputConnectorSlotNoTo, mapping.outputConnectorArrayIndexTo, mapping.maskedDofNosLocal, receivedValues, INSERT_VALUES);
    }
  }
}

//! return the data object of the timestepping scheme
template<typename FunctionSpaceType, typename NestedSolverType>
typename MapDofs<FunctionSpaceType,NestedSolverType>::Data &MapDofs<FunctionSpaceType,NestedSolverType>::
data()
{
  return data_;
}

//! Get the data that will be transferred in the operator splitting or coupling to the other term of the splitting/coupling.
//! the transfer is done by the output_connector_data_transfer class
template<typename FunctionSpaceType, typename NestedSolverType>
std::shared_ptr<typename MapDofs<FunctionSpaceType,NestedSolverType>::OutputConnectorDataType> MapDofs<FunctionSpaceType,NestedSolverType>::
getOutputConnectorData()
{
  return data_.getOutputConnectorData();
}

}  // namespace Control
