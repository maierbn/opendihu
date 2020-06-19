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
  performMappings(mappingsBeforeComputation_, nestedSolver_.startTime());

  // compute the simulation in the current time span with the nested solver
  nestedSolver_.advanceTimeSpan();

  LOG(DEBUG) << "MapDofs::performMappings afterComputation";
  // perform mapping from settings "afterComputation"
  performMappings(mappingsAfterComputation_, nestedSolver_.endTime());
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
performMappings(std::vector<DofsMappingType> &mappings, double currentTime)
{
  int ownRankNo = data_.functionSpace()->meshPartition()->rankSubset()->ownRankNo();

  for (DofsMappingType &mapping : mappings)
  {
    std::string modeString;
    switch(mapping.mode)
    {
    case DofsMappingType::modeCopyLocal:
      modeString = "modeCopyLocal";
      break;

    case DofsMappingType::modeCopyLocalIfPositive:
      modeString = "modeCopyLocalIfPositive";
      break;

    case DofsMappingType::modeLocalSetIfAboveThreshold:
      modeString = "modeLocalSetIfAboveThreshold";
      break;

    case DofsMappingType::modeCallback:
      modeString = "modeCallback";
      break;

    case DofsMappingType::modeCommunicate:
      modeString = "modeCommunicate";
      break;
    }

    LOG(DEBUG) << "MapDofs::perform mapping slots " << mapping.outputConnectorSlotNoFrom << " -> " << mapping.outputConnectorSlotNoTo << ", " << modeString;

    // static variables for input and output of values
    static std::map<int,std::vector<double>> valuesToSendToRanks;
    std::vector<double> &inputValues = valuesToSendToRanks[ownRankNo];
    static std::vector<double> valuesToSet;

    valuesToSet.clear();
    inputValues.clear();

    if (mapping.mode == DofsMappingType::modeCallback)
    {
      // get values of input dofs
      slotGetValues(mapping.outputConnectorSlotNoFrom, mapping.outputConnectorArrayIndexFrom, mapping.inputDofs, inputValues);

      // python call back signature:
      // callback(input_values, output_values, current_time, slot_nos, buffer)

      // compose callback function
      PyObject *inputValuesPy = PythonUtility::convertToPython<std::vector<double>>::get(inputValues);
      PyObject *globalParametersDict = PyDict_New();
      PyObject *arglist = Py_BuildValue("(O,O,d,O,O)", inputValuesPy, mapping.outputValuesPy, currentTime,
                                        mapping.slotNosPy, mapping.buffer);
      PyObject *returnValue = PyObject_CallObject(mapping.callback, arglist);

      // if there was an error while executing the function, print the error message
      if (returnValue == NULL)
        PythonUtility::checkForError();

      // select the dofs for which the output values were set to something different than None
      valuesToSet.clear();
      mapping.dofNosToSetLocal.clear();

      // loop over the output values
      int nOutputValues = PyList_Size(mapping.outputValuesPy);
      for (int outputValueNo = 0; outputValueNo < nOutputValues; outputValueNo++)
      {
        PyObject *outputValue = PyList_GetItem(mapping.outputValuesPy, outputValueNo);

        if (outputValue != Py_None)
        {
          double value = PythonUtility::convertFromPython<double>::get(outputValue);
          dof_no_t dofNoLocal = mapping.outputDofs[outputValueNo];
          valuesToSet.push_back(value);
          mapping.dofNosToSetLocal.push_back(dofNoLocal);
        }
      }

      // decrement reference counters for python objects
      Py_CLEAR(globalParametersDict);
      Py_CLEAR(returnValue);
      Py_CLEAR(arglist);
      Py_CLEAR(inputValuesPy);

      // set values in target field variable
      slotSetValues(mapping.outputConnectorSlotNoTo, mapping.outputConnectorArrayIndexTo, mapping.dofNosToSetLocal, valuesToSet, INSERT_VALUES);
    }
    else
    {

      // get input values from the selected slot
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

      LOG(INFO) << "get values " << valuesToSendToRanks;

      // if communication is involved
      if (mapping.mode == DofsMappingType::modeCommunicate)
      {
        mapping.valueCommunicator.communicate(valuesToSendToRanks, valuesToSet);

        // set values
        // set valuesToSet at dofs allDofNosToSetLocal
        slotSetValues(mapping.outputConnectorSlotNoTo, mapping.outputConnectorArrayIndexTo, mapping.allDofNosToSetLocal, valuesToSet, INSERT_VALUES);
      }
      else if (mapping.mode == DofsMappingType::modeCopyLocal)
      {
        // copy all values
        std::vector<double> &valuesToSet = inputValues;

        // set values
        // set valuesToSet at dofs allDofNosToSetLocal
        slotSetValues(mapping.outputConnectorSlotNoTo, mapping.outputConnectorArrayIndexTo, mapping.allDofNosToSetLocal, valuesToSet, INSERT_VALUES);
      }
      else if (mapping.mode == DofsMappingType::modeCopyLocalIfPositive)
      {
        // store dof nos for positive values in dofNosToSetLocal
        mapping.dofNosToSetLocal.clear();

        // loop over the values to be set
        int i = 0;
        for (std::vector<double>::iterator iter = inputValues.begin(); iter != inputValues.end(); iter++, i++)
        {
          // only consider value if it is positive
          double value = *iter;

          if (value > 0)
          {
            mapping.dofNosToSetLocal.push_back(mapping.allDofNosToSetLocal[i]);
            valuesToSet.push_back(value);
          }
        }

        // set values
        slotSetValues(mapping.outputConnectorSlotNoTo, mapping.outputConnectorArrayIndexTo, mapping.dofNosToSetLocal, valuesToSet, INSERT_VALUES);
      }
      else if (mapping.mode == DofsMappingType::modeLocalSetIfAboveThreshold)
      {
        // store dof nos for positive values in dofNosToSetLocal
        mapping.dofNosToSetLocal.clear();

        // loop over the values to be set
        int i = 0;
        for (std::vector<double>::iterator iter = inputValues.begin(); iter != inputValues.end(); iter++, i++)
        {
          double value = *iter;

          if (value > mapping.thresholdValue)
          {
            mapping.dofNosToSetLocal.push_back(mapping.allDofNosToSetLocal[i]);
            valuesToSet.push_back(mapping.valueToSet);
          }
        }

        LOG(DEBUG) << "values: " << inputValues << ", dofs that are above threshold " << mapping.thresholdValue
          << ", values to be set: " << valuesToSet;

        // set values
        slotSetValues(mapping.outputConnectorSlotNoTo, mapping.outputConnectorArrayIndexTo, mapping.dofNosToSetLocal, valuesToSet, INSERT_VALUES);
      }
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
