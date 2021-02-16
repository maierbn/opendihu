#include "control/map_dofs/map_dofs.h"

#include "utility/python_capture_stderr.h"

namespace Control
{

//! advance simulation by the given time span
template<typename FunctionSpaceType, typename NestedSolverType>
void MapDofs<FunctionSpaceType,NestedSolverType>::
advanceTimeSpan(bool withOutputWritersEnabled)
{
  LOG(DEBUG) << "MapDofs::advanceTimeSpan, " << mappingsBeforeComputation_.size() << " mappings before computation, "
    << mappingsAfterComputation_.size() << " mappings after computation.";
  LOG(DEBUG) << "MapDofs::performMappings beforeComputation";

  // perform mapping from settings "beforeComputation"
  performMappings(mappingsBeforeComputation_, nestedSolver_.startTime());

  // compute the simulation in the current time span with the nested solver
  nestedSolver_.advanceTimeSpan(withOutputWritersEnabled);

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

//! start time of time interval to be simulated
template<typename FunctionSpaceType, typename NestedSolverType>
double MapDofs<FunctionSpaceType,NestedSolverType>::
startTime()
{
  return nestedSolver_.startTime();
}

//! end time of simulation
template<typename FunctionSpaceType, typename NestedSolverType>
double MapDofs<FunctionSpaceType,NestedSolverType>::
endTime()
{
  return nestedSolver_.endTime();
}

//! call the output writer on the data object, output files will contain currentTime, with callCountIncrement !=1 output timesteps can be skipped
template<typename FunctionSpaceType, typename NestedSolverType>
void MapDofs<FunctionSpaceType,NestedSolverType>::
callOutputWriter(int timeStepNo, double currentTime, int callCountIncrement)
{
  // call the output writer of the nested solver
  nestedSolver_.callOutputWriter(timeStepNo, currentTime, callCountIncrement);
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

    LOG(DEBUG) << "-> MapDofs::perform mapping slots " << mapping.connectorSlotNoFrom << " -> " << mapping.connectorSlotNosTo << ", " << modeString;

    // static variables for input and output of values
    static std::map<int,std::vector<double>> valuesToSendToRanks;
    std::vector<double> &inputValues = valuesToSendToRanks[ownRankNo];
    static std::vector<double> valuesToSet;

    valuesToSet.clear();
    inputValues.clear();

    if (mapping.mode == DofsMappingType::modeCallback)
    {
      // get values of input dofs
      slotGetValues(mapping.connectorSlotNoFrom, mapping.slotConnectorArrayIndexFrom, mapping.inputDofs, inputValues);

#ifndef NDEBUG
      std::string stdoutBuffer;

      // add callback function to capture stdout buffer
      emb::stdout_write_type stdoutWrite = [&stdoutBuffer] (std::string s) {stdoutBuffer += s; };
      emb::set_stdout(stdoutWrite);
#endif

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

#ifndef NDEBUG
      LOG(DEBUG) << "callback output: " << stdoutBuffer;
      emb::reset_stdout();
#endif

      //PyObject *outputValuesPy = returnValue;
      PyObject *outputValuesPy = mapping.outputValuesPy;

      int nToSlots = PyList_Size(outputValuesPy);
      // loop over the slots to which the output values should be transferred
      for (int toSlotIndex = 0; toSlotIndex < nToSlots; toSlotIndex++)
      {
        PyObject *outputValuesSlotPy = PyList_GetItem(outputValuesPy, toSlotIndex);

        // select the dofs for which the output values were set to something different than None
        valuesToSet.clear();
        mapping.dofNosToSetLocal.clear();

        // loop over the output values for the current slot
        int nOutputValues = PyList_Size(outputValuesSlotPy);
        for (int outputValueNo = 0; outputValueNo < nOutputValues; outputValueNo++)
        {
          PyObject *outputValue = PyList_GetItem(outputValuesSlotPy, outputValueNo);

          if (outputValue != Py_None)
          {
            double value = PythonUtility::convertFromPython<double>::get(outputValue);
            dof_no_t dofNoLocal = mapping.outputDofs[toSlotIndex][outputValueNo];
            valuesToSet.push_back(value);
            mapping.dofNosToSetLocal.push_back(dofNoLocal);
          }
        }

#ifndef NDEBUG
        LOG(DEBUG) << "   slot " << mapping.connectorSlotNosTo[toSlotIndex] << " (index " << toSlotIndex << "/" << mapping.connectorSlotNosTo.size() << ")"
          << ", set values from callback: " << valuesToSet << " at dofs: " << mapping.dofNosToSetLocal;
#endif

        // set values in target field variable
        slotSetValues(mapping.connectorSlotNosTo[toSlotIndex], mapping.slotConnectorArrayIndexTo, mapping.dofNosToSetLocal, valuesToSet, INSERT_VALUES);
      }

      // decrement reference counters for python objects
      Py_CLEAR(globalParametersDict);
      Py_CLEAR(returnValue);
      Py_CLEAR(arglist);
      Py_CLEAR(inputValuesPy);
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
        slotGetValues(mapping.connectorSlotNoFrom, mapping.slotConnectorArrayIndexFrom, dofNosLocal, valuesToSendToRanks[rankNo]);
      }

      // if communication is involved
      if (mapping.mode == DofsMappingType::modeCommunicate)
      {
        mapping.valueCommunicator.communicate(valuesToSendToRanks, valuesToSet);

        // set values
        // set valuesToSet at dofs allDofNosToSetLocal
        slotSetValues(mapping.connectorSlotNosTo[0], mapping.slotConnectorArrayIndexTo, mapping.allDofNosToSetLocal, valuesToSet, INSERT_VALUES);
      }
      else if (mapping.mode == DofsMappingType::modeCopyLocal)
      {
        // copy all values
        std::vector<double> &valuesToSet = inputValues;

        // set values
        // set valuesToSet at dofs allDofNosToSetLocal
        slotSetValues(mapping.connectorSlotNosTo[0], mapping.slotConnectorArrayIndexTo, mapping.allDofNosToSetLocal, valuesToSet, INSERT_VALUES);
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
        slotSetValues(mapping.connectorSlotNosTo[0], mapping.slotConnectorArrayIndexTo, mapping.dofNosToSetLocal, valuesToSet, INSERT_VALUES);
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

#ifndef NDEBUG
        LOG(DEBUG) << "   mapping.dofNosLocalOfValuesToSendToRanks: " << mapping.dofNosLocalOfValuesToSendToRanks;
        LOG(DEBUG) << "   values: " << inputValues << ", dofs that are above threshold " << mapping.thresholdValue
          << ": " << mapping.dofNosToSetLocal
          << ", values to be set: " << valuesToSet;
#endif

        // set values
        slotSetValues(mapping.connectorSlotNosTo[0], mapping.slotConnectorArrayIndexTo, mapping.dofNosToSetLocal, valuesToSet, INSERT_VALUES);
      }

      // call setRepresentationGlobal on the field variable, not needed
      //slotSetRepresentationGlobal(mapping.connectorSlotNosTo[0], mapping.slotConnectorArrayIndexTo);
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
//! the transfer is done by the slot_connector_data_transfer class
template<typename FunctionSpaceType, typename NestedSolverType>
std::shared_ptr<typename MapDofs<FunctionSpaceType,NestedSolverType>::SlotConnectorDataType> MapDofs<FunctionSpaceType,NestedSolverType>::
getSlotConnectorData()
{
  return data_.getSlotConnectorData();
}

}  // namespace Control
