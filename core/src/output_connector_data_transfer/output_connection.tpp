#include "output_connector_data_transfer/output_connection.h"

template<typename FunctionSpaceType1, int nComponents1a, int nComponents1b, typename FunctionSpaceType2, int nComponents2a, int nComponents2b>
void OutputConnection::
initialize(const Data::OutputConnectorData<FunctionSpaceType1,nComponents1a,nComponents1b> &transferableSolutionData1,
           const Data::OutputConnectorData<FunctionSpaceType2,nComponents2a,nComponents2b> &transferableSolutionData2)
{
  if (transferDirectionTerm1To2_)
  {
    nFieldVariablesTerm1Vector1_ = transferableSolutionData1.variable1.size();
    nFieldVariablesTerm1Vector2_ = transferableSolutionData1.variable2.size();
    nFieldVariablesTerm2Vector1_ = transferableSolutionData2.variable1.size();
    nFieldVariablesTerm2Vector2_ = transferableSolutionData2.variable2.size();
  }
  else
  {
    nFieldVariablesTerm1Vector1_ = transferableSolutionData2.variable1.size();
    nFieldVariablesTerm1Vector2_ = transferableSolutionData2.variable2.size();
    nFieldVariablesTerm2Vector1_ = transferableSolutionData1.variable1.size();
    nFieldVariablesTerm2Vector2_ = transferableSolutionData1.variable2.size();
  }

  if (!fieldVariableNamesInitialized_)
  {
    fieldVariableNamesInitialized_ = true;

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

    initializeSlotInformation();
  }
}
