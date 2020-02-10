#pragma once

#include <Python.h>  // has to be the first included header

#include "control/python_config.h"
#include "data_management/output_connector_data.h"

/** This specifies the connections of data slots between two terms, Term 1 and Term 2.
 *  Each term is assumed to have data of type Data::OutputConnectorData which has the two vectors variable1 and variable2.
 *  Each variable vector can hold multiple components of field variables of the same type (number of components of the field variable).
 *  One such component is referred to as "slot". In the settings, there is the specification which slots to map to which slots by the fields
 *  "connectedSlotsTerm2To1" and "connectedSlotsTerm1To2".
 *
 */
class OutputConnection
{
public:
  //! parse connection settings from python settings
  OutputConnection(PythonConfig settings);

  //! set the number of field variable components between which data will be transferred, this has to be done initially
  template<typename FunctionSpaceType1, int nComponents1a, int nComponents1b, typename FunctionSpaceType2, int nComponents2a, int nComponents2b>
  void initialize(const Data::OutputConnectorData<FunctionSpaceType1,nComponents1a,nComponents1b> &transferableSolutionData1,
                  const Data::OutputConnectorData<FunctionSpaceType2,nComponents2a,nComponents2b> &transferableSolutionData2);

  //! set current tranfer direction that will be taken into account for mapTo
  void setTransferDirection(bool term1To2);

  //! get the information to which slot the slot (fromVectorNo, fromIndex) should be mapped, @return: if there was no error, if it returns false, do not perform this mapping as the slot is not connected
  bool getSlotInformation(int fromVectorNo, int fromVectorIndex, int &toVectorNo, int &toVectorIndex, bool &avoidCopyIfPossible, bool disableWarnings=false) const;

  /** Specifies one slot
   */
  struct Connector
  {
    int index;                  //< index in the list of slots, -1 means no connection
    bool avoidCopyIfPossible;   //< if the field variable should be reused
  };

  //! get the connectors from term 1 to term 2
  const std::vector<Connector> &connectorTerm1To2() const;

  //! get the connectors from term 2 to term 1
  const std::vector<Connector> &connectorTerm2To1() const;
private:

  //! assemble some debugging information to the mapping that will be displayed on error
  std::string getDebugInformation() const;

  //! fill the look-up table slotInformation_
  void initializeSlotInformation();

  std::vector<Connector> connectorTerm1To2_;    //< the connector information which variables to map to which for mapping from term 1 to term 2
  std::vector<Connector> connectorTerm2To1_;    //< the connector information for mapping from term 2 to term 1

  std::vector<std::string> fieldVariableNamesTerm1Vector1_;   //< only for debugging the field variable names and components nos of term 1 variable1
  std::vector<std::string> fieldVariableNamesTerm1Vector2_;   //< only for debugging the field variable names and components nos of term 1 variable2
  std::vector<std::string> fieldVariableNamesTerm2Vector1_;   //< only for debugging the field variable names and components nos of term 2 variable1
  std::vector<std::string> fieldVariableNamesTerm2Vector2_;   //< only for debugging the field variable names and components nos of term 2 variable2
  bool fieldVariableNamesInitialized_;                        //< if the fieldVAriableNamesTerm* variables were initialized yet by initialize

  int nFieldVariablesTerm1Vector1_; //< the number of slots of term 1 in vector 1
  int nFieldVariablesTerm1Vector2_; //< the number of slots of term 1 in vector 2
  int nFieldVariablesTerm2Vector1_; //< the number of slots of term 2 in vector 1
  int nFieldVariablesTerm2Vector2_; //< the number of slots of term 2 in vector 2
  bool transferDirectionTerm1To2_;   //< if the current mapping is from term 1 to 2

  struct Result
  {
    int toVectorNo;
    int toVectorIndex;
    bool avoidCopyIfPossible;
    bool successful;
  };

  std::array<std::array<std::vector<Result>,2>,2> slotInformation_;   // [transferDirectionTerm1To2_][fromVectorNo][fromVectorIndex], look-up table of getSlotInformation
  bool slotInformationInitialized_;          //< if slotInformation has been initialized

  PythonConfig settings_;         //< the settings object
};

#include "output_connector_data_transfer/output_connection.tpp"
