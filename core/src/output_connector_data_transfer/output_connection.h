#pragma once

#include <Python.h>  // has to be the first included header

#include "control/python_config/python_config.h"
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

  //! copy constructor
  OutputConnection(const OutputConnection &rhs);

  //! set the number of field variable components between which data will be transferred, this has to be done initially
  template<typename FunctionSpaceType1, int nComponents1a, int nComponents1b, typename FunctionSpaceType2, int nComponents2a, int nComponents2b>
  void initialize(const Data::OutputConnectorData<FunctionSpaceType1,nComponents1a,nComponents1b> &transferableSolutionData1,
                  const Data::OutputConnectorData<FunctionSpaceType2,nComponents2a,nComponents2b> &transferableSolutionData2,
                  int offsetSlotNoData1=0, int offsetSlotNoData2=0);

  //! set current transfer direction that will be taken into account for mapTo
  void setTransferDirection(bool term1To2);

  //! get the information to which slot the slot (fromVectorNo, fromIndex) should be mapped, @return: if there was no error, if it returns false, do not perform this mapping as the slot is not connected
  bool getSlotInformation(int fromVectorNo, int fromVectorIndex, int &toVectorNo, int &toVectorIndex, bool &avoidCopyIfPossible, bool disableWarnings=false) const;

  /** Identifies one slot where something can be connected to
   */
  struct Connector
  {
    int index;                  //< index in the list of slots, -1 means no connection
    bool avoidCopyIfPossible;   //< if the field variable should be reused
  };

  //! get the connectors from term 1 to term 2
  const std::vector<Connector> &connectorForVisualizerTerm1To2() const;

  //! get the connectors from term 2 to term 1
  const std::vector<Connector> &connectorForVisualizerTerm2To1() const;

  //! a pointer of a second output connection, used when the OutputConnectorData is a tuple of two OutputConnectorData types and therefore another output connection object is needed.
  std::shared_ptr<OutputConnection> &subOutputConnection1();

  //! one more pointer of a second output connection, used for transfer betwen two tuples and therefore a second output connection object is needed.
  std::shared_ptr<OutputConnection> &subOutputConnection2();

  //! one more pointer of a second output connection, used for transfer betwen two tuples and therefore a third output connection object is needed.
  std::shared_ptr<OutputConnection> &subOutputConnection3();

  //! one more pointer of a second output connection, used for transfer betwen two tuples and therefore a forth output connection object is needed.
  std::shared_ptr<OutputConnection> &subOutputConnection4();

  //! assemble some debugging information to the mapping that will be displayed on error
  std::string getDebugInformation() const;

private:

  //! initialize the slotInformation_ variable
  template<typename FunctionSpaceType1, int nComponents1a, int nComponents1b, typename FunctionSpaceType2, int nComponents2a, int nComponents2b>
  void initializeSlotInformation(const Data::OutputConnectorData<FunctionSpaceType1,nComponents1a,nComponents1b> &transferableSolutionData1,
                                 const Data::OutputConnectorData<FunctionSpaceType2,nComponents2a,nComponents2b> &transferableSolutionData2);


  //! fill the look-up table slotInformation_
  void initializeSlotInformation();

  std::vector<Connector> connectorTerm1To2_;    //< the connector information which variables to map to which for mapping from term 1 to term 2, this differs from connectorForVisualizerTerm1To2_ in that it respects the offsets
  std::vector<Connector> connectorTerm2To1_;    //< the connector information for mapping from term 2 to term 1

  std::vector<Connector> connectorForVisualizerTerm1To2_;    //< the connector information which variables to map to which for mapping from term 1 to term 2, without considering the offset. This is used by the solverStructureVisualizer to know which slots are connected. It cannot be used for the actual mapping.
  std::vector<Connector> connectorForVisualizerTerm2To1_;    //< the connector information for mapping from term 2 to term 1. This is used by the solverStructureVisualizer to know which slots are connected. It cannot be used for the actual mapping.

  std::vector<std::string> fieldVariableNamesTerm1Vector1_;   //< only for debugging the field variable names and components nos of term 1 variable1
  std::vector<std::string> fieldVariableNamesTerm1Vector2_;   //< only for debugging the field variable names and components nos of term 1 variable2
  std::vector<std::string> fieldVariableNamesTerm2Vector1_;   //< only for debugging the field variable names and components nos of term 2 variable1
  std::vector<std::string> fieldVariableNamesTerm2Vector2_;   //< only for debugging the field variable names and components nos of term 2 variable2
  bool fieldVariableNamesInitialized_;                        //< if the fieldVAriableNamesTerm* variables were initialized yet by initialize

  int nFieldVariablesTerm1Vector1_; //< the number of slots of term 1 in vector 1
  int nFieldVariablesTerm1Vector2_; //< the number of slots of term 1 in vector 2
  int nFieldVariablesTerm2Vector1_; //< the number of slots of term 2 in vector 1
  int nFieldVariablesTerm2Vector2_; //< the number of slots of term 2 in vector 2
  bool transferDirectionTerm1To2_;  //< if the current mapping is from term 1 to 2

  int offsetSlotNoData1_;           //< an offset value for the slot no on transferableSolutionData1, this is the first slot no
  int offsetSlotNoData2_;

  struct Result
  {
    int toVectorNo;
    int toVectorIndex;
    bool avoidCopyIfPossible;
    bool successful;
  };

  std::array<std::array<std::vector<Result>,2>,2> slotInformation_;   //< [transferDirectionTerm1To2_][fromVectorNo][fromVectorIndex], look-up table of getSlotInformation
  bool slotInformationInitialized_;                                   //< if slotInformation has been initialized

  PythonConfig settings_;                                             //< the settings object
  std::shared_ptr<OutputConnection> subOutputConnection1_;             //< a first output connection object
  std::shared_ptr<OutputConnection> subOutputConnection2_;             //< a second output connection object
  std::shared_ptr<OutputConnection> subOutputConnection3_;            //< a third output connection object
  std::shared_ptr<OutputConnection> subOutputConnection4_;            //< a forth output connection object
};

#include "output_connector_data_transfer/output_connection.tpp"
