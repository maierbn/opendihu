#pragma once

#include <Python.h>  // has to be the first included header

#include "slot_connection/slot_connection.h"

/** This class stores the data slot connections that are given under "connectedSlots".
 *
 */
class GlobalConnectionsBySlotName
{
public:

  //! parse connection settings from python settings
  GlobalConnectionsBySlotName(PythonConfig settings);

  //! add all parsed slot connections to the outputConnection_ object of a splitting scheme
  template<typename SlotConnectorDataType1, typename SlotConnectorDataType2>
  void addConnections(
    std::shared_ptr<std::tuple<std::shared_ptr<SlotConnectorDataType1>,std::shared_ptr<SlotConnectorDataType2>>> slotConnectorData,
    std::shared_ptr<SlotConnection> outputConnection
  );

  //! get a descriptive string about all stored global connections that will be used as the first section of the solver_structure file
  std::string getDescriptionForDiagram();

private:

  std::vector<std::pair<std::string,std::string>> connections_;   //< connection from a slot to another slot, identified by the slot names
};

#include "slot_connection/global_connections_by_slot_name.tpp"
