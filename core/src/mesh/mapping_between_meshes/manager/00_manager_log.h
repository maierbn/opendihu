#pragma once

#include <Python.h>  // has to be the first included header
#include <map>
#include <vector>

#include "field_variable/00_field_variable_base.h"
#include "control/python_config/python_config.h"

#include <iostream>

namespace MappingBetweenMeshes
{
  
class ManagerLog
{
public:
  //! constructor
  ManagerLog(PythonConfig specificSettings);

  //! produce the resulting file, only if solverAddingEnabled_
  static void writeLogFile();

  //! add a message
  void addLogMessage(std::string message);

protected:

  /** An event that will be written to the mapping log
   */
  struct mappingLogEntry_t
  {
    std::string fieldVariableNameFrom;
    std::string meshNameFrom;
    int componentNoFrom;
    int dimensionalityFrom;               //< dimensionality of the source mesh, i.e. 1D, 2D or 3D
    
    std::string fieldVariableNameTo;
    std::string meshNameTo;
    int componentNoTo;
    int dimensionalityTo;                 //< dimensionality of the target mesh, i.e. 1D, 2D or 3D

    int nMappedSourceMeshes;              //< number of source meshes that were mapped to the target mesh between calls to prepareMapping() and finalizeMapping()
    int nRepetitions;                     //< how often this action was repeated in sequence
    std::string message;                  //< if logEvent == eventMessage, the message to be shown
    
    enum logEvent_t
    {
      eventParseSettings,                 //< a mapping between "from" and "to" was read from the settings
      eventCreateMapping,                 //< a mapping between "from" and "to" was created
      eventMapForward,                    //< mapping was performed between "from" and "to", using the corresponding mapping
      eventMapReverse,                    //< mapping was performed between "from" and "to", using the reverse direction of the reverse mapping
      eventMessage                        //< additional message to be displayed at the current location in the log
    }
    logEvent;
  };

  //! add a log entry for two field variables
  template<typename FieldVariableSourceType, typename FieldVariableTargetType>
  void addLogEntryFieldVariable(std::shared_ptr<FieldVariableSourceType> fieldVariableSource, int componentNoSource,
                                std::shared_ptr<FieldVariableTargetType> fieldVariableTarget, int componentNoTarget, mappingLogEntry_t::logEvent_t logEvent);


  //! add a log entry when a new mapping is created or initialized
  template<typename FunctionSpaceSourceType, typename FunctionSpaceTargetType>
  void addLogEntryMapping(std::shared_ptr<FunctionSpaceSourceType> functionSpaceSource,
                          std::shared_ptr<FunctionSpaceTargetType> functionSpaceTarget, mappingLogEntry_t::logEvent_t logEvent);

  //! add a log entry when a new mapping is created or initialized
  void addLogEntryMapping(std::string meshNameSource, std::string meshNameTarget, mappingLogEntry_t::logEvent_t logEvent);

  //! create the contents of the log file as string, from the logEntries_ variable
  static std::string produceLogContents();

  PythonConfig specificSettings_;         //< python object containing the value of the python config dict with corresponding key, for meshManager
  static std::string logFilename_;        //< filename of a log file that contains information about which mappings were created and used
  int mappedSourceMeshesCounter_;         //< variable that will be reset in prepareMapping, incremented in map() for every mesh and read out in finalizeMapping() to set the value of nMappedSourceMeshes for the log
  
  static std::vector<mappingLogEntry_t> logEntries_;   //< all entries that will be written to the log file
};

}  // namespace

#include "mesh/mapping_between_meshes/manager/00_manager_log.tpp"
