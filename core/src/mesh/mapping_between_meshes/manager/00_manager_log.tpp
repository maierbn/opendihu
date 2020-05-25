#include "mesh/mapping_between_meshes/manager/00_manager_log.h"

#include "control/dihu_context.h"

namespace MappingBetweenMeshes
{

//! add a log entry for two field variables
template<typename FieldVariableSourceType, typename FieldVariableTargetType>
void ManagerLog::
addLogEntryFieldVariable(std::shared_ptr<FieldVariableSourceType> fieldVariableSource, int componentNoSource,
                         std::shared_ptr<FieldVariableTargetType> fieldVariableTarget, int componentNoTarget, mappingLogEntry_t::logEvent_t logEvent)
{
  // only log mapping actions on rank 0
  if (DihuContext::ownRankNoCommWorld() != 0)
    return;

  /*struct mappingLogEntry_t
  {
    std::string fieldVariableNameFrom;
    std::string meshNameFrom;
    int componentNoFrom;
    int dimensionalityFrom;
    
    std::string fieldVariableNameTo;
    std::string meshNameTo;
    int componentNoTo;
    int dimensionalityTo;

    int nMappedSourceMeshes;
    int nRepetitions;
    enum logEvent_t
    {
      eventParseSettings,
      eventCreateMapping,
      eventMapForward,
      eventMapReverse
    }
    logEvent;
  };*/

  mappingLogEntry_t logEntry;

  // set variables corresponding to the source field variable
  if (fieldVariableSource)
  {
    logEntry.fieldVariableNameFrom = fieldVariableSource->name();
    logEntry.meshNameFrom = fieldVariableSource->functionSpace()->meshName();
    logEntry.componentNoFrom = componentNoSource;
    logEntry.dimensionalityFrom = fieldVariableSource->functionSpace()->dim();
  }
  
  // set variables corresponding to the target field variable
  if (fieldVariableTarget)
  {
    logEntry.fieldVariableNameTo = fieldVariableTarget->name();
    logEntry.meshNameTo = fieldVariableTarget->functionSpace()->meshName();
    logEntry.componentNoTo = componentNoTarget;
    logEntry.dimensionalityTo = fieldVariableTarget->functionSpace()->dim();
  }
  
  logEntry.logEvent = logEvent;

  // check if action was already logged previously
  if (!logEntries_.empty())
  {
    if (logEntries_.back().fieldVariableNameFrom == logEntry.fieldVariableNameFrom
        && logEntries_.back().componentNoFrom == logEntry.componentNoFrom
        && logEntries_.back().fieldVariableNameTo == logEntry.fieldVariableNameTo
        && logEntries_.back().componentNoTo == logEntry.componentNoTo
        && logEntries_.back().logEvent == logEvent)
    {
      // last log entry is the same as the current action, repeat is repetitions counter and return
      logEntries_.back().nRepetitions++;
      return;
    }
  }

  logEntry.nMappedSourceMeshes = mappedSourceMeshesCounter_;
  logEntry.nRepetitions = 0;

  // add entry to logs
  logEntries_.push_back(logEntry);
}


//! add a log entry when a new mapping is created or initialized
template<typename FunctionSpaceSourceType, typename FunctionSpaceTargetType>
void ManagerLog::
addLogEntryMapping(std::shared_ptr<FunctionSpaceSourceType> functionSpaceSource,
                   std::shared_ptr<FunctionSpaceTargetType> functionSpaceTarget, mappingLogEntry_t::logEvent_t logEvent)
{
  // only log mapping actions on rank 0
  if (DihuContext::ownRankNoCommWorld() != 0)
    return;

  mappingLogEntry_t logEntry;

  // set variables corresponding to the source field variable
  if (functionSpaceSource)
  {
    logEntry.meshNameFrom = functionSpaceSource->meshName();
    logEntry.dimensionalityFrom = functionSpaceSource->dim();
  }
  
  // set variables corresponding to the target field variable
  if (functionSpaceTarget)
  {
    logEntry.meshNameTo = functionSpaceTarget->meshName();
    logEntry.dimensionalityTo = functionSpaceTarget->dim();
  }
  
  logEntry.logEvent = logEvent;

  // check if action was already logged previously
  if (!logEntries_.empty())
  {
    if (logEntries_.back().meshNameFrom == logEntry.meshNameFrom
        && logEntries_.back().meshNameTo == logEntry.meshNameTo
        && logEntries_.back().logEvent == logEvent)
    {
      // last log entry is the same as the current action, repeat is repetitions counter and return
      logEntries_.back().nRepetitions++;
      return;
    }
  }

  logEntry.nMappedSourceMeshes = 0;
  logEntry.nRepetitions = 0;

  // add entry to logs
  logEntries_.push_back(logEntry);
}

}   // namespace
