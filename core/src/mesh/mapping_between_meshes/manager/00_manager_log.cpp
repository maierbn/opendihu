#include "mesh/mapping_between_meshes/manager/00_manager_log.h"

#include "output_writer/generic.h"

namespace MappingBetweenMeshes
{


ManagerLog::ManagerLog(PythonConfig specificSettings) :
  specificSettings_(specificSettings)
{
  logFilename_ = specificSettings_.getOptionString("mappingsBetweenMeshesLogFile", "out/mappings_between_meshes.txt");
}


void ManagerLog::addLogEntryMapping(std::string meshNameSource, std::string meshNameTarget, mappingLogEntry_t::logEvent_t logEvent)
{
  // only log mapping actions on rank 0
  if (DihuContext::ownRankNoCommWorld() != 0)
    return;

  mappingLogEntry_t logEntry;

  logEntry.meshNameFrom = meshNameSource;
  logEntry.dimensionalityFrom = -1;
  logEntry.meshNameTo = meshNameTarget;
  logEntry.dimensionalityTo = -1;
  
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

std::string ManagerLog::produceLogContents()
{
  std::stringstream log;

  for (const mappingLogEntry_t &logEntry : logEntries_)
  {
    switch (logEntry.logEvent)
    {
    case mappingLogEntry_t::logEvent_t::eventParseSettings:
      log << "parse settings for mapping between meshes \"" << logEntry.meshNameFrom << "\" " 
        << "-> \"" << logEntry.meshNameTo << "\"";
      break;
    
    case mappingLogEntry_t::logEvent_t::eventCreateMapping:
      log << "create mapping between meshes \"" << logEntry.meshNameFrom << "\" (" << logEntry.dimensionalityFrom << "D) " 
        << "-> \"" << logEntry.meshNameTo << "\" (" << logEntry.dimensionalityTo << "D)";
      break;

    case mappingLogEntry_t::logEvent_t::eventMapForward:
    case mappingLogEntry_t::logEvent_t::eventMapReverse:
      log << "map from field variable \"" << logEntry.fieldVariableNameFrom << "\"";
      if (logEntry.componentNoFrom == -1)
        log << ", all components";
      else 
        log << " component " << logEntry.componentNoFrom;
      
      log << " (mesh \"" << logEntry.meshNameFrom << "\" (" << logEntry.dimensionalityFrom << "D) ";
      
      if (logEntry.nMappedSourceMeshes > 1)
        log << " [and other field variables/meshes, in total " << logEntry.nMappedSourceMeshes << " ]";
        
      log << " to field variable \"" << logEntry.fieldVariableNameTo << "\"";
      if (logEntry.componentNoTo == -1)
        log << ", all components";
      else 
        log << " component " << logEntry.componentNoTo;

      log << " (mesh \"" << logEntry.meshNameTo << "\" (" << logEntry.dimensionalityTo << "D)";
      
      if (logEntry.logEvent == mappingLogEntry_t::logEvent_t::eventMapForward)
      {
        log << " using forward mapping \"" << logEntry.meshNameFrom << "\" (" << logEntry.dimensionalityFrom << "D) -> "
          "\"" << logEntry.meshNameTo << "\" (" << logEntry.dimensionalityTo << "D)";
      }
      else 
      {
        log << " using reverse direction of mapping \"" << logEntry.meshNameTo << "\" (" << logEntry.dimensionalityTo << "D) -> "
          "\"" << logEntry.meshNameFrom << "\" (" << logEntry.dimensionalityFrom << "D)";
      }
      break;
    }
    if (logEntry.nRepetitions > 0)
    {
      log << " (Repeated " << logEntry.nRepetitions << " times)";
    }
    log << "\n";
  }

  return log.str();
}

//! produce the resulting file, only if solverAddingEnabled_
void ManagerLog::
writeLogFile()
{
  // only produce file on rank 0
  if (DihuContext::ownRankNoCommWorld() == 0)
  {
    if (!logFilename_.empty())
    {
      std::string diagram = produceLogContents();

      // output diagram to a text file
      std::ofstream file;
      OutputWriter::Generic::openFile(file, logFilename_);

      if (!file.is_open())
      {
        LOG(FATAL) << "Could not write to file \"" << logFilename_ << "\".";
      }

      // output diagram text
      file << diagram;
      file.close();

      LOG(INFO) << "File \"" << logFilename_ << "\" written.";
    }
  }
}

}   // namespace
