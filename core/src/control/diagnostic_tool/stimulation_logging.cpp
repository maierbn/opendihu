#include "control/diagnostic_tool/stimulation_logging.h"

#include <Python.h>  // has to be the first included header

#include "output_writer/generic.h"

namespace Control
{

std::string StimulationLogging::filename_;
std::vector<StimulationLogging::StimulationLogEntry> StimulationLogging::logEntries_;

StimulationLogging::StimulationLogging(PythonConfig specificSettings)
{
  // get python config
  filename_ = specificSettings.getOptionString("stimulationLogFilename", "out/stimulation.log");
}

//! register the event of starting a stimulation for given motor unit, called by FastMonodomainSolver
void StimulationLogging::logStimulationBegin(double currentTime, int motorUnitNo, int fiberNo)
{
  StimulationLogEntry logEntry;
  logEntry.time = currentTime;
  logEntry.fiberNo = fiberNo;
  logEntry.motorUnitNo = motorUnitNo;
  logEntries_.push_back(logEntry);
}

//! this will be called at the end of the simulation run
void StimulationLogging::writeLogFile()
{
  if (filename_.empty())
  {
    LOG(DEBUG) << "Don't write StimulationLogging file, because filename is empty.";
    return;
  }

  LOG(DEBUG) << "StimulationLogging::writeLogFile, local entries " << logEntries_;

  // reduce all entries to rank 0
  // communicate number of entries on the ranks
  int ownSize = logEntries_.size();
  std::vector<int> sizesOnRanks;

  int ownRankNo = DihuContext::ownRankNoCommWorld();
  int nRanks = DihuContext::nRanksCommWorld();

  sizesOnRanks.resize(nRanks);
  //std::cout << "at StimulationLogging::logStimulationBegin, ownRankNo: " << ownRankNo << ", nRanks: " << nRanks << std::endl;

  MPIUtility::handleReturnValue(MPI_Allgather(&ownSize, 1, MPI_INT, sizesOnRanks.data(), 1, MPI_INT, MPI_COMM_WORLD), "MPI_Allgather");

  // count total number of entries
  int nTotalEntries = 0;
  for (int rankNo = 0; rankNo < nRanks; rankNo++)
  {
    nTotalEntries += sizesOnRanks[rankNo];
  }

  LOG(DEBUG) << "nTotalEntries: " << nTotalEntries << ", sizesOnRanks: " << sizesOnRanks;

  if (nTotalEntries == 0)
    return;

  // communicate all values to rank 0
  std::vector<double> times(nTotalEntries);
  std::vector<int> fiberNos(nTotalEntries);
  std::vector<int> motorUnitNos(nTotalEntries);

  std::vector<double> ownTimes(logEntries_.size());
  std::vector<int> ownFiberNos(logEntries_.size());
  std::vector<int> ownMotorUnitNos(logEntries_.size());

  int i = 0;
  for (std::vector<StimulationLogEntry>::iterator iter = logEntries_.begin(); iter != logEntries_.end(); iter++, i++)
  {
    ownTimes[i] = iter->time;
    ownFiberNos[i] = iter->fiberNo;
    ownMotorUnitNos[i] = iter->motorUnitNo;
  }

  LOG(DEBUG) << "ownTimes: " << ownTimes;

  // setup offsets for MPI_Gatherv
  std::vector<int> offsets(nRanks);

  for (int rankNo = 0; rankNo < nRanks; rankNo++)
  {
    if (rankNo == 0)
    {
      offsets[rankNo] = 0;
    }
    else
    {
      offsets[rankNo] = offsets[rankNo-1] + sizesOnRanks[rankNo-1];
    }
  }

  LOG(DEBUG) << "sizesOnRanks: " << sizesOnRanks << ", offsets: " << offsets;
  LOG(DEBUG) << "ownTimes: " << ownTimes << ", ownFiberNos: " << ownFiberNos << ", ownMotorUnitNos: " << ownMotorUnitNos;

  MPIUtility::handleReturnValue(MPI_Gatherv(ownTimes.data(),        sizesOnRanks[ownRankNo], MPI_DOUBLE, times.data(),        sizesOnRanks.data(), offsets.data(), MPI_DOUBLE, 0, MPI_COMM_WORLD), "MPI_Gatherv");
  MPIUtility::handleReturnValue(MPI_Gatherv(ownFiberNos.data(),     sizesOnRanks[ownRankNo], MPI_INT,    fiberNos.data(),     sizesOnRanks.data(), offsets.data(), MPI_INT, 0, MPI_COMM_WORLD), "MPI_Gatherv");
  MPIUtility::handleReturnValue(MPI_Gatherv(ownMotorUnitNos.data(), sizesOnRanks[ownRankNo], MPI_INT,    motorUnitNos.data(), sizesOnRanks.data(), offsets.data(), MPI_INT, 0, MPI_COMM_WORLD), "MPI_Gatherv");

  // on rank 0, write log file
  if (ownRankNo == 0)
  {
    LOG(DEBUG) << "times: " << times << ", fiberNos: " << fiberNos << ", motorUnitNos: " << motorUnitNos;

    // add received entries in logEntries_;
    logEntries_.clear();
    logEntries_.reserve(nTotalEntries);

    for (int i = 0; i < nTotalEntries; i++)
    {
      StimulationLogEntry logEntry;
      logEntry.time = times[i];
      logEntry.fiberNo = fiberNos[i];
      logEntry.motorUnitNo = motorUnitNos[i];
      VLOG(1) << "add entry " << i << ": " << logEntry;
      logEntries_.push_back(logEntry);
    }

    // sort entries in logEntries_
    // sort according to motor unit no
    std::sort(logEntries_.begin(), logEntries_.end(), [](const StimulationLogEntry &a, const StimulationLogEntry &b)
    {
      if (a.fiberNo == b.fiberNo && a.motorUnitNo == b.motorUnitNo)
      {
        return a.time - b.time < 0;
      }
      else if (a.motorUnitNo == b.motorUnitNo)
      {
        return a.fiberNo - b.fiberNo < 0;
      }
      else
      {
        return a.motorUnitNo - b.motorUnitNo < 0;
      }
    });

    LOG(DEBUG) << "StimulationLogging::writeLogFile, all entries " << logEntries_;

    std::ofstream file;
    OutputWriter::Generic::openFile(file, filename_);  // open file, and create directory if necessary, truncate file

    // output all entries in logEntries_
    file << "# motor unit no; fiber no; stimulation times in ms" << std::endl;

    int currentFiberNo = -2;
    int currentMotorUnitNo = -2;
    for (int i = 0; i < logEntries_.size(); i++)
    {
      if (currentFiberNo != logEntries_[i].fiberNo || currentMotorUnitNo != logEntries_[i].motorUnitNo)
      {
        currentFiberNo = logEntries_[i].fiberNo;
        currentMotorUnitNo = logEntries_[i].motorUnitNo;

        if (i != 0)
          file << "\n";
        file << currentMotorUnitNo << ";" << currentFiberNo << ";";
      }
      file << logEntries_[i].time << ";";
    }

    LOG(INFO) << "Wrote " << logEntries_.size() << " stimulation time" << (logEntries_.size() != 1? "s" : "")
      << " to file \"" << filename_ << "\".";
    file.close();
  }
}

std::ostream &operator<<(std::ostream &stream, const StimulationLogging::StimulationLogEntry rhs)
{
  stream << "(t:" << rhs.time << ",m:" << rhs.motorUnitNo << ",f:" << rhs.fiberNo << ")";
  return stream;
}

}  // namespace
