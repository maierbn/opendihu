#include "control/diagnostic_tool/performance_measurement.h"

#include "easylogging++.h"
#include <memory>
#include <fstream>
#include <unistd.h>
#include <limits.h>
#include <ios>
#include <iostream>
#include <fstream>
#include <string>
#include <sys/param.h>
#include <iomanip>

#include "output_writer/generic.h"

namespace Control
{

void PerformanceMeasurement::writeLogFile(std::string logFileName)
{
  //LOG(DEBUG) << "PerformanceMeasurement::writeLogFile \"" << logFileName;

  parseStatusInformation();

  const bool useMPIOutput = true;   /// if the output is using MPI Output

  int ownRankNo = DihuContext::ownRankNoCommWorld();

  // determine file name
  std::stringstream filename;
  filename << logFileName;
  if (!useMPIOutput) {
    filename << "." << std::setw(7) << std::setfill('0') << ownRankNo;
  }

  // time stamp
  auto t = std::time(nullptr);
  auto tm = *std::localtime(&t);
  // host name
  char hostname[MAXHOSTNAMELEN+1];
  gethostname(hostname, MAXHOSTNAMELEN+1);

  std::stringstream header;
  std::stringstream data;

  auto logFormat = DihuContext::logFormat();

  if (logFormat == DihuContext::logFormatCsv)
  {
    filename << ".csv";
    logFileName = filename.str();

    // compose header
    header << "# timestamp;hostname;version;nRanks;rankNo;";

    std::set<std::string> measurementNames;
    for (std::pair<std::string, Measurement> measurement : measurements_)
    {
      measurementNames.insert(measurement.first);
    }

    // Add additional measurement names that could be only on some ranks and not on all.
    // The goal is to have the same header on every rank such that the log file gets consistent.
    measurementNames.insert(std::string("durationParaview1D"));
    measurementNames.insert(std::string("durationParaview1DInit"));
    measurementNames.insert(std::string("durationParaview2D"));
    measurementNames.insert(std::string("durationParaview3D"));
    measurementNames.insert(std::string("durationParaview3DInit"));
    measurementNames.insert(std::string("durationParaview3DReduction"));
    measurementNames.insert(std::string("durationParaview3DWrite"));
    measurementNames.insert(std::string("durationParaviewOutput"));
    measurementNames.insert(std::string("durationWriteOutputParaview"));

    // write measurement names
    for (std::string measurementName : measurementNames)
    {
      header << measurementName << ";n;";
    }

    // write parameter names
    for (std::pair<std::string,std::string> parameter : parameters_)
    {
      if (parameter.first == "nRanks" || parameter.first == "rankNo")
        continue;
      header << parameter.first << ";";
    }

    // write sum names
    for (std::pair<std::string,int> sum : sums_)
    {
      header << sum.first << ";";
    }

    // compose data
    data << StringUtility::timeToString(&tm) << ";";
    data << std::string(hostname) << ";";
    data << DihuContext::versionText() << ";";
    data << parameters_["nRanks"] << ";" << parameters_["rankNo"] << ";";

    // write measurement values
    for (std::string measurementName : measurementNames)
    {
      if (measurements_.find(measurementName) != measurements_.end())
      {
        data << measurements_[measurementName].totalDuration << ";"
          << measurements_[measurementName].nTimeSpans << ";";
      }
      else
      {
        data << "0.0;0;";
      }
    }

    // write parameters
    for (std::pair<std::string,std::string> parameter : parameters_)
    {
      if (parameter.first == "nRanks" || parameter.first == "rankNo")
        continue;

      // remove newlines
      StringUtility::replace(parameter.second, "\n", "");
      StringUtility::replace(parameter.second, "\r", "");
      data << parameter.second << ";";
    }

    // write sums
    for (std::pair<std::string,int> sum : sums_)
    {
      data << sum.second << ";";
    }

    data << std::endl;


  }
  else if (logFormat == DihuContext::logFormatJson)
  {
    filename << ".json";
    logFileName = filename.str();

    // compose header
    header << "#";
    header << std::endl;

    // compose data
    data << "{\"timestamp\":\"" << StringUtility::timeToString(&tm) << "\",";
    data << "\"hostname\":\"" << std::string(hostname) << "\",";
    data << "\"version\":\"" << DihuContext::versionText() << "\",";
    data << "\"nRanks\":" << parameters_["nRanks"] << ",\"rankNo\":" << parameters_["rankNo"];

    // write measurement values
    for (std::pair<std::string, Measurement> measurement : measurements_)
    {
      data << ",\"" << measurement.first << "\":" << measurement.second.totalDuration << ","
      << "\"" << measurement.first << " n\":" << measurement.second.nTimeSpans;
    }

    // write parameters
    for (std::pair<std::string,std::string> parameter : parameters_)
    {
      if (parameter.first == "nRanks" || parameter.first == "rankNo")
      continue;

      // remove newlines
      StringUtility::replace(parameter.second, "\n", "");
      StringUtility::replace(parameter.second, "\r", "");
      data << ",\"" << parameter.first << "\":\"" << parameter.second << '\"';
    }

    // write sums
    for (std::pair<std::string,int> sum : sums_)
    {
      // std::cout << "sum:  " << sum << " " << sum.first << " " << sum.second << std::endl;
      data << ",\"" << sum.first << "\":" << sum.second;
    }

    data << "}" << std::endl;
  }
  else
  {
    LOG(ERROR) << "BUG: Unknown log fromat '" << DihuContext::logFormat() << "'. This should not happen. NO log file is written!";
    return;
  }


  // check if header has to be added to file
  bool outputHeader = true;
  if (!useMPIOutput || ownRankNo == 0)
  {
    // parse header and check if it would be the same header as own
    std::ifstream inputFile(filename.str());
    if (inputFile.is_open())
    {
      std::string fileContent;
      fileContent.assign( (std::istreambuf_iterator<char>(inputFile) ),
                      (std::istreambuf_iterator<char>()    ) );
      if (fileContent.rfind("#") != std::string::npos)
      {
        std::size_t beginHeader = fileContent.rfind("#");
        std::size_t endHeader = fileContent.substr(beginHeader).find("\n");
        std::string lastHeader = fileContent.substr(beginHeader, endHeader);
        if (lastHeader == header.str())
        {
          outputHeader = false;
        }
      }
    }
  }

  // write header and data to file
  if (useMPIOutput)   // MPI output
  {
    // open log file to create directory if needed
    std::ofstream file;
    if (ownRankNo == 0)
    {
      OutputWriter::Generic::openFile(file, filename.str(), true);
      file.close();
    }

    LOG(DEBUG) << "open MPI file \"" << logFileName << "\".";

    // open file
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_File fileHandle;
    MPIUtility::handleReturnValue(MPI_File_open(MPI_COMM_WORLD, logFileName.c_str(),
                                                //MPI_MODE_WRONLY | MPI_MODE_CREATE | MPI_MODE_UNIQUE_OPEN,
                                                MPI_MODE_WRONLY | MPI_MODE_CREATE | MPI_MODE_APPEND,
                                                MPI_INFO_NULL, &fileHandle), "MPI_File_open");

    // write header
    // collective blocking write, only rank 0 writes, but afterwards all have the same shared file pointer position
    if (ownRankNo == 0)
    {
      MPI_Status status;
      MPIUtility::handleReturnValue(MPI_File_write_ordered(fileHandle, header.str().c_str(), header.str().length(), MPI_BYTE, &status), "MPI_File_write_ordered", &status);
    }
    else
    {
      char b[1];
      MPI_Status status;
      MPIUtility::handleReturnValue(MPI_File_write_ordered(fileHandle, b, 0, MPI_BYTE, &status), "MPI_File_write_ordered", &status);
    }

    // write log line for own rank
    MPIUtility::handleReturnValue(MPI_File_write_ordered(fileHandle, data.str().c_str(), data.str().length(), MPI_BYTE, MPI_STATUS_IGNORE), "MPI_File_write_ordered");

    // close file
    MPIUtility::handleReturnValue(MPI_File_close(&fileHandle), "MPI_File_close");
  }
  else  // standard POSIX output
  {
    // open log file
    std::ofstream file;
    OutputWriter::Generic::openFile(file, filename.str(), true);

    if (outputHeader)
    {
      file << header.str();
    }

    // write data to file
    file << data.str();
    file.close();
  }
}

}
