#include "control/performance_measurement.h"

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

std::map<std::string, PerformanceMeasurement::Measurement> PerformanceMeasurement::measurements_;
std::map<std::string,std::string> PerformanceMeasurement::parameters_;

PerformanceMeasurement::Measurement::Measurement() :
  start(0.0), totalDuration(0.0), nTimeSpans(0), totalError(0.0), nErrors(0)
{
}

void PerformanceMeasurement::start(std::string name)
{
  std::map<std::string, Measurement>::iterator iter = measurements_.find(name);

  // if there is no entry of name yet, create new
  if (iter == measurements_.end())
  {
    auto insertedIter = measurements_.insert(std::pair<std::string, Measurement>(name, Measurement()));
    iter = insertedIter.first;
  }

  // measure current time
  iter->second.start = MPI_Wtime();
}

void PerformanceMeasurement::stop(std::string name, int numberAccumulated)
{
  double stopTime = MPI_Wtime();

  if (measurements_.find(name) == measurements_.end())
  {
    LOG(ERROR) << "PerformanceMeasurement stop with name \"" << name << "\", a corresponding start is not present.";
  }
  else
  {
    Measurement &measurement = measurements_[name];
    double duration = stopTime - measurement.start;
    measurement.totalDuration += duration;
    measurement.nTimeSpans += numberAccumulated;

    VLOG(2) << "PerformanceMeasurement::stop(" << name << "), time span [" << measurement.start << "," << stopTime << "], duration=" << duration
      << ", now total: " << measurement.totalDuration << ", nTimeSpans: " << measurement.nTimeSpans;
  }
}

void PerformanceMeasurement::writeLogFile(std::string logFileName)
{
  //LOG(DEBUG) << "PerformanceMeasurement::writeLogFile \"" << logFileName;

  parseStatusInformation();

  const bool combined = true;   /// if the output is using MPI Output

  int ownRankNo = DihuContext::partitionManager()->rankNoCommWorld();

  // determine file name
  std::stringstream filename;
  filename << logFileName;
  if (!combined)
    filename << "." << std::setw(7) << std::setfill('0') << ownRankNo;
  filename << ".csv";
  logFileName = filename.str();

  // compose header
  std::stringstream header;
  header << "# timestamp;hostname;version;";

  // write parameter names
  for (std::pair<std::string,std::string> parameter : parameters_)
  {
    header << parameter.first << ";";
  }

  // write measurement names
  for (std::pair<std::string, Measurement> measurement : measurements_)
  {
    header << measurement.first << ";n;";
  }
  header << std::endl;

  // compose data
  std::stringstream data;

  // time stamp
  auto t = std::time(nullptr);
  auto tm = *std::localtime(&t);
  data << std::put_time(&tm, "%Y/%m/%d %H:%M:%S") << ";";

  // host name
  char hostname[MAXHOSTNAMELEN+1];
  gethostname(hostname, MAXHOSTNAMELEN+1);
  data << std::string(hostname) << ";";
  data << DihuContext::versionText() << ";";

  // write parameters
  for (std::pair<std::string,std::string> parameter : parameters_)
  {
    data << parameter.second << ";";
  }

  // write measurement values
  for (std::pair<std::string, Measurement> measurement : measurements_)
  {
    data << measurement.second.totalDuration << ";"
    << measurement.second.nTimeSpans << ";";
  }
  data << std::endl;

  // check if header has to be added to file
  bool outputHeader = true;
  if (!combined || ownRankNo == 0)
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
  if (combined)   // MPI output
  {
    // open log file to create directory
    std::ofstream file;
    if (ownRankNo == 0)
    {
      file = OutputWriter::Generic::openFile(filename.str(), true);
      file.close();
    }

    LOG(DEBUG) << "open MPI file \"" << logFileName << "\".";


    // open file
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_File fileHandle;
    MPIUtility::handleReturnValue(MPI_File_open(MPI_COMM_WORLD, logFileName.c_str(),
                                                //MPI_MODE_WRONLY | MPI_MODE_CREATE | MPI_MODE_UNIQUE_OPEN,
                                                //MPI_MODE_WRONLY | MPI_MODE_CREATE | MPI_MODE_APPEND,
                                                MPI_MODE_WRONLY | MPI_MODE_CREATE,
                                                MPI_INFO_NULL, &fileHandle), "MPI_File_open");

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

    MPIUtility::handleReturnValue(MPI_File_write_ordered(fileHandle, data.str().c_str(), data.str().length(), MPI_BYTE, MPI_STATUS_IGNORE), "MPI_File_write_ordered");

    MPIUtility::handleReturnValue(MPI_File_close(&fileHandle), "MPI_File_close");
  }
  else  // standard POSIX output
  {
    // open log file
    std::ofstream file = OutputWriter::Generic::openFile(filename.str(), true);

    if (outputHeader)
    {
      file << header.str();
    }

    // write data to file
    file << data.str();
    file.close();
  }
}

template<>
void PerformanceMeasurement::measureError<double>(std::string name, double differenceVector)
{
  std::map<std::string, Measurement>::iterator iter = measurements_.find(name);

  // if there is no entry of name yet, create new
  if (iter == measurements_.end())
  {
    auto insertedIter = measurements_.insert(std::pair<std::string, Measurement>(name, Measurement()));
    iter = insertedIter.first;
  }

  iter->second.totalError += fabs(differenceVector);
  iter->second.nErrors++;
}

void PerformanceMeasurement::parseStatusInformation()
{
   // adapted from https://stackoverflow.com/questions/669438/how-to-get-memory-usage-at-runtime-using-c

   std::string pid, comm, state, ppid, pgrp, session, ttyNr, tpgid, flags, minflt, cminflt, majflt, cmajflt, stime, cutime, cstime,
     priority, nice, o, itrealvalue, starttime, size, resident, shared, text, none;

   std::ifstream stream("/proc/self/stat", std::ios_base::in);

   unsigned long vsize;
   long long rss;
   long long utime;
   long long data;

   // parse entries in /proc/self/stat
   stream >> pid >> comm >> state >> ppid >> pgrp >> session >> ttyNr
               >> tpgid >> flags >> minflt >> cminflt >> majflt >> cmajflt
               >> utime >> stime >> cutime >> cstime >> priority >> nice
               >> o >> itrealvalue >> starttime >> vsize >> rss;

   stream.close();
   stream.open("/proc/self/statm", std::ios_base::in);

   // parse entries in /proc/self/statm
   stream >> size >> resident >> shared >> text >> none >> data;
   stream.close();

   // compute final values
   int pageSize = getpagesize(); // page size in KB
   long long virtualMemorySize = vsize;
   long long residentSetSize = rss * pageSize;
   long long dataSize = data * pageSize;

   double totalUsertime = (double)utime / sysconf(_SC_CLK_TCK);

   // http://man7.org/linux/man-pages/man5/proc.5.html
   /** Resident Set Size: number of pages the process has
                        in real memory.  This is just the pages which count
                        toward text, data, or stack space.  This does not
                        include pages which have not been demand-loaded in,
                        or which are swapped out.
   */

   // store parameters
   PerformanceMeasurement::setParameter("memoryPage", pageSize);
   PerformanceMeasurement::setParameter("memoryVirtual", virtualMemorySize);
   PerformanceMeasurement::setParameter("memoryResidentSet", residentSetSize);
   PerformanceMeasurement::setParameter("memoryData", dataSize);
   PerformanceMeasurement::setParameter("totalUsertime", totalUsertime);

   // format total user time in readable format for output
   int minTotal = int(totalUsertime/60);
   int hTotal = int(totalUsertime/3600);
   int d = int(totalUsertime/(3600*24));
   int h = hTotal - d*24;
   int min = minTotal - h*60 - d*24*60;
   double sDouble = totalUsertime - minTotal*60;
   int s = int(sDouble);
   sDouble -= s;

   std::stringstream message;
   if (d != 0)
     message << d << "d ";
   if (h != 0)
     message << h << ":";
   if (min != 0)
   {
     std::stringstream secondsFraction;
     secondsFraction << std::setprecision(4) << sDouble;
     if (h == 0)
       message << min;
     else
       message << std::setw(2) << std::setfill('0') << min;

     message << ":" << std::setw(2) << std::setfill('0') << s << secondsFraction.str().substr(1);
     if (h == 0)
       message << " min";
   }
   else
   {
     message << s << "s";
   }
   LOG(INFO) << "Total user time: " << message.str();
}

} // namespace
