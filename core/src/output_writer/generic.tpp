#include "output_writer/generic.h"

#include <iomanip>
#include <cmath>
#include <iomanip>

#include "easylogging++.h"
#include "utility/python_utility.h"

namespace OutputWriter
{

template<typename DataType>
bool Generic::prepareWrite(DataType& data, int timeStepNo, double currentTime)
{
  VLOG(2) << "Generic::prepareWrite timeStepNo=" << timeStepNo << ", currentTime=" << currentTime;

  if (!data.functionSpace())
  {
    LOG(FATAL) << "Function space is not set!";
  }

  timeStepNo_ = timeStepNo;
  currentTime_ = currentTime;
  int outputInterval = PythonUtility::getOptionInt(specificSettings_, "outputInterval", 1, PythonUtility::Positive);

  int oldWriteCallCount = writeCallCount_; // I think this one is not necessary. (Especially with the change underneath.) -- Aaron
  // @Aaron: I reverted your changes
  // example: outputInterval=10
  // output at numbers
  // Benni's version: 0, 10, 20, 30
  //  Aaron's version: 9, 19, 29, 39
  // We need output at the first time (0) for examples that don't have timestepping, e.g. Laplace eq., please reevaluate, what exactly is needed for easier times.
  
  writeCallCount_++;

  VLOG(2) << " Generic::prepareWrite, writeCallCount_=" << writeCallCount_ << ", outputInterval: " << outputInterval;
  
  // if no output should be written, because of interval, return false
  if (oldWriteCallCount % outputInterval != 0) // I changed the condition to '<--this' from '(oldWriteCallCount % outputInterval != 0)'. Now it is easyer to have data output at regular times. -- Aaron 
  {
    VLOG(2) << " do not write";
    return false;
  }

  // determine filename base
  if (filenameBase_.empty()
    && PythonUtility::getOptionString(specificSettings_, "format", "Callback") != "Callback")
  {
    filenameBase_ = PythonUtility::getOptionString(specificSettings_, "filename", "out");
  }

  // add time step number to file name base
  std::stringstream s;
  s << filenameBase_;
  if (timeStepNo != -1)
  {
    s << "_" << std::setw(7) << std::setfill('0') << outputFileNo_;   // use a continuous counter for the output file 
  }
  outputFileNo_++;
  
  // add rank no to file name base
  if (data.functionSpace()->meshPartition()->nRanks() > 1)
  {
    appendRankNo(s, data.functionSpace()->meshPartition()->nRanks(), data.functionSpace()->meshPartition()->ownRankNo());
  }

  filename_ = s.str();
  return true;
}

}; // namespace
