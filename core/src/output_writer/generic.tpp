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

  int oldWriteCallCount = writeCallCount_;
  writeCallCount_++;

  VLOG(2) << " Generic::prepareWrite, writeCallCount_=" << writeCallCount_ << ", outputInterval: " << outputInterval_;
  
  // if no output should be written, because of interval, return false
  if (oldWriteCallCount % outputInterval_ != 0)
  {
    VLOG(2) << " do not write";
    return false;
  }


  // add time step number to file name base
  std::stringstream s;
  s << filenameBase_;
  if (timeStepNo != -1)
  {
    s << "_" << std::setw(7) << std::setfill('0') << outputFileNo_;   // use a continuous counter for the output file
  }
  outputFileNo_++;
  filenameBaseWithNo_ = s.str();
  
  // add rank no to file name base
  if (data.functionSpace()->meshPartition()->nRanks() > 1)
  {
    appendRankNo(s, data.functionSpace()->meshPartition()->nRanks(), data.functionSpace()->meshPartition()->ownRankNo());
  }

  filename_ = s.str();
  return true;
}

}  // namespace
