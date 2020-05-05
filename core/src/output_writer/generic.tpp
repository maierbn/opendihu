#include "output_writer/generic.h"

#include <iomanip>
#include <cmath>
#include <iomanip>

#include "easylogging++.h"
#include "utility/python_utility.h"

namespace OutputWriter
{

template<typename DataType>
bool Generic::prepareWrite(DataType& data, int timeStepNo, double currentTime, int callCountIncrement)
{
  VLOG(2) << "Generic::prepareWrite (\"" << filenameBase_ << "\") timeStepNo=" << timeStepNo << ", currentTime=" << currentTime;

  if (!data.functionSpace())
  {
    LOG(FATAL) << "Function space is not set!";
  }

  // update timestep no and current time
  timeStepNo_ = timeStepNo;
  currentTime_ = currentTime;

  // increment the counter how often the output writer was called
  int oldWriteCallCount = writeCallCount_;
  writeCallCount_ += callCountIncrement;

  VLOG(2) << " Generic::prepareWrite (\"" << filenameBase_ << "\"), writeCallCount_=" << writeCallCount_ << ", outputInterval: " << outputInterval_;

  // if no output should be written, because of interval, return false
  int lastCallCountToWrite = writeCallCount_ - (writeCallCount_%outputInterval_);

  // if there was a call count between the last and the current where the file should have been written
  if (oldWriteCallCount < lastCallCountToWrite && lastCallCountToWrite <= writeCallCount_)
  {
    // file should be written now

    // add time step number to file name base
    std::stringstream s;
    s << filenameBase_;
    if (timeStepNo != -1)
    {
      switch (fileNumbering_) {
      case file_numbering_incremental:
        s << "_" << std::setw(7) << std::setfill('0') << outputFileNo_;   // use a continuous counter for the output file
        break;
      case file_numbering_by_time_step_index:
        s << "_" << std::setw(7) << std::setfill('0') << (writeCallCount_-1);   // 0 based: first call corresponts to 0
        break;
      default:
        LOG(ERROR) << "BUG: Unknown file numbering '" << fileNumbering_ << "'. This should not happen.";
      }
    }
    outputFileNo_++;
    filenameBaseWithNo_ = s.str();

    // add rank no to file name base
    if (data.functionSpace()->meshPartition()->nRanks() > 1)
    {
      appendRankNo(s, data.functionSpace()->meshPartition()->nRanks(), data.functionSpace()->meshPartition()->ownRankNo());
    }

    // set the current filename and return true, which indicates that the output writer will write the file
    filename_ = s.str();
    return true;
  }

  // file should not be written now
  VLOG(2) << " do not write (\"" << filenameBase_ << "\")";
  return false;
}

}  // namespace
