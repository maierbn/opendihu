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

  /*

  Example for step = .2 interval = 2
  ----------------------------------
  start time      .0    .0    .2    .4    .6
  time step       .0    .2    .4    .6    .8
  increment       0     1     1     1     1
  old call count  0     0     1     2     3
  call count      0     1     2     3     4

  last            0     0     2     2     4
  write           W     -     w     -     w


  Example for step = .1 interval = 4
  ----------------------------------
  start time      .0 .0 .1 .2 .3 .4 .5 .6 .7
  time step       .0 .1 .2 .3 .4 .5 .6 .7 .8
  increment       0  1  1  1  1  1  1  1  1
  old call count  0  0  1  2  3  4  5  6  7
  call count      0  1  2  3  4  5  6  7  8

  last            0  0  0  0  4  4  4  4  8
  write           W  -  -  -  w  -  -  -  w


  The first column is an (optional) call to write the initial values with callCountIncrement=0.
  If the call is not made, the remaining table does not change.
  '-' means no output written
  'w' means output written
  'W' means forced output for initial values

  */


  // if no output should be written, because of interval, return false
  int lastCallCountToWrite = writeCallCount_ - (writeCallCount_%outputInterval_);

  // if there was a call count between the last and the current where the file should have been written
  if ((oldWriteCallCount < lastCallCountToWrite && lastCallCountToWrite <= writeCallCount_) || callCountIncrement == 0)
  {
    // file should be written

    // add time step number to file name base
    std::stringstream s;
    s << filenameBase_;
    if (timeStepNo != -1)
    {
      switch (fileNumbering_)
      {
      case fileNumberingIncremental:
        // use a continuous counter for the output file
        s << "_" << std::setw(7) << std::setfill('0') << outputFileNo_;
        break;
      case fileNumberingByTimeStepIndex:
        // '0' corresponds to initial data (if callCountIncrement was 0)
        // '1',...,n to the normal data.
        //   For time integrators it corresponds to the data after the first,...,n-th time step
        s << "_" << std::setw(7) << std::setfill('0') << writeCallCount_;
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
  else
  {
    // file should not be written
    VLOG(2) << " do not write (\"" << filenameBase_ << "\")";
    return false;
  }
}

}  // namespace
