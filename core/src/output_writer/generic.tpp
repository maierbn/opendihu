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
  LOG(DEBUG) << "Generic::prepareWrite timeStepNo=" <<timeStepNo<< ", currentTime=" << currentTime;

  if (!data.mesh())
  {
    LOG(FATAL) << "Mesh is not set!";
  }


  timeStepNo_ = timeStepNo;
  currentTime_ = currentTime;
  int outputInterval = PythonUtility::getOptionInt(specificSettings_, "outputInterval", 1, PythonUtility::Positive);

  int oldWriteCallCount = writeCallCount_;
  writeCallCount_++;

  // if no output should be written, because of interval, return false
  if (oldWriteCallCount % outputInterval != 0)
    return false;

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
    s << "_" << std::setw(5) << std::setfill('0') << timeStepNo;
  }

  // add rank no to file name base
  if (data.mesh()->meshPartition()->nRanks() > 1)
  {
    appendRankNo(s, data.mesh()->meshPartition()->nRanks(), data.mesh()->meshPartition()->ownRankNo());
  }

  filename_ = s.str();
  return true;
}

}; // namespace
