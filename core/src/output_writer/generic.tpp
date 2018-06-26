#include "output_writer/generic.h"

#include <iomanip>

#include "easylogging++.h"
#include "utility/python_utility.h"

namespace OutputWriter
{


template<typename DataType>
bool Generic::prepareWrite(DataType& data, int timeStepNo, double currentTime)
{
  VLOG(1) << "Generic::prepareWrite timeStepNo="<<timeStepNo<<", currentTime="<<currentTime;

  if (!data.mesh())
  {
    LOG(FATAL) << "Mesh is not set!";
  }

  timeStepNo_ = timeStepNo;
  currentTime_ = currentTime;
  int outputInterval = PythonUtility::getOptionInt(specificSettings_, "outputInterval", 1);

  int oldWriteCallCount = writeCallCount_;
  writeCallCount_++;

  VLOG(1) << " Generic::prepareWrite, writeCallCount_=" << writeCallCount_ << ", outputInterval: " << outputInterval;
  
  // if no output should be written, because of interval, return false
  if (oldWriteCallCount % outputInterval != 0)
  {
    VLOG(1) << " do not write";
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
    //s << "_" << std::setw(5) << std::setfill('0') << timeStepNo;
    s << "_" << std::setw(7) << std::setfill('0') << outputFileNo_;   // use a continuous counter for the output file 
  }
  outputFileNo_++;
  
  filename_ = s.str();
  return true;
}

}; // namespace