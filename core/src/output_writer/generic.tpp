#include "output_writer/generic.h"

#include <iomanip>

#include "easylogging++.h"
#include "utility/python_utility.h"

namespace OutputWriter
{

 
template<typename DataType>
bool Generic::prepareWrite(DataType& data, int timeStepNo, double currentTime)
{
  if (!data.mesh())
  {
    LOG(FATAL) << "Mesh is not set!";
  }
  
  int outputInterval = PythonUtility::getOptionInt(specificSettings_, "outputInterval", 1);
  
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
  filename_ = s.str();
  return true;
}

}; // namespace