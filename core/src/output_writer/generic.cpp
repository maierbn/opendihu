#include "output_writer/generic.h"

#include <iomanip>

#include "easylogging++.h"
#include "control/python_utility.h"

namespace OutputWriter
{
 
void Generic::write(Data::Data& data, PyObject* specificSettings, int timeStepNo, double currentTime)
{
  LOG(DEBUG) << "Generic::write";
  
  int outputFrequency = PythonUtility::getOptionInt(specificSettings, "frequency", 1);
  
  int oldWriteCallCount = writeCallCount_;
  writeCallCount_++;
  
  // if no output should be written, because of frequency, return
  if (oldWriteCallCount % outputFrequency != 0)
    return;
    
  // determine filename base
  if (filenameBase_.empty())
  {
    filenameBase_ = PythonUtility::getOptionString(specificSettings, "filename", "out");
  }
  
  // add time step number to file name base
  std::stringstream s;
  s << filenameBase_;
  if (timeStepNo != -1)
  {
    s << "_" << std::setw(5) << std::setfill('0') << timeStepNo;
  }
  filename_ = s.str();
  
  // if data.mesh() is of class Mesh::RegularFixed<D>
  std::shared_ptr<Mesh::Mesh> mesh = data.mesh();
  if (!mesh)
  {
    LOG(FATAL) << "mesh is not set";
  }
  
  writeSolution(data);
}
  
};