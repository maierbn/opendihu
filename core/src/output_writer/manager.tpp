#include "output_writer/manager.h"

namespace OutputWriter
{
  
template<typename DataType>
void Manager::writeOutput(DataType &problemData, int timeStepNo, double currentTime) const
{
  for(auto &outputWriter : this->outputWriter_)
  {
    outputWriter->write<DataType>(problemData, timeStepNo, currentTime);
  }
}
  
};   // namespace