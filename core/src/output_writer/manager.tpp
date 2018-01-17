#include "output_writer/manager.h"

#include "output_writer/callback.h"
#include "output_writer/paraview.h"
#include "output_writer/python.h"

namespace OutputWriter
{
  
template<typename DataType>
void Manager::writeOutput(DataType &problemData, int timeStepNo, double currentTime) const
{
  for(auto &outputWriter : this->outputWriter_)
  {
    if (std::dynamic_pointer_cast<Callback>(outputWriter) != nullptr)
    {
      std::shared_ptr<Callback> writer = std::static_pointer_cast<Callback>(outputWriter);
      writer->write<DataType>(problemData, timeStepNo, currentTime);
    }
  }
}
  
};   // namespace