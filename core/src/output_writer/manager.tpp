#include "output_writer/manager.h"

#include "output_writer/callback.h"
#include "output_writer/paraview.h"
#include "output_writer/python.h"
#include "output_writer/exfile.h"

namespace OutputWriter
{
  
template<typename DataType>
void Manager::writeOutput(DataType &problemData, int timeStepNo, double currentTime) const
{
  for(auto &outputWriter : this->outputWriter_)
  {
    LOG(DEBUG) << "select outputWriter";
    if (std::dynamic_pointer_cast<Callback>(outputWriter) != nullptr)
    {
      std::shared_ptr<Callback> writer = std::static_pointer_cast<Callback>(outputWriter);
      writer->write<DataType>(problemData, timeStepNo, currentTime);
    }
    else if (std::dynamic_pointer_cast<Exfile>(outputWriter) != nullptr)
    {
      std::shared_ptr<Exfile> writer = std::static_pointer_cast<Exfile>(outputWriter);
      writer->write<DataType>(problemData, timeStepNo, currentTime);
    }
    else if (std::dynamic_pointer_cast<Paraview>(outputWriter) != nullptr)
    {
      std::shared_ptr<Paraview> writer = std::static_pointer_cast<Paraview>(outputWriter);
      writer->write<DataType>(problemData, timeStepNo, currentTime);
    }
    else if (std::dynamic_pointer_cast<Python>(outputWriter) != nullptr)
    {
      std::shared_ptr<Python> writer = std::static_pointer_cast<Python>(outputWriter);
      writer->write<DataType>(problemData, timeStepNo, currentTime);
    }
  }
}
  
};   // namespace