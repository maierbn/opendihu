#include "output_writer/manager.h"

#include "output_writer/python_callback/python_callback.h"
#include "output_writer/python_file/python_file.h"
#include "output_writer/paraview/paraview.h"
#include "output_writer/exfile/exfile.h"
#include "output_writer/megamol/megamol.h"
#include "control/performance_measurement.h"

namespace OutputWriter
{

template<typename DataType>
void Manager::writeOutput(DataType &problemData, int timeStepNo, double currentTime) const
{
  // start duration measurement
  Control::PerformanceMeasurement::start("write output");

  for (auto &outputWriter : this->outputWriter_)
  {
    if (std::dynamic_pointer_cast<Exfile>(outputWriter) != nullptr)
    {
      std::shared_ptr<Exfile> writer = std::static_pointer_cast<Exfile>(outputWriter);
      writer->write<DataType>(problemData, timeStepNo, currentTime);
    }
    else if (std::dynamic_pointer_cast<Paraview>(outputWriter) != nullptr)
    {
      std::shared_ptr<Paraview> writer = std::static_pointer_cast<Paraview>(outputWriter);
      writer->write<DataType>(problemData, timeStepNo, currentTime);
    }
    else if (std::dynamic_pointer_cast<PythonCallback>(outputWriter) != nullptr)
    {
      std::shared_ptr<PythonCallback> writer = std::static_pointer_cast<PythonCallback>(outputWriter);
      writer->write<DataType>(problemData, timeStepNo, currentTime);
    }
    else if (std::dynamic_pointer_cast<PythonFile>(outputWriter) != nullptr)
    {
      std::shared_ptr<PythonFile> writer = std::static_pointer_cast<PythonFile>(outputWriter);
      writer->write<DataType>(problemData, timeStepNo, currentTime);
    }
    else if (std::dynamic_pointer_cast<MegaMol>(outputWriter) != nullptr)
    {
      std::shared_ptr<MegaMol> writer = std::static_pointer_cast<MegaMol>(outputWriter);
      writer->write<DataType>(problemData, timeStepNo, currentTime);
    }
  }

  // stop duration measurement
  Control::PerformanceMeasurement::stop("write output");
}

}  // namespace
