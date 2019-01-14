#include "output_writer/megamol/megamol.h"

#include <fstream>
#include <iostream>
#include <thread>
#include <chrono>

#include "easylogging++.h"
#include "base64.h"

#include "output_writer/megamol/loop_output.h"

namespace OutputWriter
{

template<typename DataType>
void MegaMol::write(DataType& data, int timeStepNo, double currentTime)
{
  // check if output should be written in this timestep and prepare filename
  if (!Generic::prepareWrite(data, timeStepNo, currentTime))
  {
    return;
  }

  // get io object of adios from context
  std::shared_ptr<adios2::IO> io = context_.adiosIo();

  try
  {
    // create ADIOS writer, if it does not yet exist
    if (!writer_)
    {
      assert(io);
      adios2::Engine writer = io->Open(this->filenameBase_, adios2::Mode::Write);
      writer_ = std::make_shared<adios2::Engine>(writer);
    }

    // begin output to file for current time step
    writer_->BeginStep();

    // collect all available meshes
    std::set<std::string> meshNames;
    LoopOverTuple::loopCollectMeshNames<typename DataType::OutputFieldVariables>(data.getOutputFieldVariables(), meshNames);

    // loop over meshes
    for (std::string meshName : meshNames)
    {
      // loop over all field variables and output those that are associated with the mesh given by meshName
      MegaMolLoopOverTuple::loopOutput(data.getOutputFieldVariables(), data.getOutputFieldVariables(), meshName, specificSettings_, writer_, io);
    }

    // end output for current time step
    writer_->EndStep();
  }
  catch (std::invalid_argument &e)
  {
    LOG(ERROR) << "Invalid argument exception";
    LOG(ERROR) << e.what();
  }
  catch (std::ios_base::failure &e)
  {
    LOG(ERROR) << "System exception";
    LOG(ERROR) << e.what();
  }
  catch (std::exception &e)
  {
    LOG(ERROR) << "Exception";
    LOG(ERROR) << e.what();
  }

}

} // namespace
