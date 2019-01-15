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
#ifdef HAVE_ADIOS

  // get io object of adios from context
  std::shared_ptr<adios2::IO> adiosIo = context_.adiosIo();

  try
  {
    // create ADIOS writer, if it does not yet exist
    if (!writer_)
    {
      assert(adiosIo);
      adios2::Engine writer = adiosIo->Open(this->filenameBase_, adios2::Mode::Write);
      writer_ = std::make_shared<adios2::Engine>(writer);
    }

    // begin output to file for current time step
    writer_->BeginStep();

    // create box variable
    if (!boxVariable_)
    {
      boxVariable_ = std::make_shared<adios2::Variable<double>>(adiosIo->DefineVariable<double>("box", {6}, {0}, {6}));
    }

    // collect all available meshes
    std::set<std::string> meshNames;
    LoopOverTuple::loopCollectMeshNames<typename DataType::OutputFieldVariables>(data.getOutputFieldVariables(), meshNames);

    BoundingBox localBoundingBox, globalBoundingBox;
    localBoundingBox.initialized = false;

    // loop over meshes
    for (std::string meshName : meshNames)
    {
      // loop over all field variables and output those that are associated with the mesh given by meshName
      MegaMolLoopOverTuple::loopOutput(data.getOutputFieldVariables(), data.getOutputFieldVariables(), meshName, specificSettings_, writer_, adiosIo, localBoundingBox);
    }


    LOG(DEBUG) << "local min: " << localBoundingBox.min << ", local max: " << localBoundingBox.max;
    MPI_Reduce(localBoundingBox.min.data(), globalBoundingBox.min.data(), 3, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
    MPI_Reduce(localBoundingBox.max.data(), globalBoundingBox.max.data(), 3, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    // note: this assumes that all ranks execute the reduction calls, i.e. all have these type of meshes

    LOG(DEBUG) << "min: " << localBoundingBox.min << ", globalBoundingBox.min: " << globalBoundingBox.min
      << ", max: " << localBoundingBox.max << ", globalBoundingBox.max: " << globalBoundingBox.max;

    // define bounding box variable
    if (context_.ownRankNo() == 0)
    {
      std::array<double,6> minMaxValues({
        globalBoundingBox.min[0], globalBoundingBox.max[0],
        globalBoundingBox.min[1], globalBoundingBox.max[1],
        globalBoundingBox.min[2], globalBoundingBox.max[2]
      });
      LOG(DEBUG) << "minMaxValues: " << minMaxValues;
      assert(boxVariable_);
      assert(boxVariable_.get());
      writer_->Put<double>(*boxVariable_.get(), minMaxValues.data());
    }

    // end output for current time step
    writer_->EndStep();

    int timeStepCloseInterval = specificSettings_.getOptionInt("timeStepCloseInterval",1);

    LOG(DEBUG) << "timeStepCloseInterval: " << timeStepCloseInterval << ", writeCallCount_: " << writeCallCount_;

    if (writeCallCount_ % timeStepCloseInterval == 0 && writeCallCount_ > 0)
    {
      writer_->Close();
      writer_ = nullptr;
    }
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
#endif
}

} // namespace
