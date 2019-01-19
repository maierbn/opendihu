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
  std::shared_ptr<Partition::RankSubset> rankSubset = context_.partitionManager()->rankSubsetForCollectiveOperations();

  try
  {
    // create ADIOS writer, if it does not yet exist
    if (!writer_)
    {
      assert(adiosIo);
      adios2::Engine writer = adiosIo->Open(this->filenameBase_, adios2::Mode::Write, rankSubset->mpiCommunicator());
      writer_ = std::make_shared<adios2::Engine>(writer);

      // write python script and meta data as attribute
      adiosIo->DefineAttribute<std::string>("config", this->context_.pythonScriptText());
      adiosIo->DefineAttribute<std::string>("version", this->context_.versionText());
      adiosIo->DefineAttribute<std::string>("meta", this->context_.metaText());
    }

    // begin output to file for current time step
    writer_->BeginStep();

    // collect all available meshes
    std::set<std::string> meshNames;
    LoopOverTuple::loopCollectMeshNames<typename DataType::OutputFieldVariables>(data.getOutputFieldVariables(), meshNames);

    BoundingBox localBoundingBox, globalBoundingBox;
    localBoundingBox.initialized = false;
    std::map<std::string,int> nNodesGlobalAllMeshes;

    // loop over meshes
    for (std::string meshName : meshNames)
    {
      // loop over all field variables and output those that are associated with the mesh given by meshName
      MegaMolLoopOverTuple::loopOutput(data.getOutputFieldVariables(), data.getOutputFieldVariables(), meshName, specificSettings_, writer_, adiosIo, localBoundingBox, nNodesGlobalAllMeshes);
    }

    // reduce the bounding box values
    MPI_Reduce(localBoundingBox.min.data(), globalBoundingBox.min.data(), 3, MPI_DOUBLE, MPI_MIN, 0, rankSubset->mpiCommunicator());
    MPI_Reduce(localBoundingBox.max.data(), globalBoundingBox.max.data(), 3, MPI_DOUBLE, MPI_MAX, 0, rankSubset->mpiCommunicator());

    LOG(DEBUG) << "reduced bounding box: " << globalBoundingBox.min << ", " << globalBoundingBox.max;
    // note: this assumes that all ranks execute the reduction calls, i.e. all have these type of meshes

    // write bounding box
    if (rankSubset->ownRankNo() == 0)
    {
      if (!boxVariable_)
      {
        boxVariable_ = std::make_shared<adios2::Variable<double>>(adiosIo->DefineVariable<double>("box", {6},{0},{6}, adios2::ConstantDims));
      }
      // megamol box: (xmin, ymin, zmax, xmax, ymax, zmin)
      // opendihu box: (xmin, ymin, zmax, xmax, ymax, zmin)

      boundingBoxValues_.resize(6);
      boundingBoxValues_[0] = globalBoundingBox.min[0];
      boundingBoxValues_[1] = globalBoundingBox.min[1];
      boundingBoxValues_[2] = globalBoundingBox.max[2];
      boundingBoxValues_[3] = globalBoundingBox.max[0];
      boundingBoxValues_[4] = globalBoundingBox.max[1];
      boundingBoxValues_[5] = globalBoundingBox.min[2];
      LOG(DEBUG) << "boundingBoxValues: " << boundingBoxValues_;

      assert(boxVariable_);
      assert(boxVariable_.get());

      writer_->Put<double>(*boxVariable_.get(), boundingBoxValues_.data());
    }

    // write number of global nodes as variable p_count

    // determine global count
    int nNodesGlobalAllMeshesSumLocal = 0;
    LOG(DEBUG) << "nNodesGlobalAllMeshes: " << nNodesGlobalAllMeshes;
    for (std::map<std::string,int>::iterator iter = nNodesGlobalAllMeshes.begin(); iter != nNodesGlobalAllMeshes.end(); iter++)
    {
      nNodesGlobalAllMeshesSumLocal += iter->second;
    }

    globalNumberOfNodes_ = 0;
    MPI_Reduce(&nNodesGlobalAllMeshesSumLocal, &globalNumberOfNodes_, 1, MPI_INT, MPI_SUM, 0, rankSubset->mpiCommunicator());

    LOG(DEBUG) << "nNodesGlobalAllMeshesSumLocal: " << nNodesGlobalAllMeshesSumLocal << ", globalNumberOfNodes_: " << globalNumberOfNodes_;

    // write global number of nodes
    if (rankSubset->ownRankNo() == 0)
    {
      if (!pCountVariable_)
      {
        pCountVariable_ = std::make_shared<adios2::Variable<int>>(adiosIo->DefineVariable<int>("p_count"));
      }

      assert(*pCountVariable_.get());
      writer_->Put<int>(*pCountVariable_.get(), &globalNumberOfNodes_);
    }

    // end output for current time step
    writer_->EndStep();

    int timeStepCloseInterval = specificSettings_.getOptionInt("timeStepCloseInterval",1);

    LOG(DEBUG) << "timeStepCloseInterval: " << timeStepCloseInterval << ", writeCallCount_: " << writeCallCount_;

    if (writeCallCount_ % timeStepCloseInterval == 0 && writeCallCount_ > 0)
    {
      writer_->Close();
      //writer_ = nullptr;
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
