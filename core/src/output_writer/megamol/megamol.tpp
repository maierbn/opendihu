#include "output_writer/megamol/megamol.h"

#include <fstream>
#include <iostream>
#include <thread>
#include <chrono>

#include "easylogging++.h"
#include "base64.h"

#include "output_writer/megamol/loop_output.h"

#ifdef HAVE_MEGAMOL
#include <libzmq/zmq.hpp>
#endif

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
    if (!adiosWriters_[currentOpenWriterIndex_])
    {
      adiosWriters_[currentOpenWriterIndex_] = std::make_shared<adios_writer_t>();
      adiosWriters_[currentOpenWriterIndex_]->nOpenWriters = 0;
    }
    std::shared_ptr<adios_writer_t> writer = adiosWriters_[currentOpenWriterIndex_];
    std::vector<double> geometryFieldScalarValues;

    // if writer is not open, open it
    if (writer->nOpenWriters == 0)
    {
      assert(adiosIo);

      // determine current filename
      std::stringstream outputFileName;
      if (useFrontBackBuffer_)
      {
        outputFileName << this->filenameBase_ << "_" << currentOpenWriterIndex_;
      }
      else
      {
        outputFileName << filenameBaseWithNo_;
      }

      currentFilename_ = outputFileName.str();
      LOG(DEBUG) << "create new engine, on outputFileName: \"" << outputFileName.str() << "\"";

      // create new writer
      adios2::Engine engine = adiosIo->Open(outputFileName.str(), adios2::Mode::Write, rankSubset->mpiCommunicator());
      writer->engine = std::make_shared<adios2::Engine>(engine);

      writer->megaMolWriterContext.approximateDistanceBetweenFibers = -1;  // initialize to not set

      // write python script and meta data as attribute
      //adiosIo->DefineAttribute<std::string>("config", this->context_.pythonScriptText());
      //adiosIo->DefineAttribute<std::string>("version", this->context_.versionText());
      //adiosIo->DefineAttribute<std::string>("meta", this->context_.metaText());

      // begin output to file for current time step
      writer->engine->BeginStep();
    }
    writer->nOpenWriters++;
    LOG(DEBUG) << "in MegaMol::write, opening time step, currentOpenWriterIndex_: " << currentOpenWriterIndex_ << ", nOpenWriters: " << writer->nOpenWriters;

    // collect all available meshes
    std::set<std::string> meshNames;
    LoopOverTuple::loopCollectMeshNames<typename DataType::OutputFieldVariables>(data.getOutputFieldVariables(), meshNames);

    // loop over meshes and retrieve data
    for (std::string meshName : meshNames)
    {
      // loop over all field variables and output those that are associated with the mesh given by meshName
      MegaMolLoopOverTuple::loopOutput(data.getOutputFieldVariables(), data.getOutputFieldVariables(), meshName, specificSettings_,
                                       writer->megaMolWriterContext);
    }

    // only write the variables to file when all collecting calls have been made
    if (writer->nOpenWriters == combineNInstances_)
    {
      std::vector<Vec3> &geometryFieldValues = writer->megaMolWriterContext.geometryFieldValues;
      std::vector<double> &scalarFieldVariableValues = writer->megaMolWriterContext.scalarFieldVariableValues;
      double approximateDistanceBetweenFibers = writer->megaMolWriterContext.approximateDistanceBetweenFibers;

      // compute local bounding box
      BoundingBox localBoundingBox;
      localBoundingBox.min = geometryFieldValues[0];
      localBoundingBox.max = geometryFieldValues[0];

      // loop over geometry field entries
      for(std::vector<Vec3>::iterator iter = geometryFieldValues.begin(); iter != geometryFieldValues.end(); iter++)
      {
        for (int i = 0; i < 3; i++)
        {
          if ((*iter)[i] < localBoundingBox.min[i])
          {
            localBoundingBox.min[i] = (*iter)[i];
          }
          if ((*iter)[i] > localBoundingBox.max[i])
          {
            localBoundingBox.max[i] = (*iter)[i];
          }
        }
      }

      // reduce the bounding box values
      BoundingBox globalBoundingBox;
      MPI_Reduce(localBoundingBox.min.data(), globalBoundingBox.min.data(), 3, MPI_DOUBLE, MPI_MIN, 0, rankSubset->mpiCommunicator());
      MPI_Reduce(localBoundingBox.max.data(), globalBoundingBox.max.data(), 3, MPI_DOUBLE, MPI_MAX, 0, rankSubset->mpiCommunicator());

      LOG(DEBUG) << "reduced bounding box: " << globalBoundingBox.min << ", " << globalBoundingBox.max;

      boundingBoxValues_.resize(6);
      boundingBoxValues_[0] = globalBoundingBox.min[0];
      boundingBoxValues_[1] = globalBoundingBox.min[1];
      boundingBoxValues_[2] = globalBoundingBox.max[2];
      boundingBoxValues_[3] = globalBoundingBox.max[0];
      boundingBoxValues_[4] = globalBoundingBox.max[1];
      boundingBoxValues_[5] = globalBoundingBox.min[2];

      // reduce global number of nodes
      int nNodesLocal = geometryFieldValues.size();
      nNodesGlobal_ = 0;
      MPI_Reduce(&nNodesLocal, &nNodesGlobal_, 1, MPI_INT, MPI_SUM, 0, rankSubset->mpiCommunicator());

      // write everything
      // write geometry field

      // convert data to be send to ADIOS from vector<Vec3> to vector<double>
      geometryFieldScalarValues.resize(3*geometryFieldValues.size());
      for (int i = 0; i < geometryFieldValues.size(); i++)
      {
        for (int j = 0; j != 3; j++)
        {
          geometryFieldScalarValues[3*i + j] = geometryFieldValues[i][j];
        }
      }

      // define variable
      if (!adiosFieldVariableGeometry_)
      {
        std::string variableName = "xyz";

        // communicate offset into the global values array
        int localSize = geometryFieldScalarValues.size();
        int offset = 0;
        int globalSize = 0;

        std::shared_ptr<Partition::RankSubset> rankSubset = DihuContext::partitionManager()->rankSubsetForCollectiveOperations();

        MPI_Exscan(&localSize, &offset, 1, MPI_INT, MPI_SUM, rankSubset->mpiCommunicator());
        MPI_Allreduce(&localSize, &globalSize, 1, MPI_INT, MPI_SUM, rankSubset->mpiCommunicator());

        LOG(DEBUG) << "define variable \"" << variableName << "\", localSize: " << localSize << ", offset: " << offset << ", globalSize: " << globalSize;

        // name, global size, offset, local size
        adiosFieldVariableGeometry_ = std::make_shared<adios2::Variable<double>>(adiosIo->DefineVariable<double>(
          variableName, {(long unsigned int)globalSize}, {(long unsigned int)offset}, {(long unsigned int)localSize}, adios2::ConstantDims
        ));
      }

      // write data to file
      writer->engine->Put<double>(*adiosFieldVariableGeometry_.get(), geometryFieldScalarValues.data());


      // write scalar field
      // define variable
      if (!adiosFieldVariableScalar_)
      {
        std::string variableName = "i";

        // communicate offset and global size
        int localSize = scalarFieldVariableValues.size();
        int offset = 0;
        int globalSize = 0;

        MPI_Exscan(&localSize, &offset, 1, MPI_INT, MPI_SUM, rankSubset->mpiCommunicator());
        MPI_Allreduce(&localSize, &globalSize, 1, MPI_INT, MPI_SUM, rankSubset->mpiCommunicator());

        LOG(DEBUG) << "define variable \"" << variableName << "\", localSize: " << localSize << ", offset: " << offset << ", globalSize: " << globalSize;

        // name, global size, offset, local size
        adiosFieldVariableScalar_ = std::make_shared<adios2::Variable<double>>(adiosIo->DefineVariable<double>(
          variableName, {(long unsigned int)globalSize}, {(long unsigned int)offset}, {(long unsigned int)localSize}, adios2::ConstantDims
        ));
      }

      // write data to file
      writer->engine->Put<double>(*adiosFieldVariableScalar_.get(), scalarFieldVariableValues.data());

      // write variables that have to be written by only one rank (rank 0)
      if (rankSubset->ownRankNo() == 0)
      {
        // write box
        if (!boxVariable_)
        {
          boxVariable_ = std::make_shared<adios2::Variable<double>>(adiosIo->DefineVariable<double>("box", {6},{0},{6}, adios2::ConstantDims));
        }
        writer->engine->Put<double>(*boxVariable_.get(), boundingBoxValues_.data());

        // write radius
        globalRadius_ = 0.1;
        if (approximateDistanceBetweenFibers > 1e-3)
        {
          globalRadius_ = approximateDistanceBetweenFibers*0.1;
        }
        if (!globalRadiusVariable_)
        {
          globalRadiusVariable_ = std::make_shared<adios2::Variable<double>>(adiosIo->DefineVariable<double>("global_radius"));
        }
        writer->engine->Put<double>(*globalRadiusVariable_.get(), &globalRadius_);

        // write global number of nodes
        if (!globalNumberOfNodesVariable_)
        {
          globalNumberOfNodesVariable_ = std::make_shared<adios2::Variable<int>>(adiosIo->DefineVariable<int>("p_count"));
        }
        writer->engine->Put<int>(*globalNumberOfNodesVariable_.get(), &nNodesGlobal_);
      }

    }   // endif writer->nOpenWriters == combineNInstances_

    //LOG(DEBUG) << "timeStepCloseInterval: " << timeStepCloseInterval << ", writeCallCount_: " << writeCallCount_;

    //int timeStepCloseInterval = specificSettings_.getOptionInt("timeStepCloseInterval",1);
    //if (writeCallCount_ % timeStepCloseInterval == 0 && writeCallCount_ > 0)

    LOG(DEBUG) << "in MegaMol::write, closing time step, currentOpenWriterIndex_: " << currentOpenWriterIndex_
      << ", nOpenWriters: " << writer->nOpenWriters << ", combineNInstances_: " << combineNInstances_;

    if (writer->nOpenWriters > combineNInstances_)
    {
      LOG(FATAL) << "More open MegaMol writers than allowed.";
    }

    if (writer->nOpenWriters == combineNInstances_)
    {
      LOG(DEBUG) << "shutdown writer";

      // end output for current time step
      writer->engine->EndStep();

      writer->engine->Close();
      writer->engine = nullptr;

      // clear context
      writer->megaMolWriterContext.geometryFieldValues.clear();
      writer->megaMolWriterContext.scalarFieldVariableValues.clear();
      writer->megaMolWriterContext.approximateDistanceBetweenFibers = -1;

      lastFilename_ = currentFilename_;

      if (currentOpenWriterIndex_ == 1)
      {
        currentOpenWriterIndex_ = 0;
      }
      else
      {
        currentOpenWriterIndex_ = 1;
      }
      writer->nOpenWriters = 0;
      LOG(DEBUG) << "closing time step and writer, next currentOpenWriterIndex will be " << currentOpenWriterIndex_;

#ifdef HAVE_MEGAMOL
      notifyMegaMol();
#endif
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
