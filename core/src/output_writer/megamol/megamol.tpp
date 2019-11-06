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
void MegaMol::write(DataType& data, int timeStepNo, double currentTime, int callCountIncrement)
{
  // check if output should be written in this timestep and prepare filename
  if (!Generic::prepareWrite(data, timeStepNo, currentTime, callCountIncrement))
  {
    return;
  }
#ifdef HAVE_ADIOS

  std::shared_ptr<Partition::RankSubset> rankSubset = context_.partitionManager()->rankSubsetForCollectiveOperations();

  try
  {
    // create ADIOS writer, if it does not yet exist
    if (useFrontBackBuffer_)
    {
      // entry in map was not yet created, initialize with nullptr
      if (adiosWriters_.find(this->filenameBase_) == adiosWriters_.end())
      {
        adiosWriters_[this->filenameBase_] = std::array<std::shared_ptr<MegaMol::adios_writer_t>, 2>({nullptr,nullptr});
      }

      if (!adiosWriters_[this->filenameBase_][currentOpenWriterIndex_])
      {
        adiosWriters_[this->filenameBase_][currentOpenWriterIndex_] = std::make_shared<adios_writer_t>();
        adiosWriters_[this->filenameBase_][currentOpenWriterIndex_]->nOpenWriters = 0;
        adiosWriters_[this->filenameBase_][currentOpenWriterIndex_]->adiosIo = std::make_shared<adios2::IO>(
          context_.adios()->DeclareIO(this->filenameBase_+std::string(1,'0'+currentOpenWriterIndex_)));
        assert(adiosWriters_[this->filenameBase_][currentOpenWriterIndex_]->adiosIo);
      }
      adiosWriter_[this->filenameBase_] = adiosWriters_[this->filenameBase_][currentOpenWriterIndex_];
    }
    else
    {
      if (!(adiosWriter_[this->filenameBase_]))
      {
        adiosWriter_[this->filenameBase_] = std::make_shared<adios_writer_t>();
        adiosWriter_[this->filenameBase_]->nOpenWriters = 0;
        adiosWriter_[this->filenameBase_]->adiosIo = std::make_shared<adios2::IO>(
          context_.adios()->DeclareIO(this->filenameBase_));
        assert(adiosWriter_[this->filenameBase_]->adiosIo);
      }
    }
    std::shared_ptr<adios_writer_t> &currentWriter = adiosWriter_[this->filenameBase_];

    // if writer is not open, open it
    if (currentWriter->nOpenWriters == 0)
    {
      assert(currentWriter->adiosIo);

      // determine current filename
      std::stringstream outputFileName;
      if (useFrontBackBuffer_)
      {
        outputFileName << this->filenameBase_ << "_" << currentOpenWriterIndex_;
      }
      else
      {
        //outputFileName << filenameBaseWithNo_;
        outputFileName << filenameBase_;
      }

      currentFilename_ = outputFileName.str();

      if (!currentWriter->engine)
      {
        LOG(DEBUG) << "create new engine, on outputFileName: \"" << outputFileName.str() << "\"";

        // create new writer
        adios2::Engine engine = currentWriter->adiosIo->Open(outputFileName.str(), adios2::Mode::Write, rankSubset->mpiCommunicator());
        currentWriter->engine = std::make_shared<adios2::Engine>(engine);

        currentWriter->megaMolWriterContext.approximateDistanceBetweenFibers = -1;  // initialize to not set

        // write python script and meta data as attribute
        currentWriter->adiosIo->DefineAttribute<std::string>("config", this->context_.pythonScriptText());
        currentWriter->adiosIo->DefineAttribute<std::string>("version", this->context_.versionText());
        currentWriter->adiosIo->DefineAttribute<std::string>("meta", this->context_.metaText());
      }

      // begin output to file for current time step
      currentWriter->engine->BeginStep();
    }
    currentWriter->nOpenWriters++;
    LOG(DEBUG) << "in MegaMol::write \"" << this->filenameBase_ << "\", "
      << "opening time step, currentOpenWriterIndex_: " << currentOpenWriterIndex_ << ", nOpenWriters: " << currentWriter->nOpenWriters;

    // collect all available meshes
    std::set<std::string> meshNames;
    LoopOverTuple::loopCollectMeshNames<typename DataType::FieldVariablesForOutputWriter>(data.getFieldVariablesForOutputWriter(), meshNames);

    // loop over meshes and retrieve data
    for (std::string meshName : meshNames)
    {
      // loop over all field variables and output those that are associated with the mesh given by meshName
      MegaMolLoopOverTuple::loopOutput(data.getFieldVariablesForOutputWriter(), data.getFieldVariablesForOutputWriter(), meshName, specificSettings_,
                                       currentWriter->megaMolWriterContext);
    }

    // only write the variables to file when all collecting calls have been made
    if (currentWriter->nOpenWriters == combineNInstances_)
    {
      // write everything
      writeAdiosVariables();
    }   // endif currentWriter->nOpenWriters == combineNInstances_

    //LOG(DEBUG) << "timeStepCloseInterval: " << timeStepCloseInterval << ", writeCallCount_: " << writeCallCount_;

    //int timeStepCloseInterval = specificSettings_.getOptionInt("timeStepCloseInterval",1);
    //if (writeCallCount_ % timeStepCloseInterval == 0 && writeCallCount_ > 0)

    LOG(DEBUG) << "in MegaMol::write, \"" << this->filenameBase_ << " closing time step, currentOpenWriterIndex_: " << currentOpenWriterIndex_
      << ", nOpenWriters: " << currentWriter->nOpenWriters << ", combineNInstances_: " << combineNInstances_;

    if (currentWriter->nOpenWriters > combineNInstances_)
    {
      LOG(FATAL) << "More open MegaMol writers (" << currentWriter->nOpenWriters << ") on file \"" << this->filenameBase_
        << "\" than allowed (combineNInstances_ = " << combineNInstances_ << ")";
    }

    if (currentWriter->nOpenWriters == combineNInstances_)
    {
      LOG(DEBUG) << "shutdown writer on \"" << this->filenameBase_ << "\"";

      // end output for current time step
      currentWriter->engine->EndStep();

      /*currentWriter->engine->Close();
      currentWriter->engine = nullptr;
      */

      // clear context
      currentWriter->megaMolWriterContext.geometryFieldValues.clear();
      currentWriter->megaMolWriterContext.vmValues.clear();
      currentWriter->megaMolWriterContext.emgValues.clear();
      currentWriter->megaMolWriterContext.transmembraneFlowValues.clear();
      currentWriter->megaMolWriterContext.approximateDistanceBetweenFibers = -1;

      lastFilename_ = currentFilename_;

      if (currentOpenWriterIndex_ == 1)
      {
        currentOpenWriterIndex_ = 0;
      }
      else
      {
        currentOpenWriterIndex_ = 1;
      }
      currentWriter->nOpenWriters = 0;
      LOG(DEBUG) << "closing time step and writer \"" << this->filenameBase_ << "\""
        << ", next currentOpenWriterIndex will be " << currentOpenWriterIndex_;

#ifdef HAVE_MEGAMOL
      notifyMegaMol();
#endif
    }
  }
  catch (std::invalid_argument &e)
  {
    LOG(ERROR) << "Writer \"" << this->filenameBase_ << "\": Invalid argument exception";
    LOG(ERROR) << e.what();
  }
  catch (std::ios_base::failure &e)
  {
    LOG(ERROR) << "Writer \"" << this->filenameBase_ << "\": System exception";
    LOG(ERROR) << e.what();
  }
  catch (std::exception &e)
  {
    LOG(ERROR) << "Writer \"" << this->filenameBase_ << "\": Exception";
    LOG(ERROR) << e.what();
  }

#endif
}

} // namespace
