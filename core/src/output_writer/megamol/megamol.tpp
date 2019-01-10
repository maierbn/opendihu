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
void MegaMOL::write(DataType& data, int timeStepNo, double currentTime)
{
  // check if output should be written in this timestep and prepare filename
  if (!Generic::prepareWrite(data, timeStepNo, currentTime))
  {
    return;
  }

  // collect all available meshes
  std::set<std::string> meshNames;
  LoopOverTuple::loopCollectMeshNames<typename DataType::OutputFieldVariables>(data.getOutputFieldVariables(), meshNames);

  // loop over meshes
  for (std::string meshName : meshNames)
  {
    // loop over all field variables and output those that are associated with the mesh given by meshName
    MegaMOLLoopOverTuple::loopOutput(data.getOutputFieldVariables(), data.getOutputFieldVariables(), meshName, specificSettings_);
  }
}

};  // namespace
