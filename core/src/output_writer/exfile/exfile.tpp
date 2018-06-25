#include "output_writer/exfile/exfile.h"

#include <fstream>
#include <iostream>
#include <thread>
#include <chrono>

#include "easylogging++.h"
#include "base64.h"

#include "utility/python_utility.h"
#include "utility/petsc_utility.h"
#include "mesh/structured_regular_fixed.h"
#include "mesh/structured_deformable.h"
#include "mesh/unstructured_deformable.h"
#include "mesh/mesh.h"
#include "basis_on_mesh/basis_on_mesh.h"
#include "output_writer/exfile/exfile_writer.h"
#include "output_writer/exfile/loop_output_exelem.h"
#include "output_writer/exfile/loop_output_exnode.h"

namespace OutputWriter
{

template<typename DataType>
void Exfile::write(DataType& data, int timeStepNo, double currentTime)
{
  LOG(TRACE) << "Exfile::write";

  // check if output should be written in this timestep and prepare filename
  if (!Generic::prepareWrite(data, timeStepNo, currentTime))
  {
    return;
  }

  // collect all available meshes
  std::set<std::string> meshNames;
  LoopOverTuple::loopCollectMeshNames<typename DataType::OutputFieldVariables>(data.getOutputFieldVariables(), meshNames);
  
  LOG(DEBUG) << "collected meshNames: ";
  for (std::string meshName : meshNames)
  {
    LOG(DEBUG) << "  " << meshName;
  }
  
  // loop over meshes and create a exelem/exnode pair for each
  for (std::string meshName : meshNames)
  {
    // setup name of files
    std::stringstream filenameStart;
    if (meshNames.size() == 1)
      filenameStart << filename_;
    else
      filenameStart << filename_ << "_" << meshName;
   
    // exelem file
    // determine file name
    std::stringstream s;
    s << filenameStart.str() << ".exelem";
    std::string filenameExelem = s.str();
    

    // open file
    std::ofstream file = openFile(filenameExelem);
    // output the exelem file for all field variables that are defined on the specified meshName
    std::shared_ptr<Mesh::Mesh> mesh = nullptr;
    ExfileLoopOverTuple::loopOutputExelem(data.getOutputFieldVariables(), data.getOutputFieldVariables(), meshName, file, mesh);
    file.close();

    // exnode file
    s.str("");
    s << filenameStart.str() << ".exnode";
    std::string filenameExnode = s.str();

    // open file
    file = openFile(filenameExnode);
    // output the exnode file for all field variables that are defined on the specified meshName
    ExfileLoopOverTuple::loopOutputExnode(data.getOutputFieldVariables(), data.getOutputFieldVariables(), meshName, file);
    file.close();

    // store created filename
    FilenameWithElementAndNodeCount item;
    item.filename = filenameStart.str();
    item.nElements = mesh->nElements();
    item.nNodes = mesh->nNodes();
    item.meshName = meshName;
    item.dimensionality = mesh->dimension();
    
    //if (filenamesWithElementAndNodeCount_.find(currentTime) != filenamesWithElementAndNodeCount_.end())
    //{
    filenamesWithElementAndNodeCount_[currentTime].push_back(item);
    //}
  }
  
  
  // output visualization file
  outputComFile();
}

};
