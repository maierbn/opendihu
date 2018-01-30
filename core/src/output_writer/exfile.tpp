#include "output_writer/paraview.h"

#include <fstream>
#include <iostream>
#include <thread>
#include <chrono>

#include "easylogging++.h"
#include "base64.h"

#include "utility/python_utility.h"
#include "utility/petsc_utility.h"
#include "mesh/regular_fixed.h"
#include "mesh/structured_deformable.h"
#include "mesh/unstructured_deformable.h"
#include "mesh/mesh.h"
#include "basis_on_mesh/05_basis_on_mesh.h"
#include "output_writer/exfile_writer.h"

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
  
  LOG(DEBUG) << "output exfile";
  
  typedef typename DataType::BasisOnMesh MeshType;
  std::shared_ptr<MeshType> mesh = std::static_pointer_cast<MeshType>(data.mesh());
  
  // exelem file
  // determine file name
  std::stringstream s;
  s<<filename_<<".exelem";
  std::string filenameExelem = s.str();

  // open file
  std::ofstream file = openFile(filenameExelem); 
  ExfileWriter<typename DataType::BasisOnMesh>::outputExelem(file, data.fieldVariables()); 
  file.close();
  
  // exnode file
  s.str("");
  s<<filename_<<".exnode";
  std::string filenameExnode = s.str();

  // open file
  file = openFile(filenameExnode);
  ExfileWriter<typename DataType::BasisOnMesh>::outputExnode(file, data.fieldVariables()); 
  file.close();
}

};
