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

namespace OutputWriter
{

template<typename DataType>
void Exfile::write(DataType& data, int timeStepNo, double currentTime)
{
  // check if output should be written in this timestep and prepare filename
  if (!Generic::prepareWrite(data, timeStepNo, currentTime))
  {
    return;
  }
  
  typedef typename DataType::BasisOnMesh MeshType;
  std::shared_ptr<MeshType> mesh = std::static_pointer_cast<MeshType>(data.mesh());
  
  // exelem file
  // determine file name
  std::stringstream s;
  s<<filename_<<".exelem";
  std::string filenameExelem= s.str();

  // open file
  std::ofstream file = openFile(filenameExelem);  
  file.close();
  
  mesh->outputExelemFile(file);
  
  // exnode file
  s.str("");
  s<<filename_<<".exelem";
  std::string filenameExnode= s.str();

  // open file
  file = openFile(filenameExnode);
  
  mesh->outputExnodeFile(file);
  file.close();
}

};