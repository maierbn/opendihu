#include "data_management/data.h"

#include <iomanip>
#include <iostream>
#include <sstream>
#include <vector>
#include <numeric>

#include "easylogging++.h"

#include "utility/python_utility.h"
#include "control/dihu_context.h"

namespace Data
{

Data::Data(const DihuContext &context) : 
  context_(context), nDegreesOfFreedomPerNode_(1)
{
}

Data::~Data()
{
}

void Data::setNDegreesOfFreedomPerNode(int n)
{
  nDegreesOfFreedomPerNode_ = n;
}

int Data::nDegreesOfFreedomPerNode()
{
  return nDegreesOfFreedomPerNode_;
}

int Data::nDegreesOfFreedom()
{
  return mesh_->nNodes() * nDegreesOfFreedomPerNode_;
}


void Data::setMesh(std::shared_ptr<Mesh::Mesh> mesh)
{
  mesh_ = mesh;
  
  if (!initialized_)
  {
    createPetscObjects();
    initialized_ = true;
  }
  else
  {
    LOG(WARNING) << "Mesh is already assigned";
  }
}

std::shared_ptr<Mesh::Mesh> Data::mesh()
{
  return mesh_;
}

} // namespace Data
