#include "data_management/data.h"

#include <iomanip>
#include <iostream>
#include <sstream>
#include <vector>
#include <numeric>

#include "easylogging++.h"

#include "control/python_utility.h"
#include "control/dihu_context.h"

namespace Data
{

Data::Data(DihuContext &context) : context_(context)
{
}

Data::~Data()
{
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
