#include "mesh/mesh.h"

namespace Mesh
{

Mesh::Mesh(PythonConfig specificSettings) :
  specificSettings_(specificSettings), initialized_(false)
{

}

//! set the name of the mesh
void Mesh::setMeshName(std::string meshName)
{
  meshName_ = meshName;
}

//! get the name of the mesh
std::string Mesh::meshName()
{
  return meshName_;
}

void Mesh::initialized()
{
  return initialized_;
}

}
