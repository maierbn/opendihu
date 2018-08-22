#include "mesh/mesh.h"

namespace Mesh
{

Mesh::Mesh(PyObject *specificSettings) : 
  specificSettings_(specificSettings)
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

}