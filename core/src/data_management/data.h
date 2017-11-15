#pragma once

#include <petscmat.h>
#include <memory>

#include "control/types.h"
#include "mesh/mesh.h"

class DihuContext;

namespace Data
{

class Data
{
public:
 
  //! constructor
  Data(DihuContext &context);
  
  //! destructor
  virtual ~Data();
 
  //! initialize the mesh with e.g. number of dimensions
  void setMesh(std::shared_ptr<Mesh::Mesh> mesh);
  
  //! initialize the number of degrees of freedom per mesh node
  void setNDegreesOfFreedomPerNode(int n);
  
  //! get the stored mesh
  std::shared_ptr<Mesh::Mesh> mesh();
  
  //! return the total number of degrees of freedom, this can be a multiple of the number of nodes of the mesh
  int nDegreesOfFreedom();
  
  //! return the number of degrees of freedom per mesh node
  int nDegreesOfFreedomPerNode();
  
  //! return the solution vector that will be written by an output writer
  virtual Vec &solution() = 0;
  
protected:
 
  //! initializes the vectors and stiffness matrix with size
  virtual void createPetscObjects() = 0;
 
  DihuContext &context_;
  
  std::shared_ptr<Mesh::Mesh> mesh_;
  int nDegreesOfFreedomPerNode_;    ///< number of degrees of freedom per mesh node
  
  bool initialized_ = false;
};

}