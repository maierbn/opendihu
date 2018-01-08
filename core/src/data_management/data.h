#pragma once

#include <petscmat.h>
#include <memory>

#include "control/types.h"
#include "mesh/mesh.h"

class DihuContext;

namespace Data
{

 // TODO: template, nDegreesOfFreedom->n components, Vec -> FieldVariable
template<typename BasisOnMeshType>
class Data
{
public:
 
  //! constructor
  Data(const DihuContext &context);
  
  //! destructor
  virtual ~Data();
 
  //! initialize the mesh with e.g. number of dimensions
  void setMesh(std::shared_ptr<BasisOnMeshType> mesh);
  
  //! initialize the number of degrees of freedom per mesh node, i.e. the number of components of the field variables
  void setNComponentsPerNode(int n);
  
  //! get the stored mesh
  std::shared_ptr<BasisOnMeshType> mesh();
  
  //! return the total number of degrees of freedom, this can be a multiple of the number of nodes of the mesh
  int nDegreesOfFreedom();
  
  //! return the number of degrees of freedom per mesh node
  int nComponentsPerNode();
  
  //! return the solution vector that can be written by an output writer
  virtual Vec &solution() = 0;
  
protected:
 
  //! initializes the vectors and stiffness matrix with size
  virtual void createPetscObjects() = 0;
 
  const DihuContext &context_;
  
  std::shared_ptr<BasisOnMeshType> mesh_;
  int nComponentsPerNode_;    ///< number of degrees of freedom per mesh node
  
  bool initialized_ = false;
};

} // namespace

#include "data_management/data.tpp"