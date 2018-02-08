#pragma once

#include <Python.h>  // has to be the first included header
#include <petscmat.h>
#include <memory>
#include <vector>

#include "control/types.h"
#include "mesh/mesh.h"
#include "field_variable/field_variable.h"

class DihuContext;

namespace Data
{

template<typename BasisOnMeshType>
class Data
{
public:
 
  typedef BasisOnMeshType BasisOnMesh;
 
  //! constructor
  Data(const DihuContext &context);
  
  //! destructor
  virtual ~Data();
 
  //! initialize the mesh with e.g. number of dimensions
  virtual void setMesh(std::shared_ptr<BasisOnMeshType> mesh);
  
  //! initialize the number of degrees of freedom per mesh node, i.e. the number of components of the field variables
  void setNComponentsPerNode(int n);
  
  //! get the stored mesh
  std::shared_ptr<BasisOnMeshType> mesh();
  
  //! return the total number of degrees of freedom, this can be a multiple of the number of nodes of the mesh
  int nDegreesOfFreedom();
  
  //! return the number of degrees of freedom per mesh node
  int nComponentsPerNode();
  
  //! get all stored field variables, this is used by the output writers
  virtual std::vector<std::shared_ptr<FieldVariable::FieldVariable<BasisOnMeshType>>> fieldVariables() = 0;
  
  //! return the solution vector that can be written by an output writer
  virtual FieldVariable::FieldVariable<BasisOnMeshType> &solution() = 0;
  
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
