#pragma once

#include <Python.h>  // has to be the first included header
#include <petscmat.h>
#include <memory>
#include <vector>

#include "control/types.h"
#include "mesh/mesh.h"
#include "field_variable/field_variable.h"
#include "control/dihu_context.h"

namespace Data
{

template<typename BasisOnMeshType>
class Data
{
public:
 
  typedef BasisOnMeshType BasisOnMesh;
 
  //! constructor
  Data(DihuContext context);
  
  //! destructor
  virtual ~Data();
 
  //! initialize the mesh with e.g. number of dimensions
  virtual void setMesh(std::shared_ptr<BasisOnMeshType> mesh);
  
  //! initialize, generate petsc objects, this has to be called after setMesh
  virtual void initialize();
  
  //! get the stored mesh
  const std::shared_ptr<BasisOnMeshType> mesh() const;
  
  //! return the total number of degrees of freedom, this can be a multiple of the number of nodes of the mesh
  virtual int nDegreesOfFreedom();
  
protected:
 
  //! initializes the vectors and stiffness matrix with size
  virtual void createPetscObjects() = 0;
 
  const DihuContext context_;     ///< the context object with python config of the class that uses this data object
  std::shared_ptr<BasisOnMeshType> mesh_;
  
  bool initialized_ = false;
};

} // namespace

#include "data_management/data.tpp"
