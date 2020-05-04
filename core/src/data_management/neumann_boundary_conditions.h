#pragma once

#include <Python.h>  // has to be the first included header

#include "data_management/data.h"
#include "control/types.h"
#include "mesh/mesh.h"
#include "partition/partitioned_petsc_mat/partitioned_petsc_mat.h"

namespace Data
{

template<typename FunctionSpaceType,int nComponents>
class NeumannBoundaryConditions :
  public Data<FunctionSpaceType>
{
public:

  //! use constructor of base class
  using Data<FunctionSpaceType>::Data;

  //! destructor
  virtual ~NeumannBoundaryConditions();

  //! initialize the object, create all stored data
  virtual void initialize() override;

  //! reset the object and deallocate matrices
  virtual void reset() override;

  //! return reference to a stiffness matrix
  std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,nComponents>> rhs();

private:

  //! initializes the vectors and stiffness matrix with size
  void createPetscObjects();

  std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,nComponents>> rhs_;  //< the rhs vector contribution from Neumann BC

};
}  // namespace

#include "data_management/neumann_boundary_conditions.tpp"
