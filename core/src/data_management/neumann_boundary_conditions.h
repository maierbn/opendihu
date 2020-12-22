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

  //! return reference to the rhs field variable
  std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,nComponents>> rhs();

  //! return reference to the deformation gradient field variable
  std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,9>> &deformationGradient();

private:

  //! initializes the vectors and stiffness matrix with size
  void createPetscObjects();

  std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,nComponents>> rhs_;  //< the rhs vector contribution from Neumann BC
  std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,9>> deformationGradient_;  //< the deformation gradient in the field, which is needed to transform traction bc in current configuration to traction bc in reference configuration

};
}  // namespace

#include "data_management/neumann_boundary_conditions.tpp"
