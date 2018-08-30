#pragma once

#include <Python.h>  // has to be the first included header
#include <petscmat.h>
#include <memory>

#include "data_management/data.h"
#include "data_management/finite_elements_mixed.h"
#include "control/types.h"
#include "mesh/mesh.h"
#include "field_variable/field_variable.h"
#include "function_space/mixed_function_space.h"
#include "data_management/finite_elements_solid_mechanics.h"

namespace Data
{
/** partial specialization for mixed formulation
 */
template<typename LowOrderFunctionSpaceType,typename HighOrderFunctionSpaceType,typename Term>
class FiniteElements<FunctionSpace::Mixed<LowOrderFunctionSpaceType,HighOrderFunctionSpaceType>,Term> :
  public FiniteElementsSolidMechanics<HighOrderFunctionSpaceType,Term>
{
public:

  //! inherited constructor
  using FiniteElementsSolidMechanics<HighOrderFunctionSpaceType,Term>::FiniteElementsSolidMechanics;

  typedef FieldVariable::FieldVariable<HighOrderFunctionSpaceType,HighOrderFunctionSpaceType::dim()> HighOrderFieldVariableType;
  typedef FieldVariable::FieldVariable<LowOrderFunctionSpaceType,LowOrderFunctionSpaceType::dim()> LowOrderFieldVariableType;

  //! initialize the object, create all stored data
  virtual void initialize() override;

  //! return a reference to the pressure field
  FieldVariable::FieldVariable<LowOrderFunctionSpaceType,1> &pressure();

  //! set the mixed mesh
  void setFunctionSpace(std::shared_ptr<FunctionSpace::Mixed<LowOrderFunctionSpaceType,HighOrderFunctionSpaceType>> mixedFunctionSpace);

  //! get the mixed mesh, the high order mesh can be retrieved by mesh()
  std::shared_ptr<FunctionSpace::Mixed<LowOrderFunctionSpaceType,HighOrderFunctionSpaceType>> mixedMesh();

  //! field variables that will be output by outputWriters
  typedef std::tuple<
    std::shared_ptr<FieldVariable::FieldVariable<HighOrderFunctionSpaceType,3>>,  // geometryReference
    std::shared_ptr<FieldVariable::FieldVariable<HighOrderFunctionSpaceType,3>>,  // actual geometry (stored in mesh)
    std::shared_ptr<FieldVariable::FieldVariable<HighOrderFunctionSpaceType,HighOrderFunctionSpaceType::dim()>>,  // displacements
    std::shared_ptr<FieldVariable::FieldVariable<LowOrderFunctionSpaceType,1>>,  // displacements
    std::shared_ptr<FieldVariable::FieldVariable<HighOrderFunctionSpaceType,HighOrderFunctionSpaceType::dim()>>,   // residual
    std::shared_ptr<FieldVariable::FieldVariable<HighOrderFunctionSpaceType,HighOrderFunctionSpaceType::dim()>>   // externalVirtualWork
  > OutputFieldVariables;

  //! get pointers to all field variables that can be written by output writers
  OutputFieldVariables getOutputFieldVariables();

protected:

  //! initialize the geometryReference field variable from the geometry field
  void initializeFieldVariables();

  //! get the number of rows and columns to be used for setup of tangent stiffness matrix. This is different for mixed formulation.
  const dof_no_t getTangentStiffnessMatrixNRows() override;

  std::shared_ptr<FieldVariable::FieldVariable<LowOrderFunctionSpaceType,1>> pressure_;   //< the pressure field of the mixed formulation

  std::shared_ptr<FunctionSpace::Mixed<LowOrderFunctionSpaceType,HighOrderFunctionSpaceType>> mixedMesh_;   ///< the FunctionSpace object of Type FunctionSpace::Mixed

};
}  // namespace

#include "data_management/finite_elements_mixed.tpp"
