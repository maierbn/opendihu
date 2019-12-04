#pragma once

#include <Python.h>  // has to be the first included header
#include <petscmat.h>
#include <memory>

#include "data_management/data.h"
#include "field_variable/field_variable.h"

namespace Data
{

/**  The datastructures used for the dynamic solver of solid mechanics problems.
 *   It only contains field variables used for the timestepping (u,v,a),
 *   the static field variables are stored in a different Data class.
  */
template<typename FunctionSpaceType>
class DynamicHyperelasticitySolver :
  public Data<FunctionSpaceType>
{
public:

  typedef FieldVariable::FieldVariable<FunctionSpaceType,3> DisplacementsFieldVariableType;

  //! constructor
  DynamicHyperelasticitySolver(DihuContext context);

  //! field variable of u
  std::shared_ptr<DisplacementsFieldVariableType> displacements();

  //! field variable of v
  std::shared_ptr<DisplacementsFieldVariableType> velocity();

  //! field variable of a
  std::shared_ptr<DisplacementsFieldVariableType> acceleration();

  //! field variable of u_compressible
  std::shared_ptr<DisplacementsFieldVariableType> displacementsCompressible();

  //! field variable of v_compressible
  std::shared_ptr<DisplacementsFieldVariableType> velocityCompressible();

  //! ∂Wint_compressible, before incompressibility is enforced
  std::shared_ptr<DisplacementsFieldVariableType> internalVirtualWorkCompressible();

  //! ∂Wext
  std::shared_ptr<DisplacementsFieldVariableType> externalVirtualWork();

  //! initialize, get the displacements from the static solver
  void initialize(std::shared_ptr<DisplacementsFieldVariableType> displacements);

  //! field variables that will be output by outputWriters
  //type to use if there is no fiber direction field variable
  typedef std::tuple<
    std::shared_ptr<DisplacementsFieldVariableType>,  // current geometry field
    std::shared_ptr<DisplacementsFieldVariableType>,  // displacements
    std::shared_ptr<DisplacementsFieldVariableType>,  // velocity
    std::shared_ptr<DisplacementsFieldVariableType>,   // acceleration
    std::shared_ptr<DisplacementsFieldVariableType>,   // u_compressible
    std::shared_ptr<DisplacementsFieldVariableType>,   // v_compressible
    std::shared_ptr<DisplacementsFieldVariableType>,   // ∂Wint_compressible
    std::shared_ptr<DisplacementsFieldVariableType>   // ∂Wext
  >
  FieldVariablesForOutputWriter;

  //! get pointers to all field variables that can be written by output writers
  FieldVariablesForOutputWriter getFieldVariablesForOutputWriter();

protected:

  //! initializes the vectors with size
  void createPetscObjects() override;

  std::shared_ptr<DisplacementsFieldVariableType> displacements_;       //< the displacements, u
  std::shared_ptr<DisplacementsFieldVariableType> velocity_;     //< v, the velocity
  std::shared_ptr<DisplacementsFieldVariableType> acceleration_;     //< u, the acceleration
  std::shared_ptr<DisplacementsFieldVariableType> internalVirtualWorkCompressible_;     //< ∂Wint_compressible, before incompressibility is enforced
  std::shared_ptr<DisplacementsFieldVariableType> displacementsCompressible_;       //< u_compressible
  std::shared_ptr<DisplacementsFieldVariableType> velocityCompressible_;     //< v_compressible
  std::shared_ptr<DisplacementsFieldVariableType> externalVirtualWork_;     //< ∂Wext
};

} // namespace Data

#include "data_management/specialized_solver/dynamic_hyperelasticity_solver.tpp"
