#pragma once

#include <Python.h>  // has to be the first included header
#include <petscmat.h>
#include <memory>

#include "control/types.h"
#include "mesh/mesh.h"
#include "data_management/data.h"
#include "field_variable/field_variable.h"

namespace Data
{

/**  The datastructures for the chaste adapter where chaste computes a static hyperelastic problem.
 *   The field variables here will receive the resulting values from the chaste data structures.
  */
template<typename FunctionSpace>
class QuasiStaticNonlinearElasticityChaste :
  public Data<FunctionSpace>
{
public:

  typedef FieldVariable::FieldVariable<FunctionSpace,3> VectorFieldVariableType;
  typedef FieldVariable::FieldVariable<FunctionSpace,1> FieldVariableType;
  typedef FieldVariable::FieldVariable<FunctionSpace,9> StressFieldVariableType;

  //! constructor
  QuasiStaticNonlinearElasticityChaste(DihuContext context);

  //! return the field variable of the activation factor
  std::shared_ptr<FieldVariableType> activation();

  //! return the field variable of the active stress
  std::shared_ptr<StressFieldVariableType> activeStress();

  //! return the field variable of the displacement
  std::shared_ptr<VectorFieldVariableType> displacement();

  //! initialize
  void initialize();

  //! print all stored data to stdout
  void print();

  //! field variables that will be output by outputWriters
  typedef std::tuple<
      std::shared_ptr<FieldVariableType>,              // activation
      std::shared_ptr<StressFieldVariableType>,         // active stress
      std::shared_ptr<VectorFieldVariableType>         // displacement
    >
   FieldVariablesForOutputWriter;

  //! get pointers to all field variables that can be written by output writers
  FieldVariablesForOutputWriter getFieldVariablesForOutputWriter();

private:

  //! initializes the vectors with size
  void createPetscObjects() override;

  std::shared_ptr<FieldVariableType> activation_; ///< field variable of the activation factor field
  std::shared_ptr<StressFieldVariableType> activeStress_; ///< field variable of the active stress in the muscle
  std::shared_ptr<VectorFieldVariableType> displacement_; ///< field variable of the displacement

};

} // namespace Data

#include "data_management/specialized_solver/quasi_static_nonlinear_elasticity_chaste.tpp"
