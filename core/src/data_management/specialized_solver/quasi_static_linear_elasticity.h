#pragma once

#include <Python.h>  // has to be the first included header
#include <petscmat.h>
#include <memory>

#include "control/types.h"
#include "mesh/mesh.h"
#include "data_management/data.h"
#include "data_management/time_stepping/time_stepping.h"
#include "field_variable/field_variable.h"

namespace Data
{

/**  The datastructures used for quasi static linear elasticity solver.
  */
template<typename DataLinearElasticityType>
class QuasiStaticLinearElasticity :
  public Data<typename DataLinearElasticityType::FunctionSpace>
{
public:

  typedef typename DataLinearElasticityType::FunctionSpace FunctionSpaceType;
  typedef FieldVariable::FieldVariable<FunctionSpaceType,3> VectorFieldVariableType;
  typedef FieldVariable::FieldVariable<FunctionSpaceType,1> FieldVariableType;
  typedef FieldVariable::FieldVariable<FunctionSpaceType,9> StressFieldVariableType;

  //! constructor
  QuasiStaticLinearElasticity(DihuContext context);

  //! return the field variable of the activation factor
  std::shared_ptr<FieldVariableType> activation();

  //! return the field variable of the active stress
  std::shared_ptr<StressFieldVariableType> activeStress();

  //! return the field variable of the strain
  std::shared_ptr<StressFieldVariableType> strain();

  //! return a reference to the rhs summand vector which is needed to apply the boundary conditions, the PETSc Vec can be obtained via fieldVariable->valuesGlobal()
  std::shared_ptr<VectorFieldVariableType> fiberDirection();

  //! return the field variable of the potential for the Laplace potential flow problem
  std::shared_ptr<FieldVariableType> flowPotential();

  //! return the field variable of the active stress part of rhs, f_active
  std::shared_ptr<VectorFieldVariableType> rightHandSideActive();

  //! initialize
  void initialize();

  //! print all stored data to stdout
  void print();

  //! store the data of the finite element type
  void setData(std::shared_ptr<DataLinearElasticityType> dataLinearElasticity);

  //! field variables that will be output by outputWriters
  typedef decltype(std::tuple_cat(
    std::declval<typename DataLinearElasticityType::OutputFieldVariables>(),
    std::declval<std::tuple<
      std::shared_ptr<FieldVariableType>,              // activation
      std::shared_ptr<StressFieldVariableType>,         // active stress
      std::shared_ptr<StressFieldVariableType>,         // strain
      std::shared_ptr<VectorFieldVariableType>,      // rightHandSideActive_
      std::shared_ptr<VectorFieldVariableType>,      // fiberDirection
      std::shared_ptr<FieldVariableType>               // solution of laplace potential flow
    >>()))
   OutputFieldVariables;

  //! get pointers to all field variables that can be written by output writers
  OutputFieldVariables getOutputFieldVariables();

private:

  //! initializes the vectors with size
  void createPetscObjects() override;

  std::shared_ptr<DataLinearElasticityType> dataLinearElasticity_;   ///< data object of the linear elasticity data class

  std::shared_ptr<FieldVariableType> activation_; ///< field variable of the activation factor field
  std::shared_ptr<StressFieldVariableType> activeStress_; ///< field variable of the active stress in the muscle
  std::shared_ptr<StressFieldVariableType> strain_; ///< field variable of the strain in the muscle
  std::shared_ptr<VectorFieldVariableType> rightHandSideActive_; ///< field variable of the active stress part of rhs, f_active
  std::shared_ptr<FieldVariableType> flowPotential_; ///< solution of the laplace flow
  std::shared_ptr<VectorFieldVariableType> fiberDirection_; ///< the direction of fibers

};

} // namespace Data

#include "data_management/specialized_solver/quasi_static_linear_elasticity.tpp"
