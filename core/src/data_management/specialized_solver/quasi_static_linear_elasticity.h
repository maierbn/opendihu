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

/**  The datastructures used for static bidomain solver.
  */
template<typename DataLinearElasticityType>
class QuasiStaticLinearElasticity : public Data<FunctionSpaceType>
{
public:

  typedef typename DataLinearElasticityType::FunctionSpaceType FunctionSpaceType;
  typedef FieldVariable::FieldVariable<FunctionSpaceType,1> FieldVariableType;
  typedef FieldVariable::FieldVariable<FunctionSpaceType,9> StressFieldVariableType;

  //! constructor
  QuasiStaticLinearElasticity(DihuContext context);

  //! return the field variable of the activation factor
  std::shared_ptr<FieldVariableType> activation();

  //! return the field variable of the active stress
  std::shared_ptr<StressFieldVariableType> activeStress();

  //! initialize
  void initialize();

  //! print all stored data to stdout
  void print();

  //! store the data of the finite element type
  void setData(std::shared_ptr<DataLinearElasticityType> dataLinearElasticity);

  //! field variables that will be output by outputWriters
  typedef decltype(std::tuple_cat(
    typename DataLinearElasticityType::OutputFieldVariables,
    std::tuple<
      std::shared_ptr<FieldVariableType>,              // activation
      std::shared_ptr<StressFieldVariableType>         // active stress
    >)())
  > OutputFieldVariables;

  //! get pointers to all field variables that can be written by output writers
  OutputFieldVariables getOutputFieldVariables();

private:

  //! initializes the vectors with size
  void createPetscObjects() override;

  std::shared_ptr<DataLinearElasticityType> dataLinearElasticity_;   ///< data object of the linear elasticity data class

  std::shared_ptr<FieldVariableType> activation_; ///< field variable of the activation factor field
  std::shared_ptr<StressFieldVariableType> activeStress_; ///< field variable of the active stress in the muscle
};

} // namespace Data

#include "data_management/specialized_solver/static_bidomain.tpp"
