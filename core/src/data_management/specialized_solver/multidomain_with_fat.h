#pragma once

#include <Python.h>  // has to be the first included header
#include <petscmat.h>
#include <memory>

#include "data_management/specialized_solver/multidomain.h"

namespace Data
{

/**  The datastructures used for multidomain solver.
  */
template<typename FunctionSpaceType, typename FunctionSpaceFatType>
class MultidomainWithFat : public Data<FunctionSpaceType>
{
public:

  typedef FieldVariable::FieldVariable<FunctionSpaceFatType,1> FieldVariableFatType;
  typedef FieldVariable::FieldVariable<FunctionSpaceFatType,3> GeometryFieldVariableFatType;
  typedef typename Multidomain<FunctionSpaceType>::FieldVariableType FieldVariableType;
  typedef typename Multidomain<FunctionSpaceType>::GradientFieldVariableType GradientFieldVariableType;

  //! constructor
  MultidomainWithFat(DihuContext context);

  //! assign the data object of MultidomainSolver which contains all data without the fat layer
  void setDataMultidomain(std::shared_ptr<Multidomain<FunctionSpaceType>> dataMultidomain);

  //! return a reference to the rhs summand vector which is needed to apply the boundary conditions, the PETSc Vec can be obtained via fieldVariable->valuesGlobal()
  std::shared_ptr<FieldVariableFatType> extraCellularPotentialFat();

  //! field variables that will be output by outputWriters
  typedef std::tuple<
    std::shared_ptr<GradientFieldVariableType>,     // geometry
    std::shared_ptr<GradientFieldVariableType>,     // fiberDirection
    std::shared_ptr<FieldVariableType>,              // solution of laplace potential flow
    std::shared_ptr<FieldVariableType>,              // extra-cellular potential
    std::vector<std::shared_ptr<FieldVariableType>>,              // transmembranePotentials
    std::vector<std::shared_ptr<FieldVariableType>>,              // compartmentRelativeFactors
    std::shared_ptr<FieldVariableType>,                   // relativeFactorTotal
    std::shared_ptr<GeometryFieldVariableFatType>,   // fat layer geometry
    std::shared_ptr<FieldVariableFatType>            // fat layer phi_e
  > FieldVariablesForOutputWriter;

  //! get pointers to all field variables that can be written by output writers
  FieldVariablesForOutputWriter getFieldVariablesForOutputWriter();

private:

  //! initializes the vectors with size
  void createPetscObjects() override;

  std::shared_ptr<Multidomain<FunctionSpaceType>> dataMultidomain_;  //< data object of multidomain
  std::shared_ptr<FieldVariableFatType> extraCellularPotentialFat_; //< phi_e value on fat layer
};

} // namespace Data

#include "data_management/specialized_solver/multidomain_with_fat.tpp"
