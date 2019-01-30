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

/**  The datastructures used for multidomain solver.
  */
template<typename FunctionSpaceType>
class StaticBidomain : public Data<FunctionSpaceType>
{
public:

  typedef FieldVariable::FieldVariable<FunctionSpaceType,1> FieldVariableType;
  typedef FieldVariable::FieldVariable<FunctionSpaceType,3> GradientFieldVariableType;

  //! constructor
  StaticBidomain(DihuContext context);

  //! return a reference to the rhs summand vector which is needed to apply the boundary conditions, the PETSc Vec can be obtained via fieldVariable->valuesGlobal()
  std::shared_ptr<GradientFieldVariableType> fiberDirection();

  //! return the field variable of the potential for the Laplace potential flow problem
  std::shared_ptr<FieldVariableType> flowPotential();

  //! return the extra-cellular potential field variable
  std::shared_ptr<FieldVariableType> extraCellularPotential();

  //! return the transmembrane potential (Vm) field variable
  std::shared_ptr<FieldVariableType> transmembranePotential();

  //! return the solution vector of the transmembrane potential (Vm) field variable
  std::shared_ptr<FieldVariableType> transmembraneFlow();

  //! a field variable with constant value of zero, needed for the nested rhs vector
  std::shared_ptr<FieldVariableType> zero();

  //! get a reference to the rhs matrix
  Mat &rhsMatrix();

  //! initialize
  void initialize();

  //! print all stored data to stdout
  void print();

  //! field variables that will be output by outputWriters
  typedef std::tuple<
    std::shared_ptr<GradientFieldVariableType>,     // geometry
    std::shared_ptr<GradientFieldVariableType>,     // fiberDirection
    std::shared_ptr<FieldVariableType>,              // solution of laplace potential flow
    std::shared_ptr<FieldVariableType>,              // extra-cellular potential
    std::shared_ptr<FieldVariableType>,              // transmembranePotential
    std::shared_ptr<FieldVariableType>              // transmembraneFlow
  > OutputFieldVariables;

  //! get pointers to all field variables that can be written by output writers
  OutputFieldVariables getOutputFieldVariables();

private:

  //! initializes the vectors with size
  void createPetscObjects() override;

  Mat rhsMatrix_;   ///< the rhs premultiplication matrix
  std::shared_ptr<FieldVariableType> flowPotential_; ///< solution of the laplace flow
  std::shared_ptr<GradientFieldVariableType> fiberDirection_; ///< the direction of fibers
  std::shared_ptr<FieldVariableType> transmembraneFlow_;  ///< the Vm for the next timestep, this holds the solution in the linear solver which must be different from the rhs vector
  std::shared_ptr<FieldVariableType> transmembranePotential_;  ///< the Vm value (transmembrane potential)
  std::shared_ptr<FieldVariableType> extraCellularPotential_;  ///< the phi_e value which is the extra-cellular potential
  std::shared_ptr<FieldVariableType> zero_;  ///< a field variable with constant value of zero, needed for the nested rhs vector
};

} // namespace Data

#include "data_management/specialized_solver/static_bidomain.tpp"
