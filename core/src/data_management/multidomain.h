#pragma once

#include <Python.h>  // has to be the first included header
#include <petscmat.h>
#include <memory>

#include "control/types.h"
#include "mesh/mesh.h"
#include "data_management/data.h"
#include "data_management/time_stepping.h"
#include "field_variable/field_variable.h"

namespace Data
{

/**  The datastructures used for multidomain solver.
  */
template<typename FunctionSpaceType, int nStatesCellML>
class Multidomain : public Data<FunctionSpaceType>
{
public:

  typedef FieldVariable::FieldVariable<FunctionSpaceType,1> FieldVariableType;
  typedef FieldVariable::FieldVariable<FunctionSpaceType,3> GradientFieldVariableType;
  typedef FieldVariable::FieldVariable<FunctionSpaceType,nStatesCellML> CellMLFieldVariableType;

  //! constructor
  Multidomain(DihuContext context, int nCompartments);

  //! return a reference to the rhs summand vector which is needed to apply the boundary conditions, the PETSc Vec can be obtained via fieldVariable->valuesGlobal()
  std::shared_ptr<GradientFieldVariableType> fibreDirection();

  //! return the field variable of the potential for the Laplace potential flow problem
  std::shared_ptr<FieldVariableType> flowPotential();

  //! return the extra-cellular potential field variable
  std::shared_ptr<FieldVariableType> extraCellularPotential();

  //! return the trans-membrane potential field variable (all states) for MU compartmentNo
  std::shared_ptr<CellMLFieldVariableType> transmembranePotential(int compartmentNo);

  //! return the trans-membrane potential field variable (all states) for MU compartmentNo
  std::shared_ptr<CellMLFieldVariableType> transmembranePotentialNextTimeStep(int compartmentNo);

  //! return the increment of the trans-membrane potential field variable (all states) for MU compartmentNo
  std::shared_ptr<CellMLFieldVariableType> transmembraneIncrementNextTimeStep(int compartmentNo);

  //! return field variable for -1/Cm I_ion(Vm) for the next time step, this is a single component extracted from transmembranePotentialNextTimeStep
  std::shared_ptr<FieldVariableType> ionicCurrentNextTimeStep(int compartmentNo);

  //! return the transmembrane (Vm) increment field variable used for time stepping
  std::shared_ptr<FieldVariableType> transmembraneIncrement();


  //! print all stored data to stdout
  void print() override;

  //! field variables that will be output by outputWriters
  typedef std::tuple<
    GradientFieldVariableType,     // fibreDirection
    FieldVariableType,              // extra-cellular potential
    std::vector<FunctionSpaceType>    // transmembrane potentials
  > OutputFieldVariables;

  //! get pointers to all field variables that can be written by output writers
  OutputFieldVariables getOutputFieldVariables();

private:

  //! initializes the vectors with size
  void createPetscObjects() override;

  int nCompartments_;     ///< number of compartments i.e. motor units
  std::shared_ptr<FieldVariableType> flowPotential_; ///< the direction of fibres
  std::shared_ptr<GradientFieldVariableType> fibreDirection_; ///< the direction of fibres
  std::vector<std::shared_ptr<CellMLFieldVariableType>> transmembranePotential_;  ///< the Vm value for the compartments
  std::vector<std::shared_ptr<CellMLFieldVariableType>> transmembranePotentialNextTimeStep_;  ///< the Vm value for the compartments
  std::vector<std::shared_ptr<CellMLFieldVariableType>> transmembraneIncrementNextTimeStep_;              ///< increment helper variable
  std::vector<std::shared_ptr<FieldVariableType>> ionicCurrentNextTimestep_;  ///< -1/Cm I_ion(Vm) for the next time step
  std::shared_ptr<FieldVariableType> extraCellularPotential_;  ///< the phi_e value which is the extra-cellular potential
};

} // namespace Data

#include "data_management/multidomain.tpp"
