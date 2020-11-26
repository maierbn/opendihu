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
class Multidomain : public Data<FunctionSpaceType>
{
public:

  typedef FieldVariable::FieldVariable<FunctionSpaceType,1> FieldVariableType;
  typedef FieldVariable::FieldVariable<FunctionSpaceType,3> GradientFieldVariableType;
  typedef std::tuple<
    std::shared_ptr<std::vector<std::shared_ptr<SlotConnectorData<FunctionSpaceType,1>>>>,
    std::shared_ptr<SlotConnectorData<FunctionSpaceType,1>>
  > SlotConnectorDataType;   // contains (V_mk^(i), V_mk^(i+1), activeStress_k for every compartment k), activeStressTotal_

  //! constructor
  Multidomain(DihuContext context);

  //! return a reference to the rhs summand vector which is needed to apply the boundary conditions, the PETSc Vec can be obtained via fieldVariable->valuesGlobal()
  std::shared_ptr<GradientFieldVariableType> fiberDirection();

  //! return the field variable of the potential for the Laplace potential flow problem
  std::shared_ptr<FieldVariableType> flowPotential();

  //! return the extra-cellular potential field variable
  std::shared_ptr<FieldVariableType> extraCellularPotential();

  //! return the transmembrane potential (Vm) field variable
  std::shared_ptr<FieldVariableType> transmembranePotential(int compartmentNo);

  //! return the solution vector of the transmembrane potential (Vm) field variable
  std::shared_ptr<FieldVariableType> transmembranePotentialSolution(int compartmentNo);

  //! return the transmembrane potential (Vm) field variable as vector for all compartments
  std::vector<std::shared_ptr<FieldVariableType>> transmembranePotential();

  //! return the relative factor f_r of the given compartment, at each point
  std::shared_ptr<FieldVariableType> compartmentRelativeFactor(int compartmentNo);

  //! return the relative factor for phi_e, (1 + sum_k f_r^k), at each point
  std::shared_ptr<FieldVariableType> relativeFactorTotal();

  //! a field variable with constant value of zero, needed for the nested rhs vector
  std::shared_ptr<FieldVariableType> zero();

  //! return the field variable of the active stress for the given compartment
  std::shared_ptr<FieldVariableType> activeStress(int compartmentNo);

  //! return the field variable of the total active stress
  std::shared_ptr<FieldVariableType> activeStressTotal();

  //! initialize and set nCompartments_
  void initialize(int nCompartments);

  //! print all stored data to stdout
  void print();

  //! get the output connection da
  std::shared_ptr<SlotConnectorDataType> getSlotConnectorData();


  //! field variables that will be output by outputWriters
  typedef std::tuple<
    std::shared_ptr<GradientFieldVariableType>,       // geometry
    std::shared_ptr<GradientFieldVariableType>,       // fiberDirection
    std::shared_ptr<FieldVariableType>,               // solution of laplace potential flow
    std::shared_ptr<FieldVariableType>,               // extra-cellular potential
    std::vector<std::shared_ptr<FieldVariableType>>,  // transmembranePotentials
    std::vector<std::shared_ptr<FieldVariableType>>,  // compartmentRelativeFactors
    std::shared_ptr<FieldVariableType>,               // relativeFactorTotal
    std::shared_ptr<FieldVariableType>                // activeStressTotal
  > FieldVariablesForOutputWriter;

  //! get pointers to all field variables that can be written by output writers
  FieldVariablesForOutputWriter getFieldVariablesForOutputWriter();

private:

  //! initializes the vectors with size
  void createPetscObjects() override;

  int nCompartments_;                                                               //< number of compartments i.e. motor units
  std::shared_ptr<FieldVariableType> flowPotential_;                                //< the direction of fibers
  std::shared_ptr<GradientFieldVariableType> fiberDirection_;                       //< the direction of fibers
  std::vector<std::shared_ptr<FieldVariableType>> transmembranePotentialSolution_;  //< the Vm^(i+1) for the next timestep, this holds the solution in the linear solver which must be different from the rhs vector
  std::vector<std::shared_ptr<FieldVariableType>> transmembranePotential_;          //< the Vm^(i) value (transmembrane potential)
  std::vector<std::shared_ptr<FieldVariableType>> compartmentRelativeFactor_;       //< the relative factor f_r of the given compartment, at each point
  std::vector<std::shared_ptr<FieldVariableType>> activeStress_;                    //< input slot for the active stress
  std::shared_ptr<FieldVariableType> activeStressTotal_;                            //< output slot for the active stress, computed as sum_k f_r^k active_stress_k
  std::shared_ptr<FieldVariableType> relativeFactorTotal_;                          //< relative factor for phi_e, (1 + sum_k f_r^k), at each point
  std::shared_ptr<FieldVariableType> extraCellularPotential_;                       //< the phi_e value which is the extra-cellular potential
  std::shared_ptr<FieldVariableType> zero_;                                         //< a field variable with constant value of zero, needed for the nested rhs vector
  //std::vector<Vec> subvectorsSolution_;                                           //< a vector of the Petsc vecs that are used during computation

  std::shared_ptr<SlotConnectorDataType> slotConnectorData_;                    //< the object that stores all components of field variables that will be transferred to other solvers

};

} // namespace Data

#include "data_management/specialized_solver/multidomain.tpp"
