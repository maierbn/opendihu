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

/**  The datastructures used for MuscleContractionSolver.
  */
template<typename FunctionSpaceType>
class MuscleContractionSolver : public Data<FunctionSpaceType>
{
public:

  //! define the type of a scalar field variable (1 component), for convenience
  typedef FieldVariable::FieldVariable<FunctionSpaceType,1> ScalarFieldVariableType;

  //! define the type of a vector-valued field variable (3 component), for convenience
  typedef FieldVariable::FieldVariable<FunctionSpaceType,3> VectorFieldVariableType;
  typedef FieldVariable::FieldVariable<FunctionSpaceType,6> StressFieldVariableType;     // stress tensor in Voigt notation

  //! define the type of output connection variables, i.e. the values that will be transferred if the solver is part of a splitting or coupling scheme
  //! Two different field variables can be used: they must have the same function space but can have a different number of components. For example, for the CellMLAdapter, there are the "states" and the "algebraics" field variables.
  //! In this example, we use twice "1" as number of components, but you could, e.g. have SlotConnectorData<FunctionSpaceType,3,4>, etc.
  //! For each field variable you can transfer an abritrary subset of their components.
  typedef SlotConnectorData<FunctionSpaceType,1,3> SlotConnectorDataType;

  //! constructor
  MuscleContractionSolver(DihuContext context);

  //! initialize and create all variables
  void initialize();

  //! print all stored data to stdout
  void print();

  //! return a reference to lambda
  std::shared_ptr<ScalarFieldVariableType> lambda();

  //! return a reference to lambdaDot
  std::shared_ptr<ScalarFieldVariableType> lambdaDot();

  //! return a reference to lambdaDot
  std::shared_ptr<ScalarFieldVariableType> gamma();

  //! return a reference to materialTraction_
  std::shared_ptr<VectorFieldVariableType> materialTraction();

  //! return a reference to displacements_
  std::shared_ptr<VectorFieldVariableType> displacements();

  //! return a reference to velocities_
  std::shared_ptr<VectorFieldVariableType> velocities();

  //! return the object that will be used to transfer values between solvers, in this case this includes only Vm
  std::shared_ptr<SlotConnectorDataType> getSlotConnectorData();

  //! field variables that will be output by outputWriters
  typedef std::tuple<
    std::shared_ptr<VectorFieldVariableType>,     // geometry, this always has to be the first field variable, such that the output writer knows the geometry of the mesh
    std::shared_ptr<ScalarFieldVariableType>,     // lambda,
    std::shared_ptr<ScalarFieldVariableType>,     // lambdaDot
    std::shared_ptr<ScalarFieldVariableType>,     // gamma
    std::shared_ptr<VectorFieldVariableType>,     // displacements
    std::shared_ptr<VectorFieldVariableType>,     // velocities
    std::shared_ptr<StressFieldVariableType>,     // pK2Stress_
    std::shared_ptr<StressFieldVariableType>,     // activePK2Stress_
    std::shared_ptr<VectorFieldVariableType>,     // fiber direction
    std::shared_ptr<VectorFieldVariableType>      // material traction
  > FieldVariablesForOutputWriter;

  //! get pointers to all field variables that can be written by output writers
  FieldVariablesForOutputWriter getFieldVariablesForOutputWriter();

  //! set field variables that were created outside of this object
  //! @param setGeometryFieldForTransfer: if the geometry field should be transferred via the connector slots
  void setFieldVariables(std::shared_ptr<MuscleContractionSolver<FunctionSpaceType>::VectorFieldVariableType> displacements,
                  std::shared_ptr<MuscleContractionSolver<FunctionSpaceType>::VectorFieldVariableType> velocities,
                  std::shared_ptr<MuscleContractionSolver<FunctionSpaceType>::StressFieldVariableType> activePK2Stress,
                  std::shared_ptr<MuscleContractionSolver<FunctionSpaceType>::StressFieldVariableType> pK2Stress,
                  std::shared_ptr<MuscleContractionSolver<FunctionSpaceType>::VectorFieldVariableType> fiberDirection,
                  std::shared_ptr<MuscleContractionSolver<FunctionSpaceType>::VectorFieldVariableType> materialTraction,
                  bool setGeometryFieldForTransfer);

private:

  //! create all field variables with their respective sizes, this will be called automatically within initialize by the base class
  void createPetscObjects() override;

  std::shared_ptr<ScalarFieldVariableType> lambda_;           //< the relative sarcomere length
  std::shared_ptr<ScalarFieldVariableType> lambdaDot_;        //< the contraction velocity
  std::shared_ptr<ScalarFieldVariableType> gamma_;            //< the active stress parameter

  std::shared_ptr<VectorFieldVariableType> displacements_;    //< u, the displacements
  std::shared_ptr<VectorFieldVariableType> velocities_;       //< v, the velocities
  std::shared_ptr<StressFieldVariableType> activePK2Stress_;  //< the symmetric PK2 stress tensor of the active contribution in Voigt notation
  std::shared_ptr<StressFieldVariableType> pK2Stress_;        //< the symmetric PK2 stress tensor in Voigt notation
  std::shared_ptr<VectorFieldVariableType> fiberDirection_;   //< direction of fibers at current point
  std::shared_ptr<VectorFieldVariableType> materialTraction_; //< traction in reference configuration for top and bottom faces (z-, z+)

  std::shared_ptr<SlotConnectorDataType> slotConnectorData_;    //< the object that stores all components of field variables that will be transferred to other solvers

  // define all needed field variables or other data
};

} // namespace Data

#include "data_management/specialized_solver/muscle_contraction_solver.tpp"
