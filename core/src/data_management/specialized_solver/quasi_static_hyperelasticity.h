#pragma once

#include <Python.h>  // has to be the first included header
#include <petscmat.h>
#include <memory>

#include "data_management/data.h"
#include "field_variable/field_variable.h"

namespace Data
{

/**  The datastructures used for quasi static hyperelasticity solver.
 *   Template arguments are the two function spaces of the mixed formulation.
 *   Typically, PressureFunctionSpace is linear and DisplacementsFunctionSpace is quadratic.
  */
template<typename PressureFunctionSpace, typename DisplacementsFunctionSpace>
class QuasiStaticHyperelasticity :
  public Data<DisplacementsFunctionSpace>
{
public:

  typedef FieldVariable::FieldVariable<DisplacementsFunctionSpace,3> DisplacementsFieldVariableType;
  typedef FieldVariable::FieldVariable<PressureFunctionSpace,3> DisplacementsLinearFieldVariableType;
  typedef FieldVariable::FieldVariable<PressureFunctionSpace,1> PressureFieldVariableType;
  typedef FieldVariable::FieldVariable<DisplacementsFunctionSpace,6> StressFieldVariableType;

  //! constructor
  QuasiStaticHyperelasticity(DihuContext context);

  //! field variable of u
  std::shared_ptr<DisplacementsFieldVariableType> displacements();

  //! field variable of reference geometry
  std::shared_ptr<DisplacementsFieldVariableType> geometryReference();

  //! field variable of p
  std::shared_ptr<PressureFieldVariableType> pressure();

  //! field variable of S
  std::shared_ptr<StressFieldVariableType> pK2Stress();

  //! initialize
  void initialize();

  //! print all stored data to stdout
  void print();

  //! set the function space object that discretizes the pressure field variable
  void setPressureFunctionSpace(std::shared_ptr<PressureFunctionSpace> pressureFunctionSpace);

  //! set the function space object that discretizes the displacements field variable
  void setDisplacementsFunctionSpace(std::shared_ptr<DisplacementsFunctionSpace> displacementsFunctionSpace);

  //! get the displacements function space
  std::shared_ptr<DisplacementsFunctionSpace> displacementsFunctionSpace();

  //! get the pressure function space
  std::shared_ptr<PressureFunctionSpace> pressureFunctionSpace();

  //! add the displacements to the reference geometry to obtain the current geometry, for both function spaces (linear and quadratic)
  void updateGeometry();

  //! field variables that will be output by outputWriters
  typedef std::tuple<
      std::shared_ptr<DisplacementsFieldVariableType>,  // current geometry field
      std::shared_ptr<DisplacementsFieldVariableType>,  // displacements_
      std::shared_ptr<StressFieldVariableType>         // pK2Stress_
    >
   OutputFieldVariables;

  //! get pointers to all field variables that can be written by output writers
  OutputFieldVariables getOutputFieldVariables();

private:

  //! initializes the vectors with size
  void createPetscObjects() override;

  std::shared_ptr<PressureFunctionSpace> pressureFunctionSpace_;            //< function space object that discretizes the pressure field variable
  std::shared_ptr<DisplacementsFunctionSpace> displacementsFunctionSpace_;  //< function space object that discretizes the displacements field variable

  std::shared_ptr<DisplacementsFieldVariableType> geometryReference_;       //< the reference configuration geometry
  std::shared_ptr<DisplacementsLinearFieldVariableType> geometryReferenceLinearMesh_;            //< the reference configuration geometry in the pressure function space (linear mesh)

  std::shared_ptr<DisplacementsFieldVariableType> displacements_;     //< u, the displacements
  std::shared_ptr<PressureFieldVariableType> pressure_;     //<  p, the pressure variable
  std::shared_ptr<StressFieldVariableType> pK2Stress_;     //<  the symmetric PK2 stress tensor in Voigt notation
  std::shared_ptr<DisplacementsLinearFieldVariableType> displacementsLinearMesh_;     //<  the displacements u, but on the linear mesh not the quadratic. This is an internal helper field
};

} // namespace Data

#include "data_management/specialized_solver/quasi_static_hyperelasticity.tpp"
