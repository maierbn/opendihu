#pragma once

#include <Python.h>  // has to be the first included header
#include <petscmat.h>
#include <memory>
#include <type_traits>

#include "data_management/data.h"
#include "field_variable/field_variable.h"

namespace Data
{

/**  The datastructures used for quasi static hyperelasticity solver.
 *   Template arguments are the two function spaces of the mixed formulation.
 *   Typically, PressureFunctionSpace is linear and DisplacementsFunctionSpace is quadratic.
  */
template<typename PressureFunctionSpace, typename DisplacementsFunctionSpace, typename Term, bool withLargeOutput>
class QuasiStaticHyperelasticityBase :
  public Data<DisplacementsFunctionSpace>
{
public:

  typedef FieldVariable::FieldVariable<DisplacementsFunctionSpace,3> DisplacementsFieldVariableType;
  typedef FieldVariable::FieldVariable<PressureFunctionSpace,3> DisplacementsLinearFieldVariableType;
  typedef FieldVariable::FieldVariable<PressureFunctionSpace,1> PressureFieldVariableType;
  typedef FieldVariable::FieldVariable<DisplacementsFunctionSpace,6> StressFieldVariableType;     // Voigt notation
  typedef FieldVariable::FieldVariable<DisplacementsFunctionSpace,9> DeformationGradientFieldVariableType;  // row-major

  //! define the type of output connection variables, i.e. the values that will be transferred if the solver is part of a splitting or coupling scheme
  typedef SlotConnectorData<DisplacementsFunctionSpace,3> SlotConnectorDataType;

  //! constructor
  QuasiStaticHyperelasticityBase(DihuContext context);

  //! field variable of u (or u^(n+1) in the dynamic problem)
  std::shared_ptr<DisplacementsFieldVariableType> displacements();

  //! field variable of u^(n) for dynamic problem
  std::shared_ptr<DisplacementsFieldVariableType> displacementsPreviousTimestep();

  //! field variable of v (or v^(n+1) in the dynamic problem)
  std::shared_ptr<DisplacementsFieldVariableType> velocities();

  //! field variable of v^(n) for dynamic problem
  std::shared_ptr<DisplacementsFieldVariableType> velocitiesPreviousTimestep();

  //! field variable of p (or p^(n+1) in the dynamic problem)
  std::shared_ptr<PressureFieldVariableType> pressure();

  //! field variable of p
  std::shared_ptr<PressureFieldVariableType> pressurePreviousTimestep();

  //! field variable of reference geometry
  std::shared_ptr<DisplacementsFieldVariableType> geometryReference();

  //! field variable displacements u but on the linear mesh
  std::shared_ptr<DisplacementsLinearFieldVariableType> displacementsLinearMesh();

  //! field variable velocities v, but on the linear mesh
  std::shared_ptr<DisplacementsLinearFieldVariableType> velocitiesLinearMesh();

  //! field variable of F, all 9 values in row-major ordering
  std::shared_ptr<DeformationGradientFieldVariableType> deformationGradient();

  //! field variable of Fdot, all 9 values in row-major ordering
  std::shared_ptr<DeformationGradientFieldVariableType> deformationGradientTimeDerivative();

  //! field variable of S (Voigt notation sxx, syy, szz, sxy, syz, sxz)
  std::shared_ptr<StressFieldVariableType> pK2Stress();

  //! field variable of S_act
  std::shared_ptr<StressFieldVariableType> activePK2Stress();

  //! field variable of fiber direction
  std::shared_ptr<DisplacementsFieldVariableType> &fiberDirection();

  //! traction in current configuration, for z- and z+ surfaces
  std::shared_ptr<DisplacementsFieldVariableType> traction();

  //! traction in reference configuration, for z- and z+ surfaces
  std::shared_ptr<DisplacementsFieldVariableType> materialTraction();

  //! return the object that will be used to transfer values between solvers, in this case this includes only Vm
  std::shared_ptr<SlotConnectorDataType> getSlotConnectorData();

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
  //! \param scalingFactor factor with which to scale the displacements
  //! \param updateLinearVariables if the variables of the linear discretized mesh should be update from the quadratic mesh, this is only needed if the pressure output writer is used
  void updateGeometry(double scalingFactor = 1.0, bool updateLinearVariables = true);

  //! update the copy of the reference geometry field from the current geometry field of the displacementsFunctionSpace
  void updateReferenceGeometry();

protected:

  //! initializes the vectors with size
  void createPetscObjects() override;

  //! compute the field variable pK1Stress_ as P=F*S
  void computePk1Stress();

  std::shared_ptr<PressureFunctionSpace> pressureFunctionSpace_;            //< function space object that discretizes the pressure field variable
  std::shared_ptr<DisplacementsFunctionSpace> displacementsFunctionSpace_;  //< function space object that discretizes the displacements field variable

  std::shared_ptr<DisplacementsFieldVariableType> geometryReference_;       //< the reference configuration geometry
  std::shared_ptr<DisplacementsLinearFieldVariableType> geometryReferenceLinearMesh_;            //< the reference configuration geometry in the pressure function space (linear mesh)

  std::shared_ptr<DisplacementsFieldVariableType> displacements_;                 //< u^(n+1) or u, the displacements
  std::shared_ptr<DisplacementsFieldVariableType> displacementsPreviousTimestep_; //< u^(n), the displacements of the previous timestep in the dynamic case
  std::shared_ptr<DisplacementsFieldVariableType> velocities_;                    //< v^(n+1) or v, the velocities
  std::shared_ptr<DisplacementsFieldVariableType> velocitiesPreviousTimestep_;    //< v^(n), the velocities of the previous timestep in the dynamic case
  std::shared_ptr<PressureFieldVariableType> pressure_;                           //< p^(n+1) for dynamic case or p for static case, the pressure variable
  std::shared_ptr<PressureFieldVariableType> pressurePreviousTimestep_;           //< p^(n), the pressure variable
  std::shared_ptr<StressFieldVariableType> pK2Stress_;                            //< the symmetric PK2 stress tensor in Voigt notation
  std::shared_ptr<StressFieldVariableType> activePK2Stress_;                      //< the symmetric PK2 stress tensor of the active contribution in Voigt notation
  std::shared_ptr<DeformationGradientFieldVariableType> deformationGradient_;     //< the deformation gradient, F, all 9 values in row-major ordering
  std::shared_ptr<DeformationGradientFieldVariableType> deformationGradientTimeDerivative_;     //< the time derivative of the deformation gradient, \dot{F}, all 9 values in row-major ordering
  std::shared_ptr<DisplacementsLinearFieldVariableType> displacementsLinearMesh_; //< the displacements u, but on the linear mesh, not the quadratic. This is an internal helper field
  std::shared_ptr<DisplacementsLinearFieldVariableType> velocitiesLinearMesh_;    //< the velocities v, but on the linear mesh, not the quadratic. This is an internal helper field
  std::shared_ptr<DisplacementsFieldVariableType> fiberDirection_;                //< interpolated direction of fibers
  std::shared_ptr<DisplacementsFieldVariableType> traction_;                      //< t, the traction in current configuration
  std::shared_ptr<DisplacementsFieldVariableType> materialTraction_;              //< T, the traction in reference configuration
  std::shared_ptr<DeformationGradientFieldVariableType> pK1Stress_;               //< the unsymmetric PK1 stress tensor P=FS, this variable is only used internally for the output files
  std::shared_ptr<DeformationGradientFieldVariableType> cauchyStress_;               //< the unsymmetric Cauchy stress tensor σ=J^-1 P F^T, this variable is only used internally for the output files
  std::shared_ptr<FieldVariable::FieldVariable<DisplacementsFunctionSpace,1>> deformationGradientDeterminant_;  //< the determinant of the deformation gradient, J=det F

  std::shared_ptr<SlotConnectorDataType> slotConnectorData_;    //< the object that stores all components of field variables that will be transferred to other solvers
};

/** Helper class that outputs the field variables for the output writer.
 *  Depending on the Term if it uses fiberDirection information, also output a fiber direction field.
 *  The normal isotropic Mooney-Rivlin thus has no fiber direction output.
 *
 *  withLargeOutput = false, Term::usesFiberDirection = false
 */
template<typename PressureFunctionSpace, typename DisplacementsFunctionSpace, typename Term, bool withLargeOutput=false, typename DummyForTraits=Term>
class QuasiStaticHyperelasticity :
  public QuasiStaticHyperelasticityBase<PressureFunctionSpace, DisplacementsFunctionSpace, Term, withLargeOutput>
{
public:
  using QuasiStaticHyperelasticityBase<PressureFunctionSpace, DisplacementsFunctionSpace, Term, withLargeOutput>::QuasiStaticHyperelasticityBase;

  typedef FieldVariable::FieldVariable<DisplacementsFunctionSpace,3> DisplacementsFieldVariableType;
  typedef FieldVariable::FieldVariable<DisplacementsFunctionSpace,6> StressFieldVariableType;

  //! field variables that will be output by outputWriters
  //type to use if there is no fiber direction field variable
  typedef std::tuple<
    std::shared_ptr<DisplacementsFieldVariableType>,  // current geometry field
    std::shared_ptr<DisplacementsFieldVariableType>,  // displacements_
    std::shared_ptr<DisplacementsFieldVariableType>,  // velocities_
    std::shared_ptr<DisplacementsFieldVariableType>,  // material traction
    std::shared_ptr<StressFieldVariableType>         // pK2Stress_
  >
  FieldVariablesForOutputWriter;

  //! get pointers to all field variables that can be written by output writers
  FieldVariablesForOutputWriter getFieldVariablesForOutputWriter();
};

// withLargeOutput = true, Term::usesFiberDirection = false
template<typename PressureFunctionSpace, typename DisplacementsFunctionSpace, typename Term>
class QuasiStaticHyperelasticity<PressureFunctionSpace, DisplacementsFunctionSpace, Term, true, std::enable_if_t<!Term::usesFiberDirection,Term>> :
  public QuasiStaticHyperelasticityBase<PressureFunctionSpace, DisplacementsFunctionSpace, Term, true>
{
public:
  using QuasiStaticHyperelasticityBase<PressureFunctionSpace, DisplacementsFunctionSpace, Term, true>::QuasiStaticHyperelasticityBase;

  typedef FieldVariable::FieldVariable<DisplacementsFunctionSpace,3> DisplacementsFieldVariableType;
  typedef FieldVariable::FieldVariable<DisplacementsFunctionSpace,6> StressFieldVariableType;
  typedef FieldVariable::FieldVariable<DisplacementsFunctionSpace,9> DeformationGradientFieldVariableType;  // row-major

  //! field variables that will be output by outputWriters
  //type to use if there is no fiber direction field variable
  typedef std::tuple<
    std::shared_ptr<DisplacementsFieldVariableType>,  // current geometry field
    std::shared_ptr<DisplacementsFieldVariableType>,  // displacements_
    std::shared_ptr<DisplacementsFieldVariableType>,  // velocities_
    std::shared_ptr<DisplacementsFieldVariableType>,  // traction
    std::shared_ptr<DisplacementsFieldVariableType>,  // material traction
    std::shared_ptr<StressFieldVariableType>,         // pK2Stress_
    std::shared_ptr<DeformationGradientFieldVariableType>,  // deformationGradient_
    std::shared_ptr<DeformationGradientFieldVariableType>,  // deformationGradientTimeDerivative_
    std::shared_ptr<DeformationGradientFieldVariableType>,  // pK1Stress_
    std::shared_ptr<DeformationGradientFieldVariableType>,   // Cauchy stress
    std::shared_ptr<FieldVariable::FieldVariable<DisplacementsFunctionSpace,1>>   // determinant of the material deformation gradient
  >
  FieldVariablesForOutputWriter;

  //! get pointers to all field variables that can be written by output writers
  FieldVariablesForOutputWriter getFieldVariablesForOutputWriter();
};

// withLargeOutput = false, Term::usesFiberDirection = true
template<typename PressureFunctionSpace, typename DisplacementsFunctionSpace, typename Term>
class QuasiStaticHyperelasticity<PressureFunctionSpace, DisplacementsFunctionSpace, Term, false, std::enable_if_t<Term::usesFiberDirection,Term>> :
  public QuasiStaticHyperelasticityBase<PressureFunctionSpace, DisplacementsFunctionSpace, Term, false>
{
public:
  using QuasiStaticHyperelasticityBase<PressureFunctionSpace, DisplacementsFunctionSpace, Term, false>::QuasiStaticHyperelasticityBase;

  typedef FieldVariable::FieldVariable<DisplacementsFunctionSpace,3> DisplacementsFieldVariableType;
  typedef FieldVariable::FieldVariable<DisplacementsFunctionSpace,6> StressFieldVariableType;

  //! field variables that will be output by outputWriters
  // type to use if we have a fiber direction field variable
  typedef std::tuple<
    std::shared_ptr<DisplacementsFieldVariableType>,  // current geometry field
    std::shared_ptr<DisplacementsFieldVariableType>,  // displacements_
    std::shared_ptr<DisplacementsFieldVariableType>,  // velocities_
    std::shared_ptr<StressFieldVariableType>,         // pK2Stress_
    std::shared_ptr<StressFieldVariableType>,         // activePK2Stress_
    std::shared_ptr<DisplacementsFieldVariableType>,  // fiber direction
    std::shared_ptr<DisplacementsFieldVariableType>   // material traction
  >
  FieldVariablesForOutputWriter;

  //! get pointers to all field variables that can be written by output writers
  FieldVariablesForOutputWriter getFieldVariablesForOutputWriter();
};

// withLargeOutput = true, Term::usesFiberDirection = true
template<typename PressureFunctionSpace, typename DisplacementsFunctionSpace, typename Term>
class QuasiStaticHyperelasticity<PressureFunctionSpace, DisplacementsFunctionSpace, Term, true, std::enable_if_t<Term::usesFiberDirection,Term>> :
  public QuasiStaticHyperelasticityBase<PressureFunctionSpace, DisplacementsFunctionSpace, Term, true>
{
public:
  using QuasiStaticHyperelasticityBase<PressureFunctionSpace, DisplacementsFunctionSpace, Term, true>::QuasiStaticHyperelasticityBase;

  typedef FieldVariable::FieldVariable<DisplacementsFunctionSpace,3> DisplacementsFieldVariableType;
  typedef FieldVariable::FieldVariable<DisplacementsFunctionSpace,6> StressFieldVariableType;
  typedef FieldVariable::FieldVariable<DisplacementsFunctionSpace,9> DeformationGradientFieldVariableType;  // row-major

  //! field variables that will be output by outputWriters
  // type to use if we have a fiber direction field variable
  typedef std::tuple<
    std::shared_ptr<DisplacementsFieldVariableType>,  // current geometry field
    std::shared_ptr<DisplacementsFieldVariableType>,  // displacements_
    std::shared_ptr<DisplacementsFieldVariableType>,  // velocities_
    std::shared_ptr<StressFieldVariableType>,         // pK2Stress_
    std::shared_ptr<StressFieldVariableType>,         // activePK2Stress_
    std::shared_ptr<DisplacementsFieldVariableType>,  // fiber direction
    std::shared_ptr<DisplacementsFieldVariableType>,  // traction
    std::shared_ptr<DisplacementsFieldVariableType>,  // material traction
    std::shared_ptr<DeformationGradientFieldVariableType>,  // deformationGradient_
    std::shared_ptr<DeformationGradientFieldVariableType>,  // deformationGradientTimeDerivative_
    std::shared_ptr<DeformationGradientFieldVariableType>,  // pK1Stress_
    std::shared_ptr<DeformationGradientFieldVariableType>,  // Cauchy stress
    std::shared_ptr<FieldVariable::FieldVariable<DisplacementsFunctionSpace,1>>   // determinant of the material deformation gradient
  >
  FieldVariablesForOutputWriter;

  //! get pointers to all field variables that can be written by output writers
  FieldVariablesForOutputWriter getFieldVariablesForOutputWriter();
};

// ------------------------------------------------------------------------------------------------------
/** This is a helper class that stores copies of the pressure variables.
 *  With the normal class, QuasiStaticHyperelasticity only variables on the DisplacementsFunctionSpace can be written.
 *  Using an object of this class, the pressure values can be written.
 *  Use the following, where data_ is an object of QuasiStaticHyperelasticityPressureOutput.
 *
 *  outputWriterManager_.writeOutput(data_, 0, 0);
 */
template<typename PressureFunctionSpace>
class QuasiStaticHyperelasticityPressureOutput :
  public Data<PressureFunctionSpace>
{
public:
  typedef FieldVariable::FieldVariable<PressureFunctionSpace,3> DisplacementsLinearFieldVariableType;
  typedef FieldVariable::FieldVariable<PressureFunctionSpace,1> PressureFieldVariableType;

  //! constructor of base class
  using Data<PressureFunctionSpace>::Data;

  // code to output pressure field variable (it is not possible to output both displacements and pressure, because the function spaces are different)
  typedef std::tuple<
      std::shared_ptr<DisplacementsLinearFieldVariableType>,  // current linear geometry field
      std::shared_ptr<DisplacementsLinearFieldVariableType>,  // displacements in linear function space
      std::shared_ptr<DisplacementsLinearFieldVariableType>,  // velocities in linear function space
      std::shared_ptr<PressureFieldVariableType>  // pressure
    >
  FieldVariablesForOutputWriter;

  //! get pointers to all field variables that can be written by output writers
  FieldVariablesForOutputWriter getFieldVariablesForOutputWriter();

  //! initialize the internal pressure and displacements variables
  void initialize(std::shared_ptr<PressureFieldVariableType> pressure,
                  std::shared_ptr<DisplacementsLinearFieldVariableType> displacementsLinearMesh,
                  std::shared_ptr<DisplacementsLinearFieldVariableType> velocitiesLinearMesh
                 );

private:

  //! this class does not create any petsc objects
  void createPetscObjects(){};

  std::shared_ptr<PressureFieldVariableType> pressure_;     //<  p, the pressure variable
  std::shared_ptr<DisplacementsLinearFieldVariableType> displacementsLinearMesh_;     //<  the displacements u, but on the linear mesh not the quadratic. This is an internal helper field
  std::shared_ptr<DisplacementsLinearFieldVariableType> velocitiesLinearMesh_;     //<  the velocities v, but on the linear mesh not the quadratic. This is an internal helper field
};

} // namespace Data

#include "data_management/specialized_solver/hyperelasticity_solver.tpp"
