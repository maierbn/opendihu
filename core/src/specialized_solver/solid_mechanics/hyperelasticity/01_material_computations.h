#pragma once

#include <Python.h>  // has to be the first included header

#include "specialized_solver/solid_mechanics/hyperelasticity/expression_helper.h"
#include "specialized_solver/solid_mechanics/hyperelasticity/00_initialize.h"

namespace SpatialDiscretization
{

/** This class contains all formulas for computation of physical quantities.
  */
template<typename Term = Equation::SolidMechanics::MooneyRivlinIncompressible3D, bool withLargeOutput=true, typename MeshType = Mesh::StructuredDeformableOfDimension<3>, int nDisplacementComponents = 3>
class HyperelasticityMaterialComputations :
  public HyperelasticityInitialize<Term,withLargeOutput,MeshType,nDisplacementComponents>
{
public:
  typedef FunctionSpace::FunctionSpace<MeshType, BasisFunction::LagrangeOfOrder<2>> DisplacementsFunctionSpace;
  typedef FunctionSpace::FunctionSpace<MeshType, BasisFunction::LagrangeOfOrder<1>> PressureFunctionSpace;
  typedef FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<1>, BasisFunction::LagrangeOfOrder<1>> FiberFunctionSpace;
  typedef DisplacementsFunctionSpace FunctionSpace;

  typedef Data::QuasiStaticHyperelasticity<PressureFunctionSpace,DisplacementsFunctionSpace,Term,withLargeOutput,Term> Data;
  typedef ::Data::QuasiStaticHyperelasticityPressureOutput<PressureFunctionSpace> PressureDataCopy;

  typedef FieldVariable::FieldVariable<DisplacementsFunctionSpace,3> DisplacementsFieldVariableType;
  typedef FieldVariable::FieldVariable<PressureFunctionSpace,1> PressureFieldVariableType;
  typedef FieldVariable::FieldVariable<DisplacementsFunctionSpace,6> StressFieldVariableType;

  typedef PartitionedPetscVecForHyperelasticity<DisplacementsFunctionSpace,PressureFunctionSpace,Term,nDisplacementComponents> VecHyperelasticity;
  typedef PartitionedPetscMatForHyperelasticity<DisplacementsFunctionSpace,PressureFunctionSpace,Term,nDisplacementComponents> MatHyperelasticity;

  //! constructor
  HyperelasticityMaterialComputations(DihuContext context, std::string settingsKey = "HyperelasticitySolver");

  //! not relevant, as it does nothing. Contains a lot of commented out debugging code that is helpful to debug the analytic jacobian matrix
  template<typename double_v_t>
  void materialTesting(const double_v_t pressure,                            //< [in] pressure value p
                       const Tensor2<3,double_v_t> &rightCauchyGreen,        //< [in] C
                       const Tensor2<3,double_v_t> &inverseRightCauchyGreen, //< [in] C^{-1}
                       const std::array<double_v_t,5> reducedInvariants,     //< [in] the reduced invariants Ibar_1, Ibar_2
                       const double_v_t deformationGradientDeterminant,      //< [in] J = det(F)
                       VecD<3,double_v_t> fiberDirection,                    //< [in] a0, direction of fibers
                       Tensor2<3,double_v_t> &fictitiousPK2Stress,           //< [in] Sbar, the fictitious 2nd Piola-Kirchhoff stress tensor
                       Tensor2<3,double_v_t> &pk2StressIsochoric             //< [in] S_iso, the isochoric part of the 2nd Piola-Kirchhoff stress tensor
                      );

  //! compute δWint from given displacements, this is a wrapper to materialComputeInternalVirtualWork
  bool materialComputeInternalVirtualWork(
    std::shared_ptr<VecHyperelasticity> displacements,
    std::shared_ptr<VecHyperelasticity> internalVirtualWork
  );

  //! solve the dynamic hyperelastic problem, using a load stepping
  //! output internalVirtualWork, externalVirtualWorkDead and accelerationTerm which is the acceleration contribution to δW
  void solveDynamicProblem(std::shared_ptr<VecHyperelasticity> displacementsVelocitiesPressure, bool isFirstTimeStep,
                           Vec internalVirtualWork, Vec &externalVirtualWorkDead, Vec accelerationTerm, bool withOutputWritersEnabled = true);

 void solveQuasistaticProblem(std::shared_ptr<VecHyperelasticity> displacementsVelocitiesPressure, bool isFirstTimeStep,
                           bool withOutputWritersEnabled = true);
  //! compute u from the equation ∂W_int - externalVirtualWork = 0 and J = 1, the result will be in displacements,
  //! the previous value in displacements will be used as initial guess (call zeroEntries() of displacements beforehand, if no initial guess should be provided)
  void solveForDisplacements(
    std::shared_ptr<VecHyperelasticity> externalVirtualWork,
    std::shared_ptr<VecHyperelasticity> displacements
  );

  //! compute the resulting forces and moments at the 2- (bottom) and 2+ (top) surfaces of the volume
  //! @param elements <element_no, isTop> a list of elements that are at bottom and top of the volume and on which the integration over forces and moments is performed, !isTop means bottom element
  void computeBearingForceAndMoment(const std::vector<std::tuple<element_no_t,bool>> &elements,
                                    Vec3 &bearingForceBottom, Vec3 &bearingMomentBottom, Vec3 &bearingForceTop, Vec3 &bearingMomentTop);

protected:

  typedef HyperelasticityInitialize<Term,withLargeOutput,MeshType,nDisplacementComponents> Parent;

  //! compute δW_int, input is in this->data_.displacements() and this->data_.pressure(), output is in solverVariableResidual_
  //! @param communicateGhosts if startGhostManipulation() and finishGhostManipulation() will be called on combinedVecResidual_ inside this method, if set to false, you have to do it manually before and after this method
  //! @return true if computation was successful (i.e. no negative jacobian)
  bool materialComputeInternalVirtualWork(bool communicateGhosts=true);

  //! compute the nonlinear function F(x), x=solverVariableSolution_, F=solverVariableResidual_
  //! solverVariableResidual_[0-2] contains δW_int - δW_ext, solverVariableResidual_[3] contains int_Ω (J-1)*Ψ dV
  //! @param loadFactor: a factor with which the rhs is scaled. This is equivalent to set the body force and traction (neumann bc) to this fraction
  //! @return true if computation was successful (i.e. no negative jacobian)
  bool materialComputeResidual(double loadFactor = 1.0);

  //! compute δW_ext,dead = int_Ω B^L * phi^L * phi^M * δu^M dx + int_∂Ω T^L * phi^L * phi^M * δu^M dS
  void materialComputeExternalVirtualWorkDead();

  //! add the acceleration term that is needed for the dynamic problem to δW_int - δW_ext,dead, and secondly, add the velocity/displacement equation 1/dt (u^(n+1) - u^(n)) - v^(n+1) - v^n = 0
  //! @param communicateGhosts if startGhostManipulation() and finishGhostManipulation() will be called on combinedVecResidual_ inside this method, if set to false, you have to do it manually before and after this method
  void materialAddAccelerationTermAndVelocityEquation(bool communicateGhosts=true);

  //! compute the jacobian of the Newton scheme
  //! @return true if computation was successful (i.e. no negative jacobian)
  bool materialComputeJacobian();

  //! compute the deformation gradient, F inside the current element at position xi, the value of F is still with respect to the reference configuration,
  //! the formula is F_ij = x_i,j = δ_ij + u_i,j
  template<typename double_v_t>
  Tensor2<3,double_v_t> computeDeformationGradient(const std::array<VecD<3,double_v_t>,DisplacementsFunctionSpace::nDofsPerElement()> &displacements,
                                                 const Tensor2<3,double_v_t> &inverseJacobianMaterial,
                                                 const std::array<double,3> xi);

  //! compute the time velocity of the deformation gradient, Fdot inside the current element at position xi, the value of F is still with respect to the reference configuration,
  //! the formula is Fdot_ij = d/dt x_i,j = v_i,j
  template<typename double_v_t>
  Tensor2<3,double_v_t> computeDeformationGradientTimeDerivative(const std::array<VecD<3,double_v_t>,DisplacementsFunctionSpace::nDofsPerElement()> &velocities,
                                                 const Tensor2<3,double_v_t> &inverseJacobianMaterial,
                                                 const std::array<double,3> xi);

  //! compute the right Cauchy Green tensor, C = F^T*F. This is simply a matrix matrix multiplication
  template<typename double_v_t>
  Tensor2<3,double_v_t> computeRightCauchyGreenTensor(const Tensor2<3,double_v_t> &deformationGradient);

  //! compute the invariants, [I1,I2,I3,I4,I5] where I1 = tr(C), I2 = 1/2 * (tr(C)^2 - tr(C^2)), I3 = det(C), I4 = a0•C a0, I5 = a0•C^2 a0
  template<typename double_v_t>
  std::array<double_v_t,5> computeInvariants(const Tensor2<3,double_v_t> &rightCauchyGreen, const double_v_t rightCauchyGreenDeterminant,
                                           const VecD<3,double_v_t> fiberDirection);

  //! compute the reduced invariants, [Ibar1, Ibar2]
  template<typename double_v_t>
  std::array<double_v_t,5> computeReducedInvariants(const std::array<double_v_t,5> invariants, const double_v_t deformationGradientDeterminant);

  //! compute the 2nd Piola-Kirchhoff pressure, S
  template<typename double_v_t>
  Tensor2<3,double_v_t>
  computePK2Stress(double_v_t &pressure,                                  //< [in/out] pressure value p
                   const Tensor2<3,double_v_t> &rightCauchyGreen,         //< [in] C
                   const Tensor2<3,double_v_t> &inverseRightCauchyGreen,  //< [in] C^{-1}
                   const std::array<double_v_t,5> invariants,             //< [in] the strain invariants I_1, ..., I_5
                   const std::array<double_v_t,5> reducedInvariants,      //< [in] the reduced invariants Ibar_1, ..., Ibar_5
                   const double_v_t deformationGradientDeterminant,       //< [in] J = det(F)
                   VecD<3,double_v_t> fiberDirection,                     //< [in] a0, direction of fibers
                   dof_no_v_t elementNoLocalv,                            //< [in] the current element nos (simd vector) with unused entries set to -1, needed only as mask which entries to discard
                   Tensor2<3,double_v_t> &fictitiousPK2Stress,            //< [out] Sbar, the fictitious 2nd Piola-Kirchhoff stress tensor
                   Tensor2<3,double_v_t> &pk2StressIsochoric              //< [out] S_iso, the isochoric part of the 2nd Piola-Kirchhoff stress tensor
                  );

  //! compute the PK2 stress and the deformation gradient at every node and set value in data, for output
  void computePK2StressField();

  //! compute the material elasticity tensor
  template<typename double_v_t>
  void computeElasticityTensor(const Tensor2<3,double_v_t> &rightCauchyGreen,
                               const Tensor2<3,double_v_t> &inverseRightCauchyGreen, double_v_t deformationGradientDeterminant, double_v_t pressure,
                               std::array<double_v_t,5> invariants, std::array<double_v_t,5> reducedInvariants, const Tensor2<3,double_v_t> &fictitiousPK2Stress,
                               const Tensor2<3,double_v_t> &pk2StressIsochoric, VecD<3,double_v_t> fiberDirection,
                               Tensor4<3,double_v_t> &fictitiousElasticityTensor, Tensor4<3,double_v_t> &elasticityTensorIso,
                               Tensor4<3,double_v_t> &elasticityTensor);

  //! compute P : Sbar
  template<typename double_v_t>
  Tensor2<3,double_v_t> computePSbar(const Tensor2<3,double_v_t> &fictitiousPK2Stress, const Tensor2<3,double_v_t> &rightCauchyGreen);

  //! compute the value of Sbar : C
  template<typename double_v_t>
  double computeSbarC(const Tensor2<3,double_v_t> &Sbar, const Tensor2<3,double_v_t> &C);

  //! use Petsc to solve the nonlinear equation using the SNES solver
  virtual void nonlinearSolve() = 0;

  //! do some steps after nonlinearSolve(): close log file, copy the solution values back to this->data_.displacements() and this->data.pressure(),
  //! compute the PK2 stress at every node, update the geometry field by the new displacements, dump files containing rhs and system matrix
  virtual void postprocessSolution() = 0;

  // use variables from the base class, this avoids writing "this->" in front of those members
  using Parent::displacementsFunctionSpace_;          //< the function space with quadratic Lagrange basis functions, used for discretization of displacements
  using Parent::pressureFunctionSpace_;               //< the function space with linear Lagrange basis functions, used for discretization of pressure

  using Parent::solverVariableResidual_;              //< PETSc Vec to store the residual, equal to combinedVecResidual_->valuesGlobal()
  using Parent::solverVariableSolution_;              //< PETSc Vec to store the solution, equal to combinedVecSolution_->valuesGlobal()

  using Parent::combinedVecResidual_;                 //< the Vec for the residual and result of the nonlinear function
  using Parent::combinedVecSolution_;                 //< the Vec for the solution, combined means that ux,uy,uz and p components are combined in one vector
  using Parent::combinedVecExternalVirtualWorkDead_;  //< the Vec for the external virtual work part that does not change with u, δW_ext,dead
  using Parent::combinedMatrixJacobian_;              //< single jacobian matrix
  using Parent::combinedMatrixAdditionalNumericJacobian_;   //< only used when both analytic and numeric jacobians are computed, then this holds the numeric jacobian

  using Parent::externalVirtualWorkDead_;             //< the external virtual work resulting from the traction, this is a dead load, i.e. it does not change during deformation
  using Parent::getString;                            //< function to get a string representation of the values for debugging output
  using Parent::setUVP;                               //< function to copy the values of the vector x which contains (u and p) or (u,v and p) values to this->data_.displacements(), this->data_.velocities() and this->data_.pressure();
};

}  // namespace

#include "specialized_solver/solid_mechanics/hyperelasticity/01_material_computations.tpp"
#include "specialized_solver/solid_mechanics/hyperelasticity/01_material_computations_auxiliary.tpp"
#include "specialized_solver/solid_mechanics/hyperelasticity/01_material_computations_elasticity_tensor.tpp"
#include "specialized_solver/solid_mechanics/hyperelasticity/01_material_computations_stress.tpp"
#include "specialized_solver/solid_mechanics/hyperelasticity/01_material_computations_wrappers.tpp"
#include "specialized_solver/solid_mechanics/hyperelasticity/01_material_testing.tpp"
