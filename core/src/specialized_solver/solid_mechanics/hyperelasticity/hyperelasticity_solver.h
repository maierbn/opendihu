#pragma once

#include <Python.h>  // has to be the first included header

#include "interfaces/runnable.h"
#include "function_space/function_space.h"
#include "basis_function/basis_function.h"
#include "data_management/specialized_solver/hyperelasticity_solver.h"
#include "partition/partitioned_petsc_vec/02_partitioned_petsc_vec_for_hyperelasticity.h"
#include "partition/partitioned_petsc_mat/partitioned_petsc_mat_for_hyperelasticity.h"
#include "solver/nonlinear.h"
#include "output_writer/manager.h"
#include "spatial_discretization/boundary_conditions/dirichlet_boundary_conditions.h"
#include "spatial_discretization/boundary_conditions/neumann_boundary_conditions.h"

namespace SpatialDiscretization
{

/** This solver is for the nonlinear finite elasticity problem with Mooney-Rivlin material in 3D.
 *  It implements the static equation ∂W_int(u,p) - ∂W_ext = 0 + incompressibility (for nDisplacementComponents = 3)
 *  and the dynamic equation ∂W_int(u,p) - ∂W_ext,dead + ∂W_ext(v) = 0, dot u = v, incompressibility (see details in doc.pdf) (for nDisplacementComponents = 6)
 *
 * Further numerical improvements for solving the nonlinear system that may be implemented later:
 *   https://github.com/jedbrown/spectral-petsc/blob/master/stokes.C
 *   https://lists.mcs.anl.gov/pipermail/petsc-dev/2008-April/000711.html
 *
 * The implementation of this solver is contained in the following files:
 *
 * - "specialized_solver/solid_mechanics/hyperelasticity/hyperelasticity_solver.tpp":
 * This contains the top-level methods and the initialization of data structures.
 *
 * - "specialized_solver/solid_mechanics/hyperelasticity/material_computations.tpp":
 * This contains the implementation of equations, like the residual Wint-Wext, the stress, the analytic jacobian
 *
 * - "specialized_solver/solid_mechanics/hyperelasticity/material_testing.tpp"
 * This contains methods to test the implementations in material_computations.tpp. It is not used for production.
 *
 * - "specialized_solver/solid_mechanics/hyperelasticity/petsc_callbacks.tpp"
 * This contains plain callback functions that are passed to the PETSc nonlinear solver.
 * They call the HyperelasticitySolver object through methods that are implemented in nonlinear_solve.tpp.
 *
 * - "specialized_solver/solid_mechanics/hyperelasticity/nonlinear_solve.tpp"
 * This contains interfaces that are called by PETSc during the nonlinear solution process.
 *
 * - The material equation is given by the structs in equation/mooney_rivlin_incompressible.h, e.g.
 *    Equation::SolidMechanics::MooneyRivlinIncompressible3D or Equation::SolidMechanics::TransverselyIsotropicMooneyRivlinIncompressibleActive3D
 *
 * The template parameter @param nDisplacementComponents refers to the number of non-pressure components
 * and decides if the static (3) or the dynamic (6) case should be computed.
  */
template<typename Term = Equation::SolidMechanics::MooneyRivlinIncompressible3D, int nDisplacementComponents = 3>
class HyperelasticitySolver :
  public Runnable
{
public:
  typedef FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<3>, BasisFunction::LagrangeOfOrder<2>> DisplacementsFunctionSpace;
  typedef FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<3>, BasisFunction::LagrangeOfOrder<1>> PressureFunctionSpace;
  typedef FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<1>, BasisFunction::LagrangeOfOrder<1>> FiberFunctionSpace;

  typedef Data::QuasiStaticHyperelasticity<PressureFunctionSpace,DisplacementsFunctionSpace,Term,Term> Data;
  typedef ::Data::QuasiStaticHyperelasticityPressureOutput<PressureFunctionSpace> PressureDataCopy;

  typedef FieldVariable::FieldVariable<DisplacementsFunctionSpace,3> DisplacementsFieldVariableType;
  typedef FieldVariable::FieldVariable<PressureFunctionSpace,1> PressureFieldVariableType;
  typedef FieldVariable::FieldVariable<DisplacementsFunctionSpace,6> StressFieldVariableType;

  typedef PartitionedPetscVecForHyperelasticity<DisplacementsFunctionSpace,PressureFunctionSpace,nDisplacementComponents> VecHyperelasticity;
  typedef PartitionedPetscMatForHyperelasticity<DisplacementsFunctionSpace,PressureFunctionSpace,nDisplacementComponents> MatHyperelasticity;

  //! constructor
  HyperelasticitySolver(DihuContext context, std::string settingsKey = "HyperelasticitySolver");

  //! advance simulation by the given time span, data in solution is used, afterwards new data is in solution
  void advanceTimeSpan();

  //! initialize components of the simulation
  void initialize();

  //! dummy method, set endTime as current output time
  void setTimeSpan(double startTime, double endTime);

  //! run the simulation
  void run();

  //! reset state
  void reset();

  //! return the data object
  Data &data();

  //! this evaluates the actual nonlinear function f(x) that should be solved f(x) = 0
  void evaluateNonlinearFunction(Vec x, Vec f);

  //! this evaluates the analytic jacobian for the Newton scheme, or the material stiffness
  void evaluateAnalyticJacobian(Vec x, Mat jac);

  //! copy entries of combined vector x to u and p, if both u and p are nullptr, use this->data_.displacements() and this->data_.pressure(), if only p is nullptr, only copy to u
  void setDisplacementsAndPressureFromCombinedVec(Vec x, std::shared_ptr<DisplacementsFieldVariableType> u = nullptr, std::shared_ptr<PressureFieldVariableType> p = nullptr);

  //! copy entries of combined vector x to u, v and p, if all u,v and p are nullptr, use this->data_.displacements(), this->data_.velocities() and this->data_.pressure(), if only p is nullptr, only copy to u
  void setDisplacementsVelocitiesAndPressureFromCombinedVec(Vec x, std::shared_ptr<DisplacementsFieldVariableType> u = nullptr,
                                                            std::shared_ptr<DisplacementsFieldVariableType> v = nullptr,
                                                            std::shared_ptr<PressureFieldVariableType> p = nullptr);

  //! copy the values of the vector x which contains (u and p) or (u,v and p) values to this->data_.displacements(), this->data_.velocities() and this->data_.pressure();
  //! This simply calls setDisplacementsAndPressureFromCombinedVec or setDisplacementsVelocitiesAndPressureFromCombinedVec depending on static or dynamic problem.
  void setUVP(Vec x);

  //! copy the entries in the combined solution vector to this->data_.displacements() and this->data_.pressure()
  void setSolutionVector();

  //! get a string representation of the values for debugging output
  std::string getString(Vec x);

  //! get the PartitionedPetsVec for the residual and result of the nonlinear function
  std::shared_ptr<VecHyperelasticity> combinedVecResidual();      //< the vector that holds the computed residual

  //! get the PartitionedPetsVec for the solution
  std::shared_ptr<VecHyperelasticity> combinedVecSolution();

  //! output the jacobian matrix for debugging
  void dumpJacobianMatrix(Mat jac);

  //! not used
  void debug();

  //! not relevant, as it does nothing. Contains a lot of commented out debugging code that is helpful to debug the analytic jacobian matrix
  void materialTesting(const double pressure,                         //< [in] pressure value p
                   const Tensor2<3> &rightCauchyGreen,                //< [in] C
                   const Tensor2<3> &inverseRightCauchyGreen,         //< [in] C^{-1}
                   const std::array<double,5> reducedInvariants,      //< [in] the reduced invariants Ibar_1, Ibar_2
                   const double deformationGradientDeterminant,       //< [in] J = det(F)
                   Vec3 fiberDirection,                               //< [in] a0, direction of fibers
                   Tensor2<3> &fictitiousPK2Stress,                   //< [in] Sbar, the fictitious 2nd Piola-Kirchhoff stress tensor
                   Tensor2<3> &pk2StressIsochoric                     //< [in] S_iso, the isochoric part of the 2nd Piola-Kirchhoff stress tensor
                  );

  //! callback after each nonlinear iteration
  void monitorSolvingIteration(SNES snes, PetscInt its, PetscReal norm);

  //! create a new shared_ptr of a PartitionedPetscVecForHyperelasticity
  std::shared_ptr<VecHyperelasticity> createPartitionedPetscVec(std::string name);

  //! create a new shared_ptr of a PartitionedPetscMatForHyperelasticity
  std::shared_ptr<MatHyperelasticity> createPartitionedPetscMat(std::string name);

  //! get the precomputed external virtual work
  Vec externalVirtualWork();

  //! compute δWint from given displacements, this is a wrapper to materialComputeInternalVirtualWork
  void materialComputeInternalVirtualWork(
    std::shared_ptr<VecHyperelasticity> displacements,
    std::shared_ptr<VecHyperelasticity> internalVirtualWork
  );

  //! solve the dynamic hyperelastic problem,
  //! output internalVirtualWork, externalVirtualWorkDead and accelerationTerm which is the acceleration contribution to δW
  void solveDynamicProblem(std::shared_ptr<VecHyperelasticity> displacementsVelocitiesPressure, bool isFirstTimeStep,
                           Vec internalVirtualWork, Vec &externalVirtualWorkDead, Vec accelerationTerm);

  //! compute u from the equation ∂W_int - externalVirtualWork = 0 and J = 1, the result will be in displacements,
  //! the previous value in displacements will be used as initial guess (call zeroEntries() of displacements beforehand, if no initial guess should be provided)
  void solveForDisplacements(
    std::shared_ptr<VecHyperelasticity> externalVirtualWork,
    std::shared_ptr<VecHyperelasticity> displacements
  );

protected:

  //! initialize all Petsc Vec's and Mat's that will be used in the computation
  void initializePetscVariables();

  //! initialize the field variable which stores directions of fibers
  void initializeFiberDirections();

  //! use Petsc to solve the nonlinear equation using the SNES solver
  void nonlinearSolve();

  //! set the solution variable to zero and the initial values
  void initializeSolutionVariable();

  //! set all PETSc callback functions, e.g. for computation jacobian or the nonlinear function itself
  void initializePetscCallbackFunctions();

  //! do some steps after nonlinearSolve(): close log file, copy the solution values back to this->data_.displacements() and this->data.pressure(),
  //! compute the PK2 stress at every node, update the geometry field by the new displacements, dump files containing rhs and system matrix
  void postprocessSolution();

  //! check if the solution satisfies Dirichlet BC and the residual is zero, only for debugging output
  void checkSolution(Vec x);

  //! compute δW_int, input is in this->data_.displacements() and this->data_.pressure(), output is in solverVariableResidual_
  void materialComputeInternalVirtualWork();

  //! compute the nonlinear function F(x), x=solverVariableSolution_, F=solverVariableResidual_
  //! solverVariableResidual_[0-2] contains δW_int - δW_ext, solverVariableResidual_[3] contains int_Ω (J-1)*Ψ dV
  void materialComputeResidual();

  //! compute δW_ext,dead = int_Ω B^L * phi^L * phi^M * δu^M dx + int_∂Ω T^L * phi^L * phi^M * δu^M dS
  void materialComputeExternalVirtualWorkDead();

  //! add the acceleration term that is needed for the dynamic problem to δW_int - δW_ext,dead, and secondly, add the velocity/displacement equation 1/dt (u^(n+1) - u^(n)) - v^(n+1) - v^n = 0
  void materialAddAccelerationTermAndVelocityEquation();

  //! compute the jacobian of the Newton scheme
  void materialComputeJacobian();

  //! compute the deformation gradient, F inside the current element at position xi, the value of F is still with respect to the reference configuration,
  //! the formula is F_ij = x_i,j = δ_ij + u_i,j
  Tensor2<3> computeDeformationGradient(const std::array<Vec3,DisplacementsFunctionSpace::nDofsPerElement()> &displacements,
                                        const Tensor2<3> &inverseJacobianMaterial,
                                        const std::array<double, 3> xi);

  //! compute the time velocity of the deformation gradient, Fdot inside the current element at position xi, the value of F is still with respect to the reference configuration,
  //! the formula is Fdot_ij = d/dt x_i,j = v_i,j
  Tensor2<3> computeDeformationGradientTimeDerivative(const std::array<Vec3,DisplacementsFunctionSpace::nDofsPerElement()> &velocities,
                                        const Tensor2<3> &inverseJacobianMaterial,
                                        const std::array<double, 3> xi);

  //! compute the right Cauchy Green tensor, C = F^T*F. This is simply a matrix matrix multiplication
  Tensor2<3> computeRightCauchyGreenTensor(const Tensor2<3> &deformationGradient);

  //! compute the invariants, [I1,I2,I3,I4,I5] where I1 = tr(C), I2 = 1/2 * (tr(C)^2 - tr(C^2)), I3 = det(C), I4 = a0•C0 a0, I5 = a0•C0^2 a0
  std::array<double,5> computeInvariants(const Tensor2<3> &rightCauchyGreen, const double rightCauchyGreenDeterminant, const Vec3 fiberDirection);

  //! compute the reduced invariants, [Ibar1, Ibar2]
  std::array<double,5> computeReducedInvariants(const std::array<double,5> invariants, const double deformationGradientDeterminant);

  //! compute the 2nd Piola-Kirchhoff pressure, S
  Tensor2<3>
  computePK2Stress(const double pressure,                           //< [in] pressure value p
                   const Tensor2<3> &rightCauchyGreen,                //< [in] C
                   const Tensor2<3> &inverseRightCauchyGreen,         //< [in] C^{-1}
                   const std::array<double,5> reducedInvariants,      //< [in] the reduced invariants Ibar_1, Ibar_2
                   const double deformationGradientDeterminant,       //< [in] J = det(F)
                   Vec3 fiberDirection,                               //< [in] a0, direction of fibers
                   Tensor2<3> &fictitiousPK2Stress,                   //< [out] Sbar, the fictitious 2nd Piola-Kirchhoff stress tensor
                   Tensor2<3> &pk2StressIsochoric                     //< [out] S_iso, the isochoric part of the 2nd Piola-Kirchhoff stress tensor
                  );

  //! compute the PK2 stress and the deformation gradient at every node and set value in data, for output
  void computePK2StressField();

  //! compute the material elasticity tensor
  void computeElasticityTensor(const Tensor2<3> &rightCauchyGreen,
                               const Tensor2<3> &inverseRightCauchyGreen, double deformationGradientDeterminant, double pressure,
                               std::array<double,5> reducedInvariants, const Tensor2<3> &fictitiousPK2Stress, const Tensor2<3> &pk2StressIsochoric,
                               Vec3 fiberDirection, Tensor4<3> &fictitiousElasticityTensor, Tensor4<3> &elasticityTensorIso, Tensor4<3> &elasticityTensor);

  //! compute P : Sbar
  Tensor2<3> computePSbar(const Tensor2<3> &fictitiousPK2Stress, const Tensor2<3> &rightCauchyGreen);

  //! compute the value of Sbar : C
  double computeSbarC(const Tensor2<3> &Sbar, const Tensor2<3> &C);

  DihuContext context_;    ///< object that contains the python config for the current context and the global singletons meshManager and solverManager

  OutputWriter::Manager outputWriterManager_; ///< manager object holding all output writer for displacements based variables
  OutputWriter::Manager outputWriterManagerPressure_; ///< manager object holding all output writer for pressure based variables
  Data data_;                 ///< data object
  PressureDataCopy pressureDataCopy_;   ///< a helper object that is used to write the pressure function space based variables with the output writers
  std::shared_ptr<Solver::Nonlinear> nonlinearSolver_;  ///< the nonlinear solver object that will provide the PETSc SNES context

  std::string durationLogKey_;   ///< key with with the duration of the computation is written to the performance measurement log
  std::shared_ptr<DisplacementsFunctionSpace> displacementsFunctionSpace_;  ///< the function space with quadratic Lagrange basis functions, used for discretization of displacements
  std::shared_ptr<PressureFunctionSpace> pressureFunctionSpace_;  ///< the function space with linear Lagrange basis functions, used for discretization of pressure

  Mat solverMatrixJacobian_;           //< the jacobian matrix for the Newton solver, which in case of nonlinear elasticity is the tangent stiffness matrix
  Mat solverMatrixAdditionalNumericJacobian_;           //< only used when both analytic and numeric jacobians are computed, then this holds the numeric jacobian
  Vec solverVariableResidual_;         //< PETSc Vec to store the residual, equal to combinedVecResidual_->valuesGlobal()
  Vec solverVariableSolution_;         //< PETSc Vec to store the solution, equal to combinedVecSolution_->valuesGlobal()
  Vec zeros_;                          // a solver that contains all zeros, needed to zero the diagonal of the jacobian matrix

  std::shared_ptr<VecHyperelasticity> combinedVecResidual_;      //< the Vec for the residual and result of the nonlinear function
  std::shared_ptr<VecHyperelasticity> combinedVecSolution_;      //< the Vec for the solution, combined means that ux,uy,uz and p components are combined in one vector
  std::shared_ptr<VecHyperelasticity> combinedVecExternalVirtualWorkDead_;      //< the Vec for the external virtual work part that does not change with u, δW_ext,dead

  //std::shared_ptr<PartitionedPetscMat<FunctionSpace::Generic>> combinedMatrixJacobian_;    //< single jacobian matrix, when useNestedMat_ is false
  std::shared_ptr<MatHyperelasticity> combinedMatrixJacobian_;    //< single jacobian matrix
  std::shared_ptr<MatHyperelasticity> combinedMatrixAdditionalNumericJacobian_;   //< only used when both analytic and numeric jacobians are computed, then this holds the numeric jacobian

  Vec externalVirtualWorkDead_;     // the external virtual work resulting from the traction, this is a dead load, i.e. it does not change during deformation

  // settings variables
  bool initialized_;   ///< if this object was already initialized
  PythonConfig specificSettings_;    ///< python object containing the value of the python config dict with corresponding key
  double endTime_;     ///< end time of current time step
  std::ofstream residualNormLogFile_;   ///< ofstream of a log file that will contain the residual norm for each iteration

  std::shared_ptr<DirichletBoundaryConditions<DisplacementsFunctionSpace,nDisplacementComponents>> dirichletBoundaryConditions_ = nullptr;  ///< object that parses Dirichlet boundary conditions and applies them to rhs
  std::shared_ptr<NeumannBoundaryConditions<DisplacementsFunctionSpace,Quadrature::Gauss<3>,3>> neumannBoundaryConditions_ = nullptr;  ///< object that parses Neumann boundary conditions and applies them to the rhs

  std::vector<double> materialParameters_;    ///< material parameters, e.g. c1,c2 for Mooney-Rivlin
  double displacementsScalingFactor_;   ///< factor with which to scale the displacements
  bool dumpDenseMatlabVariables_;      ///< the current vector x, the residual, r and the jacobian, jac should be written
  Vec3 constantBodyForce_;   ///< the constant body force, if given or [0,0,0]
  double timeStepWidth_;      ///< timeStepWidth, only need for the dynamic problem
  double density_;                ///< density, only needed for the dynamic problem
  double lastNorm_;            ///< residual norm of the last iteration in the nonlinear solver
  double secondLastNorm_;            ///< residual norm of the second last iteration in the nonlinear solver

  bool useAnalyticJacobian_;   ///< if the analytically computed Jacobian of the Newton scheme should be used. Theoretically if it is correct, this is the fastest option.
  bool useNumericJacobian_;   ///< if a numerically computed Jacobian should be used, approximated by finite differences
};

}  // namespace

#include "specialized_solver/solid_mechanics/hyperelasticity/hyperelasticity_solver.tpp"
#include "specialized_solver/solid_mechanics/hyperelasticity/material_computations.tpp"
#include "specialized_solver/solid_mechanics/hyperelasticity/material_computations_wrappers.tpp"
#include "specialized_solver/solid_mechanics/hyperelasticity/nonlinear_solve.tpp"
#include "specialized_solver/solid_mechanics/hyperelasticity/material_testing.tpp"
