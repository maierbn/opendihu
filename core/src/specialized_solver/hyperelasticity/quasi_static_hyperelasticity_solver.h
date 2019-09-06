#pragma once

#include <Python.h>  // has to be the first included header

#include "interfaces/runnable.h"
#include "function_space/function_space.h"
#include "basis_function/basis_function.h"
#include "data_management/specialized_solver/quasi_static_hyperelasticity.h"
#include "partition/partitioned_petsc_vec/02_partitioned_petsc_vec_for_hyperelasticity.h"
#include "partition/partitioned_petsc_mat/partitioned_petsc_mat_for_hyperelasticity.h"
#include "output_writer/manager.h"
#include "spatial_discretization/boundary_conditions/dirichlet_boundary_conditions.h"
#include "spatial_discretization/boundary_conditions/neumann_boundary_conditions.h"

namespace TimeSteppingScheme
{

/** This solver is for the nonlinear finite elasticity problem with Mooney-Rivlin material in 3D.
  */
class QuasiStaticHyperelasticitySolver :
  public Runnable
{
public:
  typedef FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<3>, BasisFunction::LagrangeOfOrder<2>> DisplacementsFunctionSpace;
  typedef FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<3>, BasisFunction::LagrangeOfOrder<1>> PressureFunctionSpace;
  typedef Data::QuasiStaticHyperelasticity<PressureFunctionSpace,DisplacementsFunctionSpace> Data;
  typedef ::Data::QuasiStaticHyperelasticityPressureOutput<PressureFunctionSpace> PressureDataCopy;

  typedef FieldVariable::FieldVariable<DisplacementsFunctionSpace,3> DisplacementsFieldVariableType;
  typedef FieldVariable::FieldVariable<PressureFunctionSpace,1> PressureFieldVariableType;
  typedef FieldVariable::FieldVariable<DisplacementsFunctionSpace,6> StressFieldVariableType;

  //! constructor
  QuasiStaticHyperelasticitySolver(DihuContext context);

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

  //! return whether the underlying discretizableInTime object has a specified mesh type and is not independent of the mesh type
  bool knowsMeshType();

  //! return the data object
  Data &data();

  //! zero rows and columns in jac for which dirichlet values are set, set diagonal to 1
  void applyDirichletBoundaryConditionsInJacobian(Vec x, Mat jac);

  //! set f to x-x0 for values with dirichlet boundary condition
  void applyDirichletBoundaryConditionsInNonlinearFunction(Vec x, Vec f);

  //! apply the Dirichlet BCs in the input vector x, i.e. set all Dirichlet values in x
  void applyDirichletBoundaryConditionsInVector(Vec x);

  //! this evaluates the actual nonlinear function f(x) that should be solved f(x) = 0
  void evaluateNonlinearFunction(Vec x, Vec f);

  //! this evaluates the analytic jacobian for the Newton scheme, or the material stiffness
  void evaluateAnalyticJacobian(Vec x, Mat jac);

  //! if not useNestedMat_, copy entries of combined vector x to this->data_.displacements() and this->data_.pressure()
  void setInputVector(Vec x);

  //! copy the entries in the combined solution vector to this->data_.displacements() and this->data_.pressure()
  void setSolutionVector();

  //! get a string representation of the values for debugging output
  std::string getString(Vec x);

  //! get the PartitionedPetsVec for the residual and result of the nonlinear function
  std::shared_ptr<PartitionedPetscVecForHyperelasticity<DisplacementsFunctionSpace,PressureFunctionSpace>> combinedVecResidual();      //< the

  //! get the PartitionedPetsVec for the solution
  std::shared_ptr<PartitionedPetscVecForHyperelasticity<DisplacementsFunctionSpace,PressureFunctionSpace>> combinedVecSolution();      //< the Vec for the solution

  //! output the jacobian matrix for debugging
  void dumpJacobianMatrix(Mat jac);

  void debug();

  //! callback after each nonlinear iteration
  void monitorSolvingIteration(SNES snes, PetscInt its, PetscReal norm);

protected:

  //! use Petsc to solve the nonlinear equation using the SNES solver
  void nonlinearSolve();

  //! set the solution variable to zero and the initial values
  void initializeSolutionVariable();

  //! check if the solution satisfies Dirichlet BC and the residual is zero, only for debugging output
  void checkSolution(Vec x);

  //! compute the nonlinear function F(x), x=solverVariableSolution_, F=solverVariableResidual_
  //! solverVariableResidual_[0-2] contains δW_int - δW_ext, solverVariableResidual_[3] contains int_Ω (J-1)*Ψ dV
  void materialComputeResidual();

  //! compute the jacobian of the Newton scheme
  void materialComputeJacobian();

  //! compute the deformation gradient, F inside the current element at position xi, the value of F is still with respect to the reference configuration,
  //! the formula is F_ij = x_i,j = δ_ij + u_i,j
  Tensor2<3> computeDeformationGradient(const std::array<Vec3,DisplacementsFunctionSpace::nDofsPerElement()> &displacements,
                                        const Tensor2<3> &inverseJacobianMaterial,
                                        const std::array<double, 3> xi);

  //! compute the right Cauchy Green tensor, C = F^T*F. This is simply a matrix matrix multiplication
  Tensor2<3> computeRightCauchyGreenTensor(const Tensor2<3> &deformationGradient);

  //! compute the invariants, [I1,I2,I3] where I1 = tr(C), I2 = 1/2 * (tr(C)^2 - tr(C^2)), I3 = det(C)
  std::array<double,3> computeInvariants(const Tensor2<3> &rightCauchyGreen, const double rightCauchyGreenDeterminant);

  //! compute the reduced invariants, [Ibar1, Ibar2]
  std::array<double,2> computeReducedInvariants(const std::array<double,3> invariants, const double deformationGradientDeterminant);

  //! compute the 2nd Piola-Kirchhoff pressure, S
  Tensor2<3>
  computePK2Stress(const double pressure,                           //< [in] pressure value p
                   const Tensor2<3> &rightCauchyGreen,                //< [in] C
                   const Tensor2<3> &inverseRightCauchyGreen,         //< [in] C^{-1}
                   const std::array<double,2> reducedInvariants,      //< [in] the reduced invariants Ibar_1, Ibar_2
                   const double deformationGradientDeterminant,       //< [in] J = det(F)
                   Tensor2<3> &fictitiousPK2Stress,                   //< [out] Sbar, the fictitious 2nd Piola-Kirchhoff stress tensor
                   Tensor2<3> &pk2StressIsochoric                     //< [out] S_iso, the isochoric part of the 2nd Piola-Kirchhoff stress tensor
                  );

  //! compute the PK2 stress at every node and set value in data, for output
  void computePK2StressField();

  //! compute the material elasticity tensor
  void computeElasticityTensor(const Tensor2<3> &rightCauchyGreen,
                               const Tensor2<3> &inverseRightCauchyGreen, double deformationGradientDeterminant, double pressure,
                               std::array<double,2> reducedInvariants, const Tensor2<3> &fictitiousPK2Stress, const Tensor2<3> &pk2StressIsochoric,
                               Tensor4<3> &elasticityTensor);


  DihuContext context_;    ///< object that contains the python config for the current context and the global singletons meshManager and solverManager

  OutputWriter::Manager outputWriterManager_; ///< manager object holding all output writer for displacements based variables
  OutputWriter::Manager outputWriterManagerPressure_; ///< manager object holding all output writer for pressure based variables
  Data data_;                 ///< data object
  PressureDataCopy pressureDataCopy_;   ///< a helper object that is used to write the pressure function space based variables with the output writers

  std::string durationLogKey_;   ///< key with with the duration of the computation is written to the performance measurement log
  std::shared_ptr<DisplacementsFunctionSpace> displacementsFunctionSpace_;  ///< the function space with quadratic Lagrange basis functions, used for discretization of displacements
  std::shared_ptr<PressureFunctionSpace> pressureFunctionSpace_;  ///< the function space with linear Lagrange basis functions, used for discretization of pressure

  Mat solverMatrixJacobian_;           //< the jacobian matrix for the Newton solver, which in case of nonlinear elasticity is the tangent stiffness matrix
  Mat solverMatrixAdditionalNumericJacobian_;           //< only used when both analytic and numeric jacobians are computed, then this holds the numeric jacobian
  Vec solverVariableResidual_;         //< nested PETSc Vec to store the residual
  Vec solverVariableSolution_;         //< nested PETSc Vec to store the solution

  // data structures for nested matrices and vectors
  std::array<Mat,16> submatrices_;  // all submatrices of the 4x4 block jacobian matrix, solverMatrixJacobian_
  std::array<Vec,4> subvectorsSolution_;  // all subvectors of the 4 entries vector, solverVariableSolution_
  std::array<Vec,4> subvectorsResidual_;  // all subvectors of the 4 entries vector, solverVariableResidual_
  std::vector<PartitionedPetscMat<DisplacementsFunctionSpace,DisplacementsFunctionSpace>> uMatrix_;    //< upper left 3x3 blocks of matrices in the jacobian for Newton scheme
  std::vector<PartitionedPetscMat<DisplacementsFunctionSpace,PressureFunctionSpace>> upMatrix_;        //< lower left 1x3 blocks of matrices in the jacobian for Newton scheme
  std::vector<PartitionedPetscMat<PressureFunctionSpace,DisplacementsFunctionSpace>> puMatrix_;        //< upper right 3x1 blocks of matrices in the jacobian for Newton scheme
  std::vector<PartitionedPetscMat<PressureFunctionSpace,PressureFunctionSpace>> pMatrix_;                           //< lower left matrix in the jacobian for Newton scheme

  // data structures for combined matrices and vectors
  std::shared_ptr<PartitionedPetscVecForHyperelasticity<DisplacementsFunctionSpace,PressureFunctionSpace>> combinedVecResidual_;      //< the Vec for the residual and result of the nonlinear function
  std::shared_ptr<PartitionedPetscVecForHyperelasticity<DisplacementsFunctionSpace,PressureFunctionSpace>> combinedVecSolution_;      //< the Vec for the solution
  std::shared_ptr<PartitionedPetscVecForHyperelasticity<DisplacementsFunctionSpace,PressureFunctionSpace>> combinedVecExternalVirtualWork_;      //< the Vec for the external virtual work

  //std::shared_ptr<PartitionedPetscMat<FunctionSpace::Generic>> combinedMatrixJacobian_;    //< single jacobian matrix, when useNestedMat_ is false
  std::shared_ptr<PartitionedPetscMatForHyperelasticity<DisplacementsFunctionSpace,PressureFunctionSpace>> combinedMatrixJacobian_;    //< single jacobian matrix, when useNestedMat_ is false
  std::shared_ptr<PartitionedPetscMatForHyperelasticity<DisplacementsFunctionSpace,PressureFunctionSpace>> combinedMatrixAdditionalNumericJacobian_;   //< only used when both analytic and numeric jacobians are computed, then this holds the numeric jacobian

  Vec externalVirtualWork_;     // the external virtual work resulting from the traction, this is a dead load, i.e. it does not change during deformation

  // settings variables
  bool initialized_;   ///< if this object was already initialized
  PythonConfig specificSettings_;    ///< python object containing the value of the python config dict with corresponding key
  double endTime_;     ///< end time of current time step
  std::ofstream residualNormLogFile_;   ///< ofstream of a log file that will contain the residual norm for each iteration

  std::shared_ptr<SpatialDiscretization::DirichletBoundaryConditions<DisplacementsFunctionSpace,3>> dirichletBoundaryConditions_ = nullptr;  ///< object that parses Dirichlet boundary conditions and applies them to rhs
  std::shared_ptr<SpatialDiscretization::NeumannBoundaryConditions<DisplacementsFunctionSpace,Quadrature::Gauss<3>,3>> neumannBoundaryConditions_ = nullptr;  ///< object that parses Neumann boundary conditions and applies them to the rhs

  double c1_;    ///< first Mooney-Rivlin parameter
  double c2_;   ///< second Mooney-Rivlin parameter
  double displacementsScalingFactor_;   ///< factor with which to scale the displacements
  bool dumpDenseMatlabVariables_;      ///< the current vector x, the residual, r and the jacobian, jac should be written

  bool useNestedMat_ = false;   ///< if the MatNest and VecNest data structures of Petsc should be used, this avoids data copy but is harder to debug
  bool useAnalyticJacobian_;   ///< if the analytically computed Jacobian of the Newton scheme should be used. Theoretically if it is correct, this is the fastest option.
  bool useNumericJacobian_;   ///< if a numerically computed Jacobian should be used, approximated by finite differences
};

}  // namespace
