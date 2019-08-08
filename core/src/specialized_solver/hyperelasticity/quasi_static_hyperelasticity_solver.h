#pragma once

#include <Python.h>  // has to be the first included header

#include "interfaces/runnable.h"
#include "function_space/function_space.h"
#include "basis_function/basis_function.h"
#include "data_management/specialized_solver/quasi_static_hyperelasticity.h"
#include "output_writer/manager.h"

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

protected:

  //! use Petsc to solve the nonlinear equation
  void nonlinearSolve();

  //! set the solution variable to zero and the initial values
  void initializeSolutionVariable();

  //! zero rows and columns in jac for which dirichlet values are set, set diagonal to 1
  void applyDirichletBoundaryConditionsInJacobian(Vec x, Mat jac);

  //! set f to x-x0 for values with dirichlet boundary condition
  void applyDirichletBoundaryConditionsInNonlinearFunction(Vec x, Vec f);

  //! this evaluates the actual nonlinear function f(x) that should be solved f(x) = 0
  void evaluateNonlinearFunction(Vec x, Vec f);

  //! get a string representation of the values for debugging output
  std::string getString(Vec x);

  DihuContext context_;    ///< object that contains the python config for the current context and the global singletons meshManager and solverManager

  OutputWriter::Manager outputWriterManager_; ///< manager object holding all output writer
  Data data_;                 ///< data object

  std::string durationLogKey_;   ///< key with with the duration of the computation is written to the performance measurement log
  std::shared_ptr<DisplacementsFunctionSpace> displacementsFunctionSpace_;  ///< the function space with quadratic Lagrange basis functions, used for discretization of displacements
  std::shared_ptr<PressureFunctionSpace> pressureFunctionSpace_;  ///< the function space with linear Lagrange basis functions, used for discretization of pressure

  Mat solverMatrixTangentStiffness_;   //< the jacobian matrix for the Newton solver, which in case of nonlinear elasticity is the tangent stiffness matrix
  Vec solverVariableResidual_;         //< nested PETSc Vec to store the residual
  Vec solverVariableSolution_;         //< nested PETSc Vec to store the solution

  std::array<Mat,16> submatrices_;  // all submatrices of the 4x4 block jacobian matrix, solverMatrixTangentStiffness_
  std::array<Vec,4> subvectorsSolution_;  // all subvectors of the 4 entries vector, solverVariableSolution_
  std::array<Vec,4> subvectorsResidual_;  // all subvectors of the 4 entries vector, solverVariableResidual_
  std::vector<PartitionedPetscMat<DisplacementsFunctionSpace,DisplacementsFunctionSpace>> uMatrix_;    //< upper left 3x3 blocks of matrices in the jacobian for Newton scheme
  std::vector<PartitionedPetscMat<DisplacementsFunctionSpace,PressureFunctionSpace>> upMatrix_;        //< lower left 1x3 blocks of matrices in the jacobian for Newton scheme
  std::vector<PartitionedPetscMat<PressureFunctionSpace,DisplacementsFunctionSpace>> puMatrix_;        //< upper right 3x1 blocks of matrices in the jacobian for Newton scheme
  std::vector<PartitionedPetscMat<PressureFunctionSpace,PressureFunctionSpace>> pMatrix_;                           //< lower left matrix in the jacobian for Newton scheme

  std::shared_ptr<PartitionedPetscMat<FunctionSpace::Generic>> combinedJacobianMatrix_;    // single jacobian matrix, when useNestedMat_ is false
  std::shared_ptr<PartitionedPetscVec<FunctionSpace::Generic,1>> combinedVecResidual_;
  std::shared_ptr<PartitionedPetscVec<FunctionSpace::Generic,1>> combinedVecSolution_;

  // settings variables
  bool initialized_;   ///< if this object was already initialized
  PythonConfig specificSettings_;    ///< python object containing the value of the python config dict with corresponding key
  double endTime_;     ///< end time of current time step
  double c0_;    ///< first Mooney-Rivlin parameter
  double c1_;   ///< second Mooney-Rivlin parameter

  bool useNestedMat_ = false;   ///< if the MatNest and VecNest data structures of Petsc should be used, this avoids data copy but is harder to debug
};

}  // namespace
