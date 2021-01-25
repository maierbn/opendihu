#pragma once

#include <Python.h>  // has to be the first included header

#include "function_space/function_space.h"
#include "basis_function/basis_function.h"
#include "partition/partitioned_petsc_vec/02_partitioned_petsc_vec_for_hyperelasticity.h"
#include "partition/partitioned_petsc_mat/partitioned_petsc_mat_for_hyperelasticity.h"
#include "output_writer/manager.h"
#include "data_management/specialized_solver/hyperelasticity_solver.h"
#include "spatial_discretization/dirichlet_boundary_conditions/01_dirichlet_boundary_conditions.h"
#include "spatial_discretization/neumann_boundary_conditions/01_neumann_boundary_conditions.h"
#include "specialized_solver/solid_mechanics/hyperelasticity/pressure_function_space_creator.h"

namespace SpatialDiscretization
{

/** Initialization of the HyperelasticitySolver
  */
template<typename Term = Equation::SolidMechanics::MooneyRivlinIncompressible3D, bool withLargeOutput=true, typename MeshType = Mesh::StructuredDeformableOfDimension<3>, int nDisplacementComponents = 3>
class HyperelasticityInitialize
{
public:
  typedef FunctionSpace::FunctionSpace<MeshType, BasisFunction::LagrangeOfOrder<2>> DisplacementsFunctionSpace;
  typedef FunctionSpace::FunctionSpace<MeshType, BasisFunction::LagrangeOfOrder<1>> PressureFunctionSpace;
  typedef FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<1>, BasisFunction::LagrangeOfOrder<1>> FiberFunctionSpace;
  typedef DisplacementsFunctionSpace FunctionSpace;

  typedef Data::QuasiStaticHyperelasticity<PressureFunctionSpace,DisplacementsFunctionSpace,Term,withLargeOutput,Term> Data;
  typedef ::Data::QuasiStaticHyperelasticityPressureOutput<PressureFunctionSpace> PressureDataCopy;

  //! Define the type of data that will be transferred between solvers when there is a coupling scheme.
  //! Usually you define this type in the "Data" class and reuse it here.
  typedef typename Data::SlotConnectorDataType SlotConnectorDataType;

  typedef FieldVariable::FieldVariable<DisplacementsFunctionSpace,3> DisplacementsFieldVariableType;
  typedef FieldVariable::FieldVariable<PressureFunctionSpace,1> PressureFieldVariableType;
  typedef FieldVariable::FieldVariable<DisplacementsFunctionSpace,6> StressFieldVariableType;

  typedef PartitionedPetscVecForHyperelasticity<DisplacementsFunctionSpace,PressureFunctionSpace,Term,nDisplacementComponents> VecHyperelasticity;
  typedef PartitionedPetscMatForHyperelasticity<DisplacementsFunctionSpace,PressureFunctionSpace,Term,nDisplacementComponents> MatHyperelasticity;

  //! constructor
  HyperelasticityInitialize(DihuContext context, std::string settingsKey = "HyperelasticitySolver");

  //! initialize components of the simulation
  void initialize();

  //! reset state
  void reset();

  //! return the data object
  Data &data();

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

  //! create a new shared_ptr of a PartitionedPetscVecForHyperelasticity
  std::shared_ptr<VecHyperelasticity> createPartitionedPetscVec(std::string name);

  //! create a new shared_ptr of a PartitionedPetscMatForHyperelasticity
  std::shared_ptr<MatHyperelasticity> createPartitionedPetscMat(std::string name);

  //! get the precomputed external virtual work
  Vec externalVirtualWork();

  //! get a pointer to the dirichlet boundary conditions object
  std::shared_ptr<DirichletBoundaryConditions<DisplacementsFunctionSpace,nDisplacementComponents>> dirichletBoundaryConditions();

  //! get a pointer to the neumann boundary conditions object
  std::shared_ptr<NeumannBoundaryConditions<DisplacementsFunctionSpace,Quadrature::Gauss<3>,3>> neumannBoundaryConditions();

  //! set new neumann bc's = traction for the next solve
  void updateNeumannBoundaryConditions(std::shared_ptr<NeumannBoundaryConditions<DisplacementsFunctionSpace,Quadrature::Gauss<3>,3>> newNeumannBoundaryConditions);

  //! set new dirichlet boundary condition values for existing dofs
  void updateDirichletBoundaryConditions(std::vector<std::pair<global_no_t,std::array<double,3>>> newDirichletBCValues);

  //! add new dirichlet bc's, it is also possible to set new dofs that were not prescribed beforehand
  //! this calles addBoundaryConditions() of the dirichletBoundaryConditions_ object
  //! @param overwriteBcOnSameDof if existing bc dofs that are also in the ones to set newly should be overwritten, else they are not touched
  void addDirichletBoundaryConditions(std::vector<typename DirichletBoundaryConditions<DisplacementsFunctionSpace,nDisplacementComponents>::ElementWithNodes> &boundaryConditionElements, bool overwriteBcOnSameDof);

  //! this evaluates the analytic jacobian for the Newton scheme, or the material stiffness
  //! @return if computation was successful
  virtual bool evaluateAnalyticJacobian(Vec x, Mat jac) = 0;

  //! get the Petsc Vec of the current state (uvp vector), this is needed to save and restore checkpoints from the PreciceAdapter
  Vec currentState();

  //! Get the data that will be transferred in the operator splitting or coupling to the other term of the splitting/coupling.
  //! the transfer is done by the slot_connector_data_transfer class
  std::shared_ptr<SlotConnectorDataType> getSlotConnectorData();

protected:

  //! initialize all Petsc Vec's and Mat's that will be used in the computation
  void initializePetscVariables();

  //! initialize the field variable which stores directions of fibers
  void initializeFiberDirections();

  //! set the solution variable to zero and the initial values
  void initializeSolutionVariable();

  //! set all PETSc callback functions, e.g. for computation jacobian or the nonlinear function itself
  virtual void initializePetscCallbackFunctions() = 0;

  //! determine the number of local non-zero entries in the jacobian
  unsigned int materialDetermineNumberNonzerosInJacobian();

  //! compute the elemental coordinate frame (elementalX,elementalY,elementalZ) at the node (i,j,k) in the element with node positions given by geometry
  void getElementalBasis(int i, int j, int k,
                         const std::array<Vec3,27> &geometry,
                         Vec3 &elementalX, Vec3 &elementalY, Vec3 &elementalZ);

  //! compute δW_ext,dead = int_Ω B^L * phi^L * phi^M * δu^M dx + int_∂Ω T^L * phi^L * phi^M * δu^M dS
  virtual void materialComputeExternalVirtualWorkDead() = 0;

  DihuContext context_;                                     //< object that contains the python config for the current context and the global singletons meshManager and solverManager

  OutputWriter::Manager outputWriterManager_;               //< manager object holding all output writer for displacements based variables
  OutputWriter::Manager outputWriterManagerPressure_;       //< manager object holding all output writer for pressure based variables
  OutputWriter::Manager outputWriterManagerLoadIncrements_; //< manager object holding all output writer that write during computation of several load increments
  Data data_;                                               //< data object
  PressureDataCopy pressureDataCopy_;                       //< a helper object that is used to write the pressure function space based variables with the output writers
  std::shared_ptr<Solver::Nonlinear> nonlinearSolver_;      //< the nonlinear solver object that will provide the PETSc SNES context

  std::string durationLogKey_;                              //< key with with the duration of the computation is written to the performance measurement log
  std::shared_ptr<DisplacementsFunctionSpace> displacementsFunctionSpace_;  //< the function space with quadratic Lagrange basis functions, used for discretization of displacements
  std::shared_ptr<PressureFunctionSpace> pressureFunctionSpace_;  //< the function space with linear Lagrange basis functions, used for discretization of pressure

  Mat solverMatrixJacobian_;                                //< the jacobian matrix for the Newton solver, which in case of nonlinear elasticity is the tangent stiffness matrix
  Mat solverMatrixAdditionalNumericJacobian_;               //< only used when both analytic and numeric jacobians are computed, then this holds the numeric jacobian
  Vec solverVariableResidual_;                              //< PETSc Vec to store the residual, equal to combinedVecResidual_->valuesGlobal()
  Vec solverVariableSolution_;                              //< PETSc Vec to store the solution, equal to combinedVecSolution_->valuesGlobal()
  Vec zeros_;                                               //< a solver that contains all zeros, needed to zero the diagonal of the jacobian matrix
  Vec lastSolution_;                                        //< a temporary variable to hold the previous solution in the nonlinear solver, to be used to reset the nonlinear scheme if it diverged
  Vec bestSolution_;                                        //< a temporary variable to hold the best solution so, the one with the lowest residual norm

  std::shared_ptr<VecHyperelasticity> combinedVecResidual_; //< the Vec for the residual and result of the nonlinear function
  std::shared_ptr<VecHyperelasticity> combinedVecSolution_; //< the Vec for the solution, combined means that ux,uy,uz and p components are combined in one vector
  std::shared_ptr<VecHyperelasticity> combinedVecExternalVirtualWorkDead_;      //< the Vec for the external virtual work part that does not change with u, δW_ext,dead
  std::shared_ptr<MatHyperelasticity> combinedMatrixJacobian_;                  //< single jacobian matrix
  std::shared_ptr<MatHyperelasticity> combinedMatrixAdditionalNumericJacobian_; //< only used when both analytic and numeric jacobians are computed, then this holds the numeric jacobian

  Vec externalVirtualWorkDead_;                             // the external virtual work resulting from the traction, this is a dead load, i.e. it does not change during deformation

  // settings variables
  bool initialized_;                                        //< if this object was already initialized
  PythonConfig specificSettings_;                           //< python object containing the value of the python config dict with corresponding key
  double endTime_;                                          //< end time of current time step
  std::shared_ptr<std::ofstream> residualNormLogFile_;      //< ofstream of a log file that will contain the residual norm for each iteration

  std::shared_ptr<DirichletBoundaryConditions<DisplacementsFunctionSpace,nDisplacementComponents>> dirichletBoundaryConditions_ = nullptr;  //< object that parses Dirichlet boundary conditions and applies them to rhs
  std::shared_ptr<NeumannBoundaryConditions<DisplacementsFunctionSpace,Quadrature::Gauss<3>,3>> neumannBoundaryConditions_ = nullptr;  //< object that parses Neumann boundary conditions and applies them to the rhs

  std::vector<double> materialParameters_;                  //< material parameters, e.g. c1,c2 for Mooney-Rivlin
  double displacementsScalingFactor_;                       //< factor with which to scale the displacements
  bool dumpDenseMatlabVariables_;                           //< the current vector x, the residual, r and the jacobian, jac should be written
  Vec3 constantBodyForce_;                                  //< the constant body force, if given or [0,0,0]
  double timeStepWidth_;                                    //< timeStepWidth, only need for the dynamic problem
  double density_;                                          //< density, only needed for the dynamic problem
  double lastNorm_;                                         //< residual norm of the last iteration in the nonlinear solver
  double secondLastNorm_;                                   //< residual norm of the second last iteration in the nonlinear solver
  double bestResidualNorm_;                                 //< best residual norm for load factor 1.0 achieved so far
  double currentLoadFactor_;                                //< current value of the load factor, this value is passed to materialComputeResidual(), 1.0 means normal computation, any lower value reduces the right hand side (scales body and traction forces)
  double previousLoadFactor_;                               //< previous value of the load factor
  int nNonlinearSolveCalls_;                                //< how often the nonlinear solve should be called in sequence
  bool lastSolveSucceeded_;                                 //< if the last computation of the residual or jacobian succeeded, if this is false, it indicates that there was a negative jacobian
  double loadFactorGiveUpThreshold_;                        //< a threshold for the load factor, if it is below, the solve is aborted
  unsigned int nNonZerosJacobian_;                          //< number of nonzero entries in the material jacobian on the local domain, used for preallocation of the matrix

  std::vector<double> loadFactors_;                         //< vector of load factors, 1.0 means normal computation, any lower value reduces the right hand side (scales body and traction forces)
  std::vector<double> norms_;                               //< vector that collects the norms in every iteration, it will be cleared for every new load factor

  bool useAnalyticJacobian_;                                //< if the analytically computed Jacobian of the Newton scheme should be used. Theoretically if it is correct, this is the fastest option.
  bool useNumericJacobian_;                                 //< if a numerically computed Jacobian should be used, approximated by finite differences
  bool extrapolateInitialGuess_;                            //< if the initial values for the dynamic nonlinear problem should be computed by extrapolating the previous displacements and velocities
};

}  // namespace

#include "specialized_solver/solid_mechanics/hyperelasticity/00_initialize.tpp"
