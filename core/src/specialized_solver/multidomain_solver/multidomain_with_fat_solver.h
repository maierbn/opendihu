#pragma once

#include <Python.h>  // has to be the first included header

#include "control/dihu_context.h"
#include "data_management/specialized_solver/multidomain_with_fat.h"

namespace TimeSteppingScheme
{

/** A specialized solver for the multidomain equation with fat layer 
 * 
*/
template<typename FiniteElementMethodPotentialFlow,typename FiniteElementMethodDiffusionMuscle,typename FiniteElementMethodDiffusionFat>
class MultidomainWithFatSolver :
  public MultidomainSolver<FiniteElementMethodPotentialFlow,FiniteElementMethodDiffusionMuscle>
{
public:
  typedef ::Data::MultidomainWithFat<
    typename FiniteElementMethodDiffusionMuscle::FunctionSpace,
    typename FiniteElementMethodDiffusionFat::FunctionSpace
  > DataFat;
  typedef typename FiniteElementMethodDiffusionMuscle::FunctionSpace FunctionSpace;
  typedef DataFat Data;

  //! constructor
  MultidomainWithFatSolver(DihuContext context);

  //! initialize components of the simulation
  void initialize();

  //! call the output writer on the data object
  virtual void callOutputWriter(int timeStepNo, double currentTime, int callCountIncrement = 1);

  //! return the data object
  DataFat &data();

protected:

  //! update the system matrix after the geometry has changed, this is done in advanceTimeSpan, if the option "updateSystemMatrixEveryTimestep" is True
  virtual void updateSystemMatrix() override;

  //! assemble the system matrix which is a block matrix containing stiffness matrices of the diffusion sub problems
  void setSystemMatrixSubmatrices(double timeStepWidth);

  //! solve the linear system of equations of the implicit scheme with rightHandSide_ and solution_
  void solveLinearSystem();

  //! initialize sharedNodes_, which contains the node nos that are shared between the muscle mesh (key) and the fat mesh (value)
  void findSharedNodesBetweenMuscleAndFat();

  //! initialize the last entry of the rhs and the matrices B,C,D,E which are special with regard to the boundary between muscle and fat mesh
  void initializeBoundaryVariables();

  //! compute the entries of the matrices B,C,D and E from the stiffness matrix of the muscle mesh
  void setEntriesBoundaryMatrices(Mat originalMatrixB, Mat originalMatrixC, Mat matrixB, Mat matrixC, Mat matrixD, Mat matrixE);

  //! compute th entries of the matrices B,C,D and E from the stiffness matrix of the muscle mesh again, the matrices have to exist already, this is used after geometry update
  void updateBoundaryMatrices();

  //! transform the numbering from normal global dof nos to a numbering that skips the shared dofs between muscle and fat mesh
  PetscInt getDofNoGlobalFatWithoutSharedDofs(PetscInt dofNoGlobal);
  
  //! transform global dof no from a shared dof in the fat mesh to the corresponding shared dof in the muscle mesh
  PetscInt getDofNoGlobalMuscleFromDofNoGlobalFat(PetscInt dofNoGlobal);
  
  //! transform global dof no from a shared dof in the muscle mesh to the corresponding shared dof in the fat mesh
  PetscInt getDofNoGlobalFatFromDofNoGlobalMuscle(PetscInt dofNoGlobal);

  //! set in the Petsc preconditioner object the information about matrix blocks for block jacobi and the node positions (PCSetCoordinates)
  virtual void setInformationToPreconditioner() override;

  //! substract the average value of phi_e and phi_b, respectively, if enabled by resetToAverageZeroPhiE and resetToAverageZeroPhiB
  void resetToAverageZero();
  
  //! set the zero Dirichlet BCs for phi_e and phi_b in the system matrix, if enabled by setDirichletBoundaryConditionPhiE_ and setDirichletBoundaryConditionPhiB_
  void setDirichletBoundaryConditionsInSystemMatrix();
  
  //! set the zero Dirichlet BCs for phi_e and phi_b in the rhs, if enabled by setDirichletBoundaryConditionPhiE_ and setDirichletBoundaryConditionPhiB_
  void setDirichletBoundaryConditionsInRightHandSide();

  // ! copy the incoming data in the phiB field variable from the dataFat_ object to the nested Vec solution which only contains not-shared dofs
  void copyPhiBToSolution();

  // ! copy the results from the linear solve in solution, which contains non-shared dofs, to the phi_b field variable in dataFat_, the missing values for the shared dofs are taken from phi_e
  void copySolutionToPhiB();

  DataFat dataFat_;  //< the data object of the multidomain solver with fat, which stores all field variables and matrices
  FiniteElementMethodDiffusionFat finiteElementMethodFat_;   //< the finite element object that is used for the Laplace problem of the potential flow, needed for the fiber directions

  std::map<node_no_t,node_no_t> sharedNodes_;   //< the node nos that are shared between the muscle mesh (key) and the fat mesh (value)
  std::set<dof_no_t> boundaryDofsFat_;            //< all local dofs without ghosts with no. in the fat mesh on the boundary ΓM
  std::set<PetscInt> boundaryDofsGlobalFat_;      //< same as boundaryDofsFat_ but for all global dofs
  std::map<PetscInt,PetscInt> fatDofToMuscleDofGlobal_;   //< the dof nos that are shared between the muscle mesh (key) and the fat mesh (value)
  std::map<PetscInt,PetscInt> muscleDofToFatDofGlobal_;   //< the dof nos that are shared between the muscle mesh (key) and the fat mesh (value)
  
  PetscInt nSharedDofsLocal_;                   //< number of shared dofs between fat and muscle mesh on the local domain
  std::vector<Mat> b1_;          //< b1^k = ((θ-1)*1/(Am^k*Cm^k)*K_sigmai^k - 1/dt*M), first factor matrix for rhs entry b, for compartment k, total: b = b1_ * Vm^(i) + b2_ * phi_e^(i)
  std::vector<Mat> b2_;          //< b2^k = (θ-1)*K_sigmai^k, second factor matrix for rhs entry b, for compartment k, total: b = b1_ * Vm^(i) + b2_ * phi_e^(i) 
  Vec temporary_;                //< temporary vector that can be multiplied right by b1_

  double theta_;                 //< θ value for Crank-Nicolson scheme
  bool useLumpedMassMatrix_;     //< if the formulation with lumped mass matrix should be used
  bool enableFatComputation_;    //< if the computation of the fat layer should be enabled, if not the values are set to 0
  bool isFirstTimestep_;         //< if the current call to solveLinearSystem is in the first timestep
};

}  // namespace

#include "specialized_solver/multidomain_solver/multidomain_with_fat_solver.tpp"
#include "specialized_solver/multidomain_solver/multidomain_with_fat_boundary_matrices.tpp"
