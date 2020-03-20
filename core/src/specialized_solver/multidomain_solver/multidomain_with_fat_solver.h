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

  //! constructor
  MultidomainWithFatSolver(DihuContext context);

  //! initialize components of the simulation
  void initialize();

  //! return the data object
  DataFat &data();

protected:

  //! call the output writer on the data object
  virtual void callOutputWriter(int timeStepNo, double currentTime);

  //! assemble the system matrix which is a block matrix containing stiffness matrices of the diffusion sub problems
  void setSystemMatrixSubmatrices(double timeStepWidth);

  //! solve the linear system of equations of the implicit scheme with rightHandSide_ and solution_
  void solveLinearSystem();

  //! initialize sharedNodes_, which contains the node nos that are shared between the muscle mesh (key) and the fat mesh (value)
  void findSharedNodesBetweenMuscleAndFat();

  //! initialize the last entry of the rhs and the two matrices I_ΓM and -I_ΓM that contains 1's for dof on the border between muscle and fat mesh
  void initializeBorderVariables();

  DataFat dataFat_;  //< the data object of the multidomain solver with fat, which stores all field variables and matrices
  FiniteElementMethodDiffusionFat finiteElementMethodFat_;   //< the finite element object that is used for the Laplace problem of the potential flow, needed for the fiber directions

  std::map<node_no_t,node_no_t> sharedNodes_;   //< the node nos that are shared between the muscle mesh (key) and the fat mesh (value)

  std::vector<Mat> b1_;          //< b1^k = ((θ-1)*1/(Am^k*Cm^k)*K_sigmai^k - 1/dt*M), first factor matrix for rhs entry b, for compartment k, total: b = b1_ * Vm^(i) + b2_ * phi_e^(i)
  std::vector<Mat> b2_;          //< b2^k = (θ-1)*K_sigmai^k, second factor matrix for rhs entry b, for compartment k, total: b = b1_ * Vm^(i) + b2_ * phi_e^(i) 
  Vec temporary_;                //< temporary vector that can be multiplied right by b1_

  double theta_; // θ value for Crank-Nicolson scheme
};

}  // namespace

#include "specialized_solver/multidomain_solver/multidomain_with_fat_solver.tpp"
