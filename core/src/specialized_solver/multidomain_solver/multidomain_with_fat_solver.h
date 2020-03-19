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
  void setSystemMatrix(double timeStepWidth);

  //! solve the linear system of equations of the implicit scheme with rightHandSide_ and solution_
  void solveLinearSystem();

  DataFat dataFat_;  //< the data object of the multidomain solver with fat, which stores all field variables and matrices

  FiniteElementMethodDiffusionFat finiteElementMethodFat_;   //< the finite element object that is used for the Laplace problem of the potential flow, needed for the fiber directions
};

}  // namespace

#include "specialized_solver/multidomain_solver/multidomain_with_fat_solver.tpp"
