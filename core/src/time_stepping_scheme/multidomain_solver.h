#pragma once

#include <Python.h>  // has to be the first included header
#include "time_stepping_scheme/time_stepping_scheme_ode.h"
#include "interfaces/runnable.h"
#include "data_management/time_stepping/multidomain.h"
#include "control/dihu_context.h"
#include "partition/rank_subset.h"

namespace TimeSteppingScheme
{

/** A specialized solver for the multidomain equation, as formulated by Thomas Klotz (2017)
  */
template<typename FiniteElementMethodPotentialFlow,typename FiniteElementMethodDiffusion>
class MultidomainSolver :
  public TimeSteppingScheme, public Runnable
{
public:
  typedef typename FiniteElementMethodDiffusion::FunctionSpace FunctionSpace;
  typedef typename Data::Multidomain<typename FiniteElementMethodDiffusion::FunctionSpace>::FieldVariableType FieldVariableType;
  typedef std::pair<std::vector<Vec>,std::vector<std::shared_ptr<FieldVariableType>>> TransferableSolutionDataType;

  //! constructor
  MultidomainSolver(DihuContext context);

  //! advance simulation by the given time span, data in solution is used, afterwards new data is in solution
  void advanceTimeSpan();

  //! initialize components of the simulation
  void initialize();

  //! run the simulation
  void run();

  //! return whether the underlying discretizableInTime object has a specified mesh type and is not independent of the mesh type
  bool knowsMeshType();

  //! get the data that will be transferred in the operator splitting to the other term of the splitting
  //! the transfer is done by the solution_vector_mapping class
  TransferableSolutionDataType getSolutionForTransferInOperatorSplitting();

protected:

  //! assemble the system matrix which is a block matrix containing stiffness matrices of the diffusion sub problems
  void setSystemMatrix(double timeStepWidth);

  //! solve the linear system of equations of the implicit scheme with rightHandSide_ and solution_
  void solveLinearSystem();

  //! initialize the relative factors fr_k
  void initializeCompartmentRelativeFactors();

  Data::Multidomain<typename FiniteElementMethodDiffusion::FunctionSpace> dataMultidomain_;  ///< the data object of the multidomain solver which stores all field variables and matrices

  FiniteElementMethodPotentialFlow finiteElementMethodPotentialFlow_;   ///< the finite element object that is used for the Laplace problem of the potential flow, needed for the fiber directions
  std::vector<FiniteElementMethodDiffusion> finiteElementMethodDiffusionCompartment_;   ///< the finite element object that is used for the diffusion of the compartments, with prefactor f_r
  FiniteElementMethodDiffusion finiteElementMethodDiffusion_;   ///< the finite element object that is used for the diffusion with diffusion tensor sigma_i, without prefactor
  FiniteElementMethodDiffusion finiteElementMethodDiffusionTotal_;   ///< the finite element object that is used for the diffusion with diffusion tensor (sigma_i + sigma_e), bottom right block of system matrix

  std::shared_ptr<Solver::Linear> linearSolver_;   ///< the linear solver used for solving the system
  std::shared_ptr<Partition::RankSubset> rankSubset_;  ///< the rankSubset for all involved ranks

  int nCompartments_;   ///< the number of instances of the diffusion problem, or the number of motor units
  Mat systemMatrix_;    ///< for now, the system matrix which has more components than dofs, later this should be placed inside the data object
  Vec solution_;        ///< nested solution vector
  Vec rightHandSide_;             ///< distributed rhs
  std::vector<Vec> subvectorsRightHandSide_; ///< the sub vectors that are used in the nested vector rightHandSide_
  std::vector<Vec> subvectorsSolution_; ///< the sub vectors that are used for the solution nested vector

  std::vector<double> am_, cm_;  ///< the Am and Cm prefactors for the compartments, Am = surface-volume ratio, Cm = capacitance
  int lastNumberOfIterations_;   ///< the number of iterations that were needed the last time to solve the linear system
};

}  // namespace

#include "time_stepping_scheme/multidomain_solver.tpp"
