#pragma once

#include <Python.h>  // has to be the first included header
#include "time_stepping_scheme/time_stepping_scheme_ode.h"
#include "interfaces/runnable.h"
#include "data_management/specialized_solver/static_bidomain.h"
#include "control/dihu_context.h"
#include "partition/rank_subset.h"

namespace TimeSteppingScheme
{

/** A specialized solver for the bidomain equation, div((sigma_i+sigma_e)*grad(phi_e)) + div(sigma_i*grad(Vm)) = 0
  */
template<typename FiniteElementMethodPotentialFlow,typename FiniteElementMethodDiffusion>
class StaticBidomainSolver :
  public Runnable
{
public:
  typedef typename FiniteElementMethodDiffusion::FunctionSpace FunctionSpace;
  typedef typename Data::StaticBidomain<typename FiniteElementMethodDiffusion::FunctionSpace>::FieldVariableType FieldVariableType;
  typedef std::shared_ptr<FieldVariableType> TransferableSolutionDataType;
  typedef typename Data::StaticBidomain<typename FiniteElementMethodDiffusion::FunctionSpace> Data;

  //! constructor
  StaticBidomainSolver(DihuContext context);

  //! advance simulation by the given time span, data in solution is used, afterwards new data is in solution
  void advanceTimeSpan();

  //! initialize components of the simulation
  void initialize();

  //! dummy method
  void setTimeSpan(double startTime, double endTime){}

  //! run the simulation
  void run();

  //! reset state
  void reset();

  //! return whether the underlying discretizableInTime object has a specified mesh type and is not independent of the mesh type
  bool knowsMeshType();

  //! return the data object
  Data &data();

  //! get the data that will be transferred in the operator splitting to the other term of the splitting
  //! the transfer is done by the solution_vector_mapping class
  TransferableSolutionDataType getSolutionForTransfer();

protected:

  //! assemble the system matrix which is a block matrix containing stiffness matrices of the diffusion sub problems
  void setSystemMatrix();

  //! solve the linear system of equations of the implicit scheme with rightHandSide_ and solution_
  void solveLinearSystem();

  DihuContext context_;    ///< object that contains the python config for the current context and the global singletons meshManager and solverManager
  Data data_;  ///< the data object of the multidomain solver which stores all field variables and matrices

  OutputWriter::Manager outputWriterManager_; ///< manager object holding all output writer

  FiniteElementMethodPotentialFlow finiteElementMethodPotentialFlow_;   ///< the finite element object that is used for the Laplace problem of the potential flow, needed for the fiber directions
  FiniteElementMethodDiffusion finiteElementMethodDiffusionTransmembrane_;   ///< the finite element object that is used for the diffusion with diffusion tensor sigma
  FiniteElementMethodDiffusion finiteElementMethodDiffusionExtracellular_;   ///< the finite element object that is used for the diffusion with diffusion tensor (sigma_i + sigma_e), bottom right block of system matrix

  std::shared_ptr<Solver::Linear> linearSolver_;   ///< the linear solver used for solving the system
  std::shared_ptr<Partition::RankSubset> rankSubset_;  ///< the rankSubset for all involved ranks

  Mat systemMatrix_;    ///< for now, the system matrix which has more components than dofs, later this should be placed inside the data object
  Vec solution_;        ///< nested solution vector
  Vec rightHandSide_;             ///< distributed rhs
  std::vector<Vec> subvectorsRightHandSide_; ///< the sub vectors that are used in the nested vector rightHandSide_
  std::vector<Vec> subvectorsSolution_; ///< the sub vectors that are used for the solution nested vector

  int lastNumberOfIterations_;   ///< the number of iterations that were needed the last time to solve the linear system
  std::string durationLogKey_;   ///< key with with the duration of the computation is written to the performance measurement log

  bool initialized_;   ///< if this object was already initialized
  PythonConfig specificSettings_;    ///< python object containing the value of the python config dict with corresponding key
};

}  // namespace

#include "specialized_solver/static_bidomain_solver.tpp"
