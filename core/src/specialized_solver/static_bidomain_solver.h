#pragma once

#include <Python.h>  // has to be the first included header
#include "time_stepping_scheme/02_time_stepping_scheme_ode.h"
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
  typedef typename Data::StaticBidomain<typename FiniteElementMethodDiffusion::FunctionSpace> Data;
  typedef typename Data::SlotConnectorDataType SlotConnectorDataType;

  //! constructor
  StaticBidomainSolver(DihuContext context);

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

  //! get the data that will be transferred in the operator splitting to the other term of the splitting
  //! the transfer is done by the slot_connector_data class
  std::shared_ptr<SlotConnectorDataType> getSlotConnectorData();

  //! output the given data for debugging
  std::string getString(std::shared_ptr<SlotConnectorDataType> data);

protected:

  //! solve the linear system of equations of the implicit scheme with rightHandSide_ and solution_
  void solveLinearSystem();

  //! dump rhs vector
  void debugDumpData();

  DihuContext context_;    //< object that contains the python config for the current context and the global singletons meshManager and solverManager
  Data data_;              //< the data object of the multidomain solver which stores all field variables and matrices

  OutputWriter::Manager outputWriterManager_;           //< manager object holding all output writer

  FiniteElementMethodPotentialFlow finiteElementMethodPotentialFlow_;        //< the finite element object that is used for the Laplace problem of the potential flow, needed for the fiber directions
  FiniteElementMethodDiffusion finiteElementMethodDiffusionTransmembrane_;   //< the finite element object that is used for the diffusion with diffusion tensor sigma
  FiniteElementMethodDiffusion finiteElementMethodDiffusionExtracellular_;   //< the finite element object that is used for the diffusion with diffusion tensor (sigma_i + sigma_e), bottom right block of system matrix

  std::shared_ptr<Solver::Linear> linearSolver_;        //< the linear solver used for solving the system
  std::shared_ptr<Partition::RankSubset> rankSubset_;   //< the rankSubset for all involved ranks

  std::string durationLogKey_;    //< key with with the duration of the computation is written to the performance measurement log

  bool initialized_;              //< if this object was already initialized
  PythonConfig specificSettings_; //< python object containing the value of the python config dict with corresponding key
  double endTime_;                //< end time of current time step
  bool initialGuessNonzero_;      //< if the initial guess for the linear solver is set to the previous solution
};

}  // namespace

#include "specialized_solver/static_bidomain_solver.tpp"
