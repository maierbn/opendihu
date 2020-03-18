#pragma once

#include <Python.h>  // has to be the first included header
#include "time_stepping_scheme/02_time_stepping_scheme_ode.h"
#include "interfaces/runnable.h"
#include "data_management/specialized_solver/multidomain.h"
#include "control/dihu_context.h"
#include "partition/rank_subset.h"

namespace TimeSteppingScheme
{

/** A specialized solver for the multidomain equation, as formulated by Thomas Klotz (2017)
 *
 * The system to be solved here is
 *
 * [ -dt/(a_mk*c_mk)*M^{-1}*K_ik + I   ...   -dt/(a_mk*c_mk)*M^{-1}*K   ] [V_mk^(i+1)  ]   [V_mk^(i)]
 * [  ...                                     ...                       ]*[...         ] = [       ]
 * [ f_rk * K_ik                       ...    K_ei                      ] [phi_e^(i+1) ]   [       ]
 *
 * V_mk^(i) is computed by the 0D part.
 * The two output connection slots are V_mk^(i) and V_mk^(i+1),
 * where V_mk^(i) is the input and should be connected to the output of the reaction term.
 * V_mk^(i+1) is the output and should be connected to the input of the reaction term.
 */
template<typename FiniteElementMethodPotentialFlow,typename FiniteElementMethodDiffusion>
class MultidomainSolver :
  public TimeSteppingScheme, public Runnable
{
public:
  typedef typename FiniteElementMethodDiffusion::FunctionSpace FunctionSpace;
  typedef typename Data::Multidomain<typename FiniteElementMethodDiffusion::FunctionSpace>::FieldVariableType FieldVariableType;
  typedef typename Data::Multidomain<typename FiniteElementMethodDiffusion::FunctionSpace> Data;
  typedef typename Data::OutputConnectorDataType OutputConnectorDataType;

  //! constructor
  MultidomainSolver(DihuContext context);

  //! advance simulation by the given time span, data in solution is used, afterwards new data is in solution
  void advanceTimeSpan();

  //! initialize components of the simulation
  void initialize();

  //! run the simulation
  void run();

  //! return the data object
  Data &data();

  //! get the data that will be transferred in the operator splitting to the other term of the splitting
  //! the transfer is done by the output_connector_data_transfer class
  std::shared_ptr<OutputConnectorDataType> getOutputConnectorData();

protected:

  //! assemble the system matrix which is a block matrix containing stiffness matrices of the diffusion sub problems
  void setSystemMatrix(double timeStepWidth);

  //! solve the linear system of equations of the implicit scheme with rightHandSide_ and solution_
  void solveLinearSystem();

  //! initialize the relative factors fr_k
  void initializeCompartmentRelativeFactors();

  Data dataMultidomain_;  //< the data object of the multidomain solver which stores all field variables and matrices

  FiniteElementMethodPotentialFlow finiteElementMethodPotentialFlow_;   //< the finite element object that is used for the Laplace problem of the potential flow, needed for the fiber directions
  std::vector<FiniteElementMethodDiffusion> finiteElementMethodDiffusionCompartment_;   //< the finite element object that is used for the diffusion of the compartments, with prefactor f_r
  FiniteElementMethodDiffusion finiteElementMethodDiffusion_;   //< the finite element object that is used for the diffusion with diffusion tensor sigma_i, without prefactor
  FiniteElementMethodDiffusion finiteElementMethodDiffusionTotal_;   //< the finite element object that is used for the diffusion with diffusion tensor (sigma_i + sigma_e), bottom right block of system matrix

  std::shared_ptr<Solver::Linear> linearSolver_;   //< the linear solver used for solving the system
  std::shared_ptr<Partition::RankSubset> rankSubset_;  //< the rankSubset for all involved ranks

  int nCompartments_;                         //< the number of instances of the diffusion problem, or the number of motor units
  Mat nestedSystemMatrix_;                    //< nested Petsc Mat, the system matrix which has more components than dofs, later this should be placed inside the data object
  Vec nestedSolution_;                        //< nested Petsc Vec, solution vector
  Vec nestedRightHandSide_;                   //< nested Petsc Vec, rhs
  Mat singleSystemMatrix_;                    //< non-nested Petsc Mat that contains all entries, the system matrix
  Vec singleSolution_;                        //< non-nested Petsc Vec, solution vector
  Vec singleRightHandSide_;                   //< non-nested Petsc Vec, distributed rhs

  std::vector<Vec> subvectorsRightHandSide_;  //< the sub vectors that are used in the nested vector nestedRightHandSide_
  std::vector<Vec> subvectorsSolution_;       //< the sub vectors that are used in the nested vector nestedSolution_

  std::vector<double> am_, cm_;  //< the Am and Cm prefactors for the compartments, Am = surface-volume ratio, Cm = capacitance
  bool initialGuessNonzero_;     //< if the initial guess should be set to the last solution
  int lastNumberOfIterations_;   //< the number of iterations that were needed the last time to solve the linear system
};

}  // namespace

#include "specialized_solver/multidomain_solver/multidomain_solver.tpp"
