#pragma once

#include <Python.h>  // has to be the first included header
#include "time_stepping_scheme/time_stepping_scheme_ode.h"
#include "control/runnable.h"
#include "data_management/time_stepping.h"
#include "control/dihu_context.h"
#include "model_order_reduction/pod.h"

namespace TimeSteppingScheme
{

/** A specialized solver for the multidomain equation, as formulated by Thomas Klotz (2017)
  */
template<typename FiniteElementMethodPotentialFlow,typename CellMLAdapter,typename FiniteElementMethodDiffusion>
class MultidomainSolver :
  public TimeSteppingImplicit<FiniteElementMethodDiffusion>, public Runnable
{
public:

  //! constructor
  MultidomainSolver(DihuContext context);

  //! advance simulation by the given time span, data in solution is used, afterwards new data is in solution
  void advanceTimeSpan();

  //! initialize components of the simulation
  void initialize();

  //! run the simulation
  void run();

protected:

  //! assemble the system matrix which is a block matrix containing stiffness matrices of the diffusion sub problems
  void setSystemMatrix(double timeStepWidth);

  //! solve the linear system of equations of the implicit scheme with rightHandSide_ and solution_
  void solveLinearSystem();

  Data::Multidomain<typename DiscretizableInTimeType::FunctionSpace> dataMultidomain_;  ///< the data object of the multidomain solver which stores all field variables and matrices
  FiniteElementMethodPotentialFlow finiteElementMethodPotentialFlow_;   ///< the finite element object that is used for the Laplace problem of the potential flow, needed for the fibre directions
  FiniteElementMethodPotentialFlow finiteElementMethodDiffusion_;   ///< the finite element object that is used for the diffusion, only the stiffness matrix is computed by this object
  FiniteElementMethodPotentialFlow finiteElementMethodDiffusionTotal_;   ///< the finite element object that is used for the diffusion with diffusion tensor (sigma_i + sigma_e), bottom right block of system matrix
  std::vector<CellMLAdapter> cellMLAdapters_;   ///< the cellml adapter objects that solves the cellml rhs, e.g. Hodgkin-Huxley model

  int nCompartments_;   ///< the number of instances of the diffusion problem, or the number of motor units
  Mat systemMatrix_;    ///< for now, the system matrix which has more components than dofs, later this should be placed inside the data object
  Vec solution_;        ///< nested solution vector
  Vec rightHandSide_;             ///< distributed rhs

  std::vector<double> am_, cm_;  ///< the Am and Cm prefactors for the compartments, Am = surface-volume ratio, Cm = capacitance
};

}  // namespace

#include "time_stepping_scheme/multidomain_solver.tpp"
