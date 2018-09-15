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

  Data::Multidomain<typename DiscretizableInTimeType::FunctionSpace> data_;  ///< the data object of the multidomain solver which stores all field variables and matrices
  FiniteElementMethodPotentialFlow finiteElementMethodPotentialFlow_;   ///< the finite element object that is used for the initial Laplace problem that determines the fibre direction.
  CellMLAdapter cellMLAdapter_;   ///< the cellml adapter object that solves the cellml rhs, e.g. Hodgkin-Huxley model

};

}  // namespace

#include "time_stepping_scheme/multidomain_solver.tpp"
