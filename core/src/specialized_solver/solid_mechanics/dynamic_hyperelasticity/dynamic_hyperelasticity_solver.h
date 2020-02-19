#pragma once

#include <Python.h>  // has to be the first included header

#include "time_stepping_scheme/00_time_stepping_scheme.h"
#include "specialized_solver/solid_mechanics/hyperelasticity/hyperelasticity_solver.h"
#include "data_management/specialized_solver/dynamic_hyperelasticity_solver.h"

namespace TimeSteppingScheme
{

/** This is a solver for a dynamic nonlinear finite elasticity problem.
 *  It uses the hyperelasticity solver for the static computations.
 *  The dynamic problem includes the static problem, only the right hand side changes to account for inertia effects.
 *  Furthermore, an integration scheme for 2nd order ODEs is needed. This is given by leap frog integration.
 *
 */
template<typename Term = Equation::SolidMechanics::MooneyRivlinIncompressible3D>
class DynamicHyperelasticitySolver :
  public TimeSteppingScheme
{
public:

  typedef SpatialDiscretization::HyperelasticitySolver<Term,6> HyperelasticitySolverType;    // the hyperelasticity solver that solves the nonlinear problem, 6 non-pressure components (u and v)
  typedef typename HyperelasticitySolverType::DisplacementsFunctionSpace DisplacementsFunctionSpace;
  typedef typename HyperelasticitySolverType::PressureFunctionSpace PressureFunctionSpace;

  typedef PartitionedPetscVecForHyperelasticity<DisplacementsFunctionSpace,PressureFunctionSpace,6> VecHyperelasticity;
  typedef PartitionedPetscMatForHyperelasticity<DisplacementsFunctionSpace,PressureFunctionSpace,6> MatHyperelasticity;

  typedef Data::DynamicHyperelasticitySolver<DisplacementsFunctionSpace> Data;

  //! constructor
  DynamicHyperelasticitySolver(DihuContext context);

  //! advance simulation by the given time span
  void advanceTimeSpan();

  //! initialize everything for the simulation
  void initialize();

  //! run the whole simulation, repeatedly calls advanceTimeSpan
  void run();

  //! return the data object, with the call to this method the output writers get the data to create their output files
  Data &data();

  //! get a reference to the underlying HyperelasticitySolver which has the material formulation and the nonlinear solver
  HyperelasticitySolverType &hyperelasticitySolver();


private:

  //! set initial values for u and v from settings
  void setInitialValues();

  //! call the callback function to update dirichlet boundary condition values
  void callUpdateDirichletBoundaryConditionsFunction(double t);

  HyperelasticitySolverType hyperelasticitySolver_;  //< hyperelasticity solver that solver the static problem
  Data data_;

  double density_;   //< density rho, used for inertia

  std::shared_ptr<VecHyperelasticity> uvp_;     //< combined vector of u,v and p values
  Vec internalVirtualWork_;                     //< internal virtual work, computed by the hyperelasticity solver
  Vec accelerationTerm_;                        //< contribution to virtual work from acceleration, computed by the hyperelasticity solver
  Vec externalVirtualWorkDead_;                 //< external virtual work, computed by the hyperelasticity solver, dead load i.e. constant over time

  bool inputMeshIsGlobal_;                      //< value of the setting "inputMeshIsGlobal", if the new dirichletBC values are given in global or local numbering

  PyObject *pythonUpdateDirichletBoundaryConditionsFunction_;       //< the callback function
  int updateDirichletBoundaryConditionsFunctionCallInterval_;       //< the interval with which the function will be called
  int updateDirichletBoundaryConditionsFunctionCallCount_ = 0;      //< the counter of number of call to the updateDirichletBoundaryConditionsFunction
};

}  // namespace

#include "specialized_solver/solid_mechanics/dynamic_hyperelasticity/dynamic_hyperelasticity_solver.tpp"
