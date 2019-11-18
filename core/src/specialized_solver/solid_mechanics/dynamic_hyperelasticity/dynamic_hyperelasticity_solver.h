#pragma once

#include <Python.h>  // has to be the first included header

#include <time_stepping_scheme/time_stepping_scheme.h>

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

  //! constructor
  DynamicHyperelasticitySolver(DihuContext context);

  //! advance simulation by the given time span
  void advanceTimeSpan();

  //! initialize everything for the simulation
  void initialize();

  //! run the whole simulation, repeatedly calls advanceTimeSpan
  void run();

private:

  HyperelasticitySolver<Term> staticSolver_;  //< hyperelasticity solver that solver the static problem

};

}  // namespace

#include "specialized_solver/solid_mechanics/dynamic_hyperelasticity/dynamic_hyperelasticity_solver.tpp"
