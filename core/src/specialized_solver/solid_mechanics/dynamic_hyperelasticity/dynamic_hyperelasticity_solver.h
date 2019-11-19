#pragma once

#include <Python.h>  // has to be the first included header

#include "time_stepping_scheme/time_stepping_scheme.h"
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

  typedef SpatialDiscretization::HyperelasticitySolver<Term> StaticSolverType;
  typedef typename StaticSolverType::DisplacementsFunctionSpace DisplacementsFunctionSpace;
  typedef typename StaticSolverType::PressureFunctionSpace PressureFunctionSpace;

  typedef PartitionedPetscVecForHyperelasticity<DisplacementsFunctionSpace,PressureFunctionSpace> VecHyperelasticity;
  typedef PartitionedPetscMatForHyperelasticity<DisplacementsFunctionSpace,PressureFunctionSpace> MatHyperelasticity;

  //! constructor
  DynamicHyperelasticitySolver(DihuContext context);

  //! advance simulation by the given time span
  void advanceTimeSpan();

  //! initialize everything for the simulation
  void initialize();

  //! run the whole simulation, repeatedly calls advanceTimeSpan
  void run();

private:

  //! compute the mass matrix and the inverse lumped mass matrix
  void initializeMassMatrix();

  //! compute the acceleration variable as a = M^-1(f - Ku - Cv), u and v are input, a is output
  void computeAcceleration(std::shared_ptr<VecHyperelasticity> u, std::shared_ptr<VecHyperelasticity> v, std::shared_ptr<VecHyperelasticity> a);

  //! compute damping d_LaMb = ∫φ_Lb*μ*φ_Mb*v_Mb dV (no sum over b, but over M)
  void addDamping(std::shared_ptr<VecHyperelasticity> v, std::shared_ptr<VecHyperelasticity> damping);

  //! compute the next displacements and velocities by the RK-4 scheme
  void computeRungeKutta4();

  StaticSolverType staticSolver_;  //< hyperelasticity solver that solver the static problem
  Data::DynamicHyperelasticitySolver<DisplacementsFunctionSpace> data_;

  double density_;   //< density rho, used for inertia
  double viscosity_;  //< viscosity mu, used for damping, set to 0 to disable damping

  std::shared_ptr<MatHyperelasticity> massMatrix_;    //< mass matrix
  std::shared_ptr<MatHyperelasticity> inverseLumpedMassMatrix_;    //< mass matrix with inverse row sums on diagonal as combined matrix for u and p, only the u part contains the lumped mass matrix

  std::shared_ptr<VecHyperelasticity> u_;     //< helper variable for displacements u, combined vector of u and p values
  std::shared_ptr<VecHyperelasticity> v_;     //< helper variable for velocities v, combined vector of u and p values
  std::shared_ptr<VecHyperelasticity> a_;     //< helper variable for acceleration a, combined vector of u and p values
  std::array<std::shared_ptr<VecHyperelasticity>,4> k_;   //< intermediate values for the RK4 scheme
  std::array<std::shared_ptr<VecHyperelasticity>,4> l_;   //< intermediate values for the RK4 scheme
  std::array<std::shared_ptr<VecHyperelasticity>,9> temp_;  //< temporary values for the RK4 scheme
};

}  // namespace

#include "specialized_solver/solid_mechanics/dynamic_hyperelasticity/dynamic_hyperelasticity_solver.tpp"
