#pragma once

#include <Python.h>  // has to be the first included header

#include "interfaces/runnable.h"
#include "solver/nonlinear.h"
#include "specialized_solver/solid_mechanics/hyperelasticity/01_material_computations.h"

namespace SpatialDiscretization
{

/** This solver is for the nonlinear finite elasticity problem with Mooney-Rivlin material in 3D.
 *  It implements the static equation ∂W_int(u,p) - ∂W_ext = 0 + incompressibility (for nDisplacementComponents = 3)
 *  and the dynamic equation ∂W_int(u,p) - ∂W_ext,dead + ∂W_ext(v) = 0, dot u = v, incompressibility (see details in doc.pdf) (for nDisplacementComponents = 6)
 *
 * Further numerical improvements for solving the nonlinear system that may be implemented later:
 *   https://github.com/jedbrown/spectral-petsc/blob/master/stokes.C
 *   https://lists.mcs.anl.gov/pipermail/petsc-dev/2008-April/000711.html
 *
 * The implementation of this solver is contained in the following files:
 *
 * - "specialized_solver/solid_mechanics/hyperelasticity/02_hyperelasticity_solver.tpp":
 * This contains the top-level methods and the initialization of data structures.
 *
 * - "specialized_solver/solid_mechanics/hyperelasticity/material_computations.tpp":
 * This contains the implementation of equations, like the residual Wint-Wext, the stress, the analytic jacobian
 *
 * - "specialized_solver/solid_mechanics/hyperelasticity/material_testing.tpp"
 * This contains methods to test the implementations in material_computations.tpp. It is not used for production.
 *
 * - "specialized_solver/solid_mechanics/hyperelasticity/02_petsc_callbacks.tpp"
 * This contains plain callback functions that are passed to the PETSc nonlinear solver.
 * They call the HyperelasticitySolver object through methods that are implemented in nonlinear_solve.tpp.
 *
 * - "specialized_solver/solid_mechanics/hyperelasticity/02_nonlinear_solve.tpp"
 * This contains interfaces that are called by PETSc during the nonlinear solution process.
 *
 * - The material equation is given by the structs in equation/mooney_rivlin_incompressible.h, e.g.
 *    Equation::SolidMechanics::MooneyRivlinIncompressible3D or Equation::SolidMechanics::TransverselyIsotropicMooneyRivlinIncompressibleActive3D
 *
 * The template parameter @param nDisplacementComponents refers to the number of non-pressure components
 * and decides if the static (3) or the dynamic (6) case should be computed.
  */
template<typename Term = Equation::SolidMechanics::MooneyRivlinIncompressible3D, bool withLargeOutput=true, typename MeshType = Mesh::StructuredDeformableOfDimension<3>, int nDisplacementComponents = 3>
class HyperelasticitySolver :
  public HyperelasticityMaterialComputations<Term,withLargeOutput,MeshType,nDisplacementComponents>,
  public Runnable
{
public:
  typedef FunctionSpace::FunctionSpace<MeshType, BasisFunction::LagrangeOfOrder<2>> DisplacementsFunctionSpace;
  typedef FunctionSpace::FunctionSpace<MeshType, BasisFunction::LagrangeOfOrder<1>> PressureFunctionSpace;
  typedef FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<1>, BasisFunction::LagrangeOfOrder<1>> FiberFunctionSpace;
  typedef DisplacementsFunctionSpace FunctionSpace;

  typedef Data::QuasiStaticHyperelasticity<PressureFunctionSpace,DisplacementsFunctionSpace,Term,withLargeOutput,Term> Data;
  typedef ::Data::QuasiStaticHyperelasticityPressureOutput<PressureFunctionSpace> PressureDataCopy;

  typedef FieldVariable::FieldVariable<DisplacementsFunctionSpace,3> DisplacementsFieldVariableType;
  typedef FieldVariable::FieldVariable<PressureFunctionSpace,1> PressureFieldVariableType;
  typedef FieldVariable::FieldVariable<DisplacementsFunctionSpace,6> StressFieldVariableType;

  typedef PartitionedPetscVecForHyperelasticity<DisplacementsFunctionSpace,PressureFunctionSpace,Term,nDisplacementComponents> VecHyperelasticity;
  typedef PartitionedPetscMatForHyperelasticity<DisplacementsFunctionSpace,PressureFunctionSpace,Term,nDisplacementComponents> MatHyperelasticity;

  //! constructor
  HyperelasticitySolver(DihuContext context, std::string settingsKey = "HyperelasticitySolver");

  //! advance simulation by the given time span, data in solution is used, afterwards new data is in solution
  void advanceTimeSpan();

  //! dummy method, set endTime as current output time
  void setTimeSpan(double startTime, double endTime);

  //! run the simulation
  void run();

  //! this evaluates the actual nonlinear function f(x) that should be solved f(x) = 0
  //! @return if computation was successful
  bool evaluateNonlinearFunction(Vec x, Vec f);

  //! this evaluates the analytic jacobian for the Newton scheme, or the material stiffness
  //! @return if computation was successful
  bool evaluateAnalyticJacobian(Vec x, Mat jac);

  //! callback after each nonlinear iteration
  void monitorSolvingIteration(SNES snes, PetscInt its, PetscReal norm);

protected:

  typedef HyperelasticityMaterialComputations<Term,withLargeOutput,MeshType,nDisplacementComponents> Parent;

  //! use Petsc to solve the nonlinear equation using the SNES solver
  void nonlinearSolve();

  //! do some steps after nonlinearSolve(): close log file, copy the solution values back to this->data_.displacements() and this->data.pressure(),
  //! compute the PK2 stress at every node, update the geometry field by the new displacements, dump files containing rhs and system matrix
  void postprocessSolution();

  //! check if the solution satisfies Dirichlet BC and the residual is zero, only for debugging output
  void checkSolution(Vec x);

  //! set all PETSc callback functions, e.g. for computation jacobian or the nonlinear function itself
  void initializePetscCallbackFunctions();

  using Parent::lastSolution_;           //< a temporary variable to hold the previous solution in the nonlinear solver, to be used to reset the nonlinear scheme if it diverged
  using Parent::bestSolution_;           //< a temporary variable to hold the best solution so, the one with the lowest residual norm

  using Parent::loadFactors_;            //< vector of load factors, 1.0 means normal computation, any lower value reduces the right hand side (scales body and traction forces)

  using Parent::lastNorm_;               //< residual norm of the last iteration in the nonlinear solver
  using Parent::secondLastNorm_;         //< residual norm of the second last iteration in the nonlinear solver
  using Parent::bestResidualNorm_;       //< best residual norm for load factor 1.0 achieved so far
  using Parent::currentLoadFactor_;      //< current value of the load factor, this value is passed to materialComputeResidual(), 1.0 means normal computation, any lower value reduces the right hand side (scales body and traction forces)
  using Parent::previousLoadFactor_;     //< previous value of the load factor
  using Parent::lastSolveSucceeded_;     //< if the last computation of the residual or jacobian succeeded, if this is false, it indicates that there was a negative jacobian
  using Parent::loadFactorGiveUpThreshold_;   //< a threshold for the load factor, if it is below, the solve is aborted

  using Parent::endTime_;                //< end time of the simulation
  using Parent::combinedVecResidual_;    //< the Vec for the residual and result of the nonlinear function
  using Parent::combinedVecSolution_;    //< the Vec for the solution, combined means that ux,uy,uz and p components are combined in one vector

  using Parent::solverVariableResidual_; //< PETSc Vec to store the residual, equal to combinedVecResidual_->valuesGlobal()
  using Parent::solverVariableSolution_; //< PETSc Vec to store the solution, equal to combinedVecSolution_->valuesGlobal()

};

}  // namespace

#include "specialized_solver/solid_mechanics/hyperelasticity/02_hyperelasticity_solver.tpp"
#include "specialized_solver/solid_mechanics/hyperelasticity/02_nonlinear_solve.tpp"
