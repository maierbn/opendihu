#pragma once

#include "spatial_discretization/finite_element_method/04_rhs.h"

#include "mesh/mesh.h"
#include "discretizable_in_time/discretizable_in_time.h"

namespace SpatialDiscretization
{

/** class used for timestepping as for diffusion equation
 */
template<typename FunctionSpaceType, typename QuadratureType, typename Term>
class FiniteElementMethodTimeStepping :
  public AssembleRightHandSide<FunctionSpaceType, QuadratureType, Term>,
  public DiscretizableInTime
{
public:
  FiniteElementMethodTimeStepping(DihuContext context);

  //! return the compile-time constant number of variable components of the solution field variable
  static constexpr int nComponents();

  //! proceed time stepping by computing output = stiffnessMatrix*input, output back in strong form
  void evaluateTimesteppingRightHandSideExplicit(Vec &input, Vec &output, int timeStepNo, double currentTime);
  
  //! timestepping rhs of equation Au^(t+1)=Rhs^(t), used for the case (M/dt-K)u^(t+1)=M/dtu^(t)
  //void evaluateTimesteppingRightHandSideImplicit(Vec &input, Vec &output, int timeStepNo, double currentTime);

  //! initialize for use with timestepping
  void initialize();
  
  //! initialize for use with timestepping
  void initialize(double timeStepWidth); 

  //! reset the object to uninitialized state
  void reset();

  //! hook to set initial values for a time stepping from this FiniteElement context, return true if it has set the values or don't do anything and return false
  bool setInitialValues(std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,1>> initialValues);

  //! set the subset of ranks that will compute the work
  void setRankSubset(Partition::RankSubset rankSubset);

  //! return true because the object has a specified mesh type
  bool knowsMeshType();

  //! enable or disable boundary condition handling on initialization, set to false to not care for boundary conditions
  void setBoundaryConditionHandlingEnabled(bool boundaryConditionHandlingEnabled);

  //! return the mesh that is stored in the data class
  std::shared_ptr<Mesh::Mesh> mesh();

  typedef FunctionSpaceType FunctionSpace;   ///< the FunctionSpace type needed for time stepping scheme

  friend class StiffnessMatrixTester;    ///< a class used for testing
protected:

  //! do nothing, needed for initialize of base class that is overridden anyway
  void setRightHandSide(){};

  //! Compute from the rhs in weak formulation the rhs vector in strong formulation
  void computeInverseMassMatrixTimesRightHandSide(Vec &result);

  //! compute the inverse of the lumped mass matrix
  void setInverseLumpedMassMatrix();

  //! compute the system matrix for implicit timestepping, A=I-dt*M^(-1)K from the inverse of the mass matrix M^(-1) and stiffness matrix K
  //void setSystemMatrix(double timeStepWidth);
/*
  //! precomputes the system matrix A=M/dt-K for variant 1 of the implicit Euler
  void precomputeSystemMatrix1();
*/
  //! initialize the linear solve that is needed for the solution of the implicit timestepping system
  void initializeLinearSolver();

  //! solves the linear system of equations resulting from the Implicit Euler method time discretization
  void solveLinearSystem(Vec &input, Vec &output);
  
  std::shared_ptr<Solver::Linear> linearSolver_;   ///< the linear solver used for inverting the mass matrix
  std::shared_ptr<KSP> ksp_;                       ///< the linear solver context
  
};

};  // namespace

#include "spatial_discretization/finite_element_method/05_time_stepping.tpp"
#include "spatial_discretization/finite_element_method/05_time_stepping_explicit.tpp"
#include "spatial_discretization/finite_element_method/05_time_stepping_implicit.tpp"
