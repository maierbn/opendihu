#pragma once

#include <vector>

#include "spatial_discretization/finite_element_method/00_base.h"
#include "spatial_discretization/finite_element_method/solid_mechanics/solid_mechanics_utility.h"
#include "spatial_discretization/finite_element_method/solid_mechanics/solid_mechanics_boundary_conditions.h"
#include "spatial_discretization/finite_element_method/01_matrix.h"
#include "equation/mooney_rivlin_incompressible.h"
#include "equation/type_traits.h"
#include "mesh/face_t.h"

namespace SpatialDiscretization
{

/** This class performs the nonlinear (Newton) solution process.
 * The nonlinear function to be solved is given by object->computeInternalMinusExternalVirtualWork(f);
 */
template<typename FunctionSpaceType, typename QuadratureType, typename Term>
class SolidMechanicsNonlinearSolve :
  public FiniteElementMethodBase<FunctionSpaceType, QuadratureType, Term>,
  public SolidMechanicsBoundaryConditions<FunctionSpaceType, Term>
{
public:
  // use constructor of base class
  using FiniteElementMethodBase<FunctionSpaceType, QuadratureType, Term>::FiniteElementMethodBase;

protected:
  //! solve nonlinear system
  virtual void solve() override;

  //! update the current configuration geometry from reference geometry and displacements
  void updateGeometryActual();

  //! debugging method
  void debug();

  //--- methods that are used by this class but are defined further down in the class hierarchy ----
  //! set the internal displacement variable as copy of the given values
  virtual void setFromSolverVariableSolution(Vec &solverSolutionVariable) = 0;

  //! compute and return the appropriate analytical stiffness matrix
  virtual void computeAnalyticStiffnessMatrix(std::shared_ptr<PartitionedPetscMat<FunctionSpaceType>> solverStiffnessMatrix) = 0;

  //! get the number of unknowns, also counting displacement values for which Dirichlet BC are set as unknown, for mixed formulation sum of u and p unknowns
  virtual int nLocalUnknowns() = 0;

  //! write output of current state
  virtual void writeOutput() = 0;

  //! compute the nonlinear function
  virtual void evaluateNonlinearFunction(Vec &result) = 0;

};

};  // namespace

#include "spatial_discretization/finite_element_method/solid_mechanics/solid_mechanics_nonlinear_solve.tpp"
