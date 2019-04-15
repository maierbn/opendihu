#pragma once

#include <Python.h>  // has to be the first included header
#include <vector>
#include <memory>

#include "spatial_discretization/boundary_conditions/boundary_conditions_base.h"

namespace SpatialDiscretization
{

/** A class that handles Dirichlet type boundary conditions.
  */
template<typename FunctionSpaceType, int nComponents>
class DirichletBoundaryConditionsBase :
  public BoundaryConditionsBase<FunctionSpaceType, nComponents>
{
public:

  //! use constructor of base class
  using BoundaryConditionsBase<FunctionSpaceType, nComponents>::BoundaryConditionsBase;

  //! set the Dirichlet boundary condition dofs to the values
  void applyInVector(std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,nComponents>> fieldVariable);
};

template<typename FunctionSpaceType, int nComponents>
class DirichletBoundaryConditions :
  public DirichletBoundaryConditionsBase<FunctionSpaceType,nComponents>
{
  //! use constructor of base class
  using DirichletBoundaryConditionsBase<FunctionSpaceType,nComponents>::DirichletBoundaryConditionsBase;
};

/* Specialization for variables with 1 component as used in normal finite element computations.
 */
template<typename FunctionSpaceType>
class DirichletBoundaryConditions<FunctionSpaceType,1> :
  public DirichletBoundaryConditionsBase<FunctionSpaceType,1>
{
public:
  //! use constructor of base class
  using DirichletBoundaryConditionsBase<FunctionSpaceType,1>::DirichletBoundaryConditionsBase;

  //! set the boundary conditions to system matrix, i.e. zero rows and columns of Dirichlet BC dofs and set diagonal to 1. Store the cleared matrix values in boundaryConditionsRightHandSideSummand such that they can be used for adjusting the rhs vector afterwards
  void applyInSystemMatrix(std::shared_ptr<PartitionedPetscMat<FunctionSpaceType>> systemMatrix,
                           std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,1>> boundaryConditionsRightHandSideSummand);

  //! set the boundary conditions to the right hand side using the information in boundaryConditionsRightHandSideSummand that was obtained by applyInSystemMatrix
  //! if the argument boundaryConditionsRightHandSideSummand is also set as rightHandSide, only write prescribed values into rightHandSide, do not add anything
  void applyInRightHandSide(std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,1>> rightHandSide,
                            std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,1>> boundaryConditionsRightHandSideSummand);

protected:

  //! fill auxiliary ghost element data structures
  void initializeGhostElements();

  struct GhostElement
  {
    std::vector<global_no_t> nonBoundaryConditionDofsOfRankGlobalPetsc;    ///< the non-BC dofs of this element, as global petsc no. that are owned by the rank with no neighbouringRankNo
    std::vector<global_no_t> boundaryConditionDofsGlobalPetsc;      ///< the Dirichlet BC dofs of this element
    std::vector<double> boundaryConditionValues;   ///< the prescribed value, corresponding to boundaryConditionDofsGlobalPetsc
  };
  std::map<int,std::vector<GhostElement>> foreignGhostElements_;   ///< ghost elements that are normal elements on this rank, key is the rankNo of the rank to send them to
  std::vector<GhostElement> ownGhostElements_;   ///< the ghost elements for this rank
};

} // namespace

#include "spatial_discretization/boundary_conditions/dirichlet_boundary_conditions.tpp"
