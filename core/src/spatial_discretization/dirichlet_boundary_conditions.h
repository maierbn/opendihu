#pragma once

#include <Python.h>  // has to be the first included header
#include <vector>
#include <memory>

namespace SpatialDiscretization
{

/** A class that handles Dirichlet type boundary conditions for scalar fields.
 *  The static method parseBoundaryConditions parses BC for field variables with any number of components.
  */
template<typename FunctionSpaceType, int nComponents>
class DirichletBoundaryConditionsBase
{
public:
  typedef std::array<double,nComponents> ValueType;   ///< the type of value of one dirichlet boundary condition

  struct ElementWithNodes
  {
    element_no_t elementNoLocal;   ///< local element no
    std::vector<std::pair<int,ValueType>> elementalDofIndex;   ///< the element-local dof index and the value of the boundary condition on this dof
  };

  //! initialize data structures by parsing boundary conditions from python config, key "dirichletBoundaryConditions"
  void initialize(PyObject *specificSettings, std::shared_ptr<FunctionSpaceType> functionSpace);

  //! initialize directly from given vectors
  void initialize(std::shared_ptr<FunctionSpaceType> functionSpace, std::vector<ElementWithNodes> &boundaryConditionElements_,
                  std::vector<dof_no_t> &boundaryConditionNonGhostDofLocalNos_, std::vector<ValueType> &boundaryConditionValues_);

  //! set the boundary condition dofs to the values
  void applyInVector(std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,nComponents>> fieldVariable);


protected:

  //! fill auxiliary ghost element data structures
  virtual void initializeGhostElements(){};

  //! parse config and extract boundary conditions specified under key "dirichletBoundaryConditions", store in boundaryConditions
  void parseBoundaryConditions(PyObject *settings, std::shared_ptr<FunctionSpaceType> functionSpace,
                               std::vector<std::pair<int,std::array<double,nComponents>>> &boundaryConditions);

  //! parse boundary conditions from config and store them in boundaryConditionElements_, boundaryConditionNonGhostDofLocalNos_ and boundaryConditionValues_
  void parseBoundaryConditionsForElements();

  //! print the data of Dirichlet Boundary Conditions to VLOG(1)
  void printDebuggingInfo();

  PyObject *specificSettings_;            ///< the python config that contains the boundary conditions
  std::shared_ptr<FunctionSpaceType> functionSpace_;     ///< function space for which boundary conditions are specified

  std::vector<ElementWithNodes> boundaryConditionElements_;   ///< nodes grouped by elements on which boundary conditions are specified
  std::vector<dof_no_t> boundaryConditionNonGhostDofLocalNos_;        ///< vector of all local (non-ghost) boundary condition dofs
  std::vector<ValueType> boundaryConditionValues_;               ///< vector of the local prescribed values, related to boundaryConditionNonGhostDofLocalNos_
};

template<typename FunctionSpaceType, int nComponents>
class DirichletBoundaryConditions :
  public DirichletBoundaryConditionsBase<FunctionSpaceType,nComponents>
{
};

/* Specialization for variables with 1 component as used in normal finite element computations.
 */
template<typename FunctionSpaceType>
class DirichletBoundaryConditions<FunctionSpaceType,1> :
  public DirichletBoundaryConditionsBase<FunctionSpaceType,1>
{
public:

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

#include "spatial_discretization/dirichlet_boundary_conditions.tpp"
