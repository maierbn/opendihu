#pragma once

#include <Python.h>  // has to be the first included header
#include <vector>
#include <memory>

#include "output_writer/output_surface/output_points.h"

namespace SpatialDiscretization
{

/** A class that parses boundary conditions for scalar fields. This is the base class for both Dirichlet and Neumann BC.
 *  The static method parseBoundaryConditions parses BC for field variables with any number of components.
  */
template<typename FunctionSpaceType, int nComponents>
class DirichletBoundaryConditionsBase
{
public:
  typedef std::array<double,nComponents> ValueType;   //< the type of value of one boundary condition

  struct ElementWithNodes
  {
    element_no_t elementNoLocal;   //< local element no
    std::map<int,ValueType> elementalDofIndex;   //< the element-local dof index and the value of the boundary condition on this dof
  };

  struct BoundaryConditionsForComponent
  {
    std::vector<dof_no_t> dofNosLocal;    // the local dof no (without ghost dofs), both vectors are sorted according to this
    std::vector<double> values;           // the values corresponding to the dofs
  };

  //! constructor
  DirichletBoundaryConditionsBase(DihuContext context);

  //! initialize data structures by parsing boundary conditions from python config, keys "dirichletBoundaryConditions" or "neumannBoundaryConditions"
  virtual void initialize(PythonConfig specificSettings, std::shared_ptr<FunctionSpaceType> functionSpace,
                          std::string boundaryConditionsConfigKey);

  //! initialize directly from given vectors
  void initialize(std::shared_ptr<FunctionSpaceType> functionSpace, std::vector<ElementWithNodes> &boundaryConditionElements,
                  std::vector<dof_no_t> &boundaryConditionNonGhostDofLocalNos, std::vector<ValueType> &boundaryConditionValues);

  //! add boundary conditions to the currently present boundary conditions
  //! @param overwriteBcOnSameDof if existing bc dofs that are also in the ones to set newly should be overwritten, else they are not touched
  void addBoundaryConditions(std::vector<ElementWithNodes> &boundaryConditionElements, bool overwriteBcOnSameDof);

  //! get a reference to the vector of bc local dofs
  const std::vector<dof_no_t> &boundaryConditionNonGhostDofLocalNos() const;

  //! get a reference to the vector of bc local dofs
  const std::vector<ValueType> &boundaryConditionValues() const;

  //! get the boundary conditions data organized by component
  const std::array<BoundaryConditionsForComponent, nComponents> &boundaryConditionsByComponent() const;

protected:

  //! fill auxiliary ghost element data structures, this is only needed for Dirichlet boundary conditions on scalar fields
  virtual void initializeGhostElements(){};

  //! parse config and extract boundary conditions specified under key "dirichletBoundaryConditions", store in boundaryConditions
  void parseBoundaryConditions(PythonConfig settings, std::shared_ptr<FunctionSpaceType> functionSpace, std::string boundaryConditionsConfigKey,
                               std::vector<std::pair<int,std::array<double,nComponents>>> &boundaryConditions);

  //! parse boundary conditions from config and store them in boundaryConditionElements_, boundaryConditionNonGhostDofLocalNos_ and boundaryConditionValues_
  void parseBoundaryConditionsForElements(std::string boundaryConditionsConfigKey);

  //! print the data to VLOG(1)
  void printDebuggingInfo();

  //! create the boundaryConditionsByComponent_ data structure from boundaryConditionNonGhostDofLocalNos_ and boundaryConditionValues_
  void generateBoundaryConditionsByComponent();

  //! write the constraint dofs and the prescribed values as points in a vtk file
  void writeOutput();

  PythonConfig specificSettings_;                               //< the python config that contains the boundary conditions
  std::shared_ptr<FunctionSpaceType> functionSpace_;            //< function space for which boundary conditions are specified
  OutputWriter::OutputPoints outputPoints_;                     //< output writer that writes points of the Dirichlet BCs

  std::vector<ElementWithNodes> boundaryConditionElements_;     //< nodes grouped by elements on which boundary conditions are specified, this includes ghost nodes
  std::vector<dof_no_t> boundaryConditionNonGhostDofLocalNos_;  //< vector of all local (non-ghost) boundary condition dofs, sorted
  std::vector<ValueType> boundaryConditionValues_;              //< vector of the local (non-ghost) prescribed values, related to boundaryConditionNonGhostDofLocalNos_

  std::array<BoundaryConditionsForComponent, nComponents> boundaryConditionsByComponent_;   //< the local boundary condition data organized by component, entries are sorted by dofNoLocal, without ghost dofs

  std::string filenameOutput_;                                  //< output filename for the vtk file that will contain all Dirichlet boundary conditions for visualization
};

template<typename FunctionSpaceType, int nComponents>
std::ostream &operator<<(std::ostream &stream, const typename DirichletBoundaryConditionsBase<FunctionSpaceType,nComponents>::BoundaryConditionsForComponent &rhs);

template<typename FunctionSpaceType, int nComponents>
std::ostream &operator<<(std::ostream &stream, const typename DirichletBoundaryConditionsBase<FunctionSpaceType,nComponents>::ElementWithNodes &rhs);

} // namespace

#include "spatial_discretization/dirichlet_boundary_conditions/00_dirichlet_boundary_conditions_base.tpp"
