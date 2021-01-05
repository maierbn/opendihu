#pragma once

#include <Python.h>  // has to be the first included header

#include <array>
#include <memory>
#include "control/types.h"

#include "function_space/05_function_space_geometry.h"
#include "function_space/06_function_space_dofs_nodes_structured.h"
#include "field_variable/00_field_variable_base.h"
#include "mesh/type_traits.h"

namespace FunctionSpace
{

template<typename MeshType,typename BasisFunctionType>
class FunctionSpace;

/** class to get node positions
 */
template<typename MeshType,typename BasisFunctionType>
class FunctionSpaceDofsNodes
{
};

/** Partial specialization for RegularFixed mesh
 */
template<int D,typename BasisFunctionType>
class FunctionSpaceDofsNodes<Mesh::StructuredRegularFixedOfDimension<D>,BasisFunctionType> :
  public FunctionSpaceDofsNodesStructured<Mesh::StructuredRegularFixedOfDimension<D>,BasisFunctionType>,
  public std::enable_shared_from_this<FunctionSpaceDofsNodes<Mesh::StructuredRegularFixedOfDimension<D>,BasisFunctionType>>
{
public:

  //! constructor from python settings
  FunctionSpaceDofsNodes(std::shared_ptr<Partition::Manager> partitionManager, PythonConfig specificSettings);

  //! construct from element numbers and physical extent
  //! @param inputMeshIsGlobal if the number of elements in nElements is considered to be the global number of elements
  FunctionSpaceDofsNodes(std::shared_ptr<Partition::Manager> partitionManager, std::array<element_no_t, D> nElements,
                         std::array<double, D> physicalExtent, const std::array<int,D> nRanksPerCoordinateDirection, bool inputMeshIsGlobal=true);

  //! constructor from python settings, null argument is ignored
  FunctionSpaceDofsNodes(std::shared_ptr<Partition::Manager> partitionManager, std::vector<double> &null, PythonConfig specificSettings);

  //! construct from element numbers and physical extent, null argument is ignored
  FunctionSpaceDofsNodes(std::shared_ptr<Partition::Manager> partitionManager, std::vector<double> &null,
                         std::array<element_no_t, D> nElements, std::array<double, D> physicalExtent,
                         const std::array<int,D> nRanksPerCoordinateDirection, bool inputMeshIsGlobal=true);

  typedef FieldVariable::FieldVariable<FunctionSpace<Mesh::StructuredRegularFixedOfDimension<D>,BasisFunctionType>,3> GeometryFieldType;  //< the class typename of the geometry field variable

  //! get mesh width (=distance between nodes) of the given coordinate direction
  double meshWidth() const;

  //! initialize geometry
  virtual void initialize();
  
protected:
  
  //! compute the meshWidth_ from physicalExtent_
  void computeMeshWidth();
  
  std::array<double,D> physicalExtent_;   //< geometrical "size" of the mesh, i.e. length x width x height
  double meshWidth_;   //< uniform mesh width, i.e. distance between nodes (not elements for quadratic element), this is a copy of the value which is stored in this->geometryField_
};

/** Partial specialization for StructuredDeformable mesh
 */
template<int D,typename BasisFunctionType>
class FunctionSpaceDofsNodes<Mesh::StructuredDeformableOfDimension<D>,BasisFunctionType> :
  public FunctionSpaceDofsNodesStructured<Mesh::StructuredDeformableOfDimension<D>,BasisFunctionType>,
  public std::enable_shared_from_this<FunctionSpaceDofsNodes<Mesh::StructuredDeformableOfDimension<D>,BasisFunctionType>>
{
public:
  //! constructor from python settings, it is possible to create a basisOnMesh object without geometry field, e.g. for the lower order mesh of a mixed formulation
  FunctionSpaceDofsNodes(std::shared_ptr<Partition::Manager> partitionManager, PythonConfig specificSettings, bool noGeometryField=false);

  //! constructor from python settings with additionally given node positions
  FunctionSpaceDofsNodes(std::shared_ptr<Partition::Manager> partitionManager, std::vector<double> &nodePositions, PythonConfig specificSettings,
                         bool noGeometryField=false);

  //! constructor from node positions, nElementsPerCoordinateDirection are the local elements
  FunctionSpaceDofsNodes(std::shared_ptr<Partition::Manager> partitionManager, const std::vector<Vec3> &nodePositions,
                         const std::array<element_no_t,D> nElementsPerCoordinateDirection,
                         const std::array<int,D> nRanksPerCoordinateDirection, bool hasTriangleCorners);

  //! constructor from node positions (dummy, needed for creation from node positions from file)
  FunctionSpaceDofsNodes(std::shared_ptr<Partition::Manager> partitionManager, const std::vector<double> &nodePositionsFromBinaryFile, const std::vector<Vec3> &nodePositions,
                         const std::array<element_no_t,D> nElementsPerCoordinateDirection,
                         const std::array<int,D> nRanksPerCoordinateDirection, bool hasTriangleCorners);

  typedef FieldVariable::FieldVariable<FunctionSpace<Mesh::StructuredDeformableOfDimension<D>,BasisFunctionType>,3> GeometryFieldType;  //< the class typename of the geometry field variable

  //! initialize geometry
  virtual void initialize();

  //! refine the mesh by given factor, create new node positions
  void refineMesh(std::array<int,D> refinementFactors);

  //! set the dependent dofs in the given field variable by interpolating the independent dofs of the triangle basis, this implementation does nothing
  virtual void interpolateNonDofValuesInFieldVariable(
    std::shared_ptr<FieldVariable::FieldVariableBaseFunctionSpace<FunctionSpace<Mesh::StructuredDeformableOfDimension<D>,BasisFunctionType>>> fieldVariable,
    int componentNo) const = 0;

protected:

  //! parse the node from python config into a vector
  void parseNodePositionsFromSettings(PythonConfig specificSettings);

  //! set the values of the geometry field
  void setGeometryFieldValues();
  
  std::vector<double> localNodePositions_; //< Node positions to be inserted into geometry field, for local domain
 
};

/** Partial specialization for UnstructuredDeformable mesh
 */
template<int D,typename BasisFunctionType>
class FunctionSpaceDofsNodes<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType> :
  public FunctionSpaceGeometry<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>
{
public:
  //! inherit constructor
  using FunctionSpaceGeometry<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>::FunctionSpaceGeometry;

  //! fill a vector with the node position entries, nodes will contain consecutively the (x,y,z) values of just all nodes, i.e. for Hermite not the derivatives
  void getNodePositions(std::vector<double> &nodes) const;

  //! return number of nodes including ghost nodes, i.e. these nodes are known locally but some of them are owned by other ranks
  node_no_t nNodesLocalWithGhosts() const;

  //! return number of nodes that are owned by this partition
  node_no_t nNodesLocalWithoutGhosts() const;

  //! return number of dofs
  dof_no_t nDofsLocalWithGhosts() const;

  //! return number of dofs
  dof_no_t nDofsLocalWithoutGhosts() const;

  //! return global number of nodes
  global_no_t nNodesGlobal() const;

  //! return global number of dofs
  global_no_t nDofsGlobal() const;

};

/** Partial specialization for composite mesh
 */
template<int D,typename BasisFunctionType>
class FunctionSpaceDofsNodes<Mesh::CompositeOfDimension<D>,BasisFunctionType> :
  public FunctionSpaceGeometry<Mesh::CompositeOfDimension<D>,BasisFunctionType>,
  public std::enable_shared_from_this<FunctionSpaceDofsNodes<Mesh::CompositeOfDimension<D>,BasisFunctionType>>
{
public:

  using FunctionSpaceGeometry<Mesh::CompositeOfDimension<D>,BasisFunctionType>::FunctionSpaceGeometry;

  //! constructor from sub function spaces, it is possible to create a basisOnMesh object without geometry field, e.g. for the lower order mesh of a mixed formulation
  FunctionSpaceDofsNodes(std::shared_ptr<Partition::Manager> partitionManager,
                         std::vector<std::shared_ptr<FunctionSpace<Mesh::StructuredDeformableOfDimension<D>,BasisFunctionType>>> subFunctionSpaces, bool noGeometryField=false);

  //! constructor with PythonConfig, this is needed such that mesh manager compiles but never called
  FunctionSpaceDofsNodes(std::shared_ptr<Partition::Manager> partitionManager,
                         PythonConfig settings, bool noGeometryField=false);

  //! constructor with additionally given node positions, this is needed such that mesh manager compiles but never called
  FunctionSpaceDofsNodes(std::shared_ptr<Partition::Manager> partitionManager, std::vector<double> &nodePositions,
                         PythonConfig settings, bool noGeometryField=false);

  //! initialize geometry
  virtual void initialize();

  //! fill a vector with the node position entries, nodes will contain consecutively the (x,y,z) values of just all nodes, i.e. for Hermite not the derivatives
  void getNodePositions(std::vector<double> &nodes) const;

  //! return number of nodes including ghost nodes, i.e. these nodes are known locally but some of them are owned by other ranks
  node_no_t nNodesLocalWithGhosts() const;

  //! return number of nodes that are owned by this partition
  node_no_t nNodesLocalWithoutGhosts() const;

  //! return number of dofs
  dof_no_t nDofsLocalWithGhosts() const;

  //! return number of dofs
  dof_no_t nDofsLocalWithoutGhosts() const;

  //! return global number of nodes
  global_no_t nNodesGlobal() const;

  //! return global number of dofs
  global_no_t nDofsGlobal() const;

};

}  // namespace

#include "function_space/06_function_space_dofs_nodes_regular_fixed.tpp"
#include "function_space/06_function_space_dofs_nodes_structured_deformable.tpp"
#include "function_space/06_function_space_dofs_nodes_unstructured_deformable.tpp"
#include "function_space/06_function_space_dofs_nodes_composite.tpp"
