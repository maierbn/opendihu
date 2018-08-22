#pragma once

#include <Python.h>  // has to be the first included header

#include <array>
#include "control/types.h"

#include "basis_on_mesh/05_basis_on_mesh_geometry.h"
#include "basis_on_mesh/06_basis_on_mesh_dofs_nodes_structured.h"
#include "mesh/type_traits.h"

namespace BasisOnMesh
{

template<typename MeshType,typename BasisFunctionType>
class BasisOnMesh;

/** class to get node positions
 */
template<typename MeshType,typename BasisFunctionType>
class BasisOnMeshDofsNodes
{
protected:
};

/** Partial specialization for RegularFixed mesh
 */
template<int D,typename BasisFunctionType>
class BasisOnMeshDofsNodes<Mesh::StructuredRegularFixedOfDimension<D>,BasisFunctionType> :
  public BasisOnMeshDofsNodesStructured<Mesh::StructuredRegularFixedOfDimension<D>,BasisFunctionType>,
  public std::enable_shared_from_this<BasisOnMeshDofsNodes<Mesh::StructuredRegularFixedOfDimension<D>,BasisFunctionType>>
{
public:

  //! constructor from python settings
  BasisOnMeshDofsNodes(std::shared_ptr<Partition::Manager> partitionManager, PyObject *specificSettings);

  //! construct from element numbers and physical extent
  BasisOnMeshDofsNodes(std::shared_ptr<Partition::Manager> partitionManager, std::array<element_no_t, D> nElements, std::array<double, D> physicalExtent);

  typedef FieldVariable::FieldVariable<BasisOnMesh<Mesh::StructuredRegularFixedOfDimension<D>,BasisFunctionType>,3> GeometryFieldType;  ///< the class typename of the geometry field variable

  //! get mesh width (=distance between nodes) of the given coordinate direction
  double meshWidth() const;

  //! initialize geometry
  virtual void initialize();
  
protected:
  
  //! compute the meshWidth_ from physicalExtent_
  void computeMeshWidth();
  
  std::array<double,D> physicalExtent_;   ///< geometrical "size" of the mesh, i.e. length x width x height
  double meshWidth_;   ///< uniform mesh width, i.e. distance between nodes (not elements for quadratic element), this is a copy of the value which is stored in this->geometryField_
};

/** Partial specialization for StructuredDeformable mesh
 */
template<int D,typename BasisFunctionType>
class BasisOnMeshDofsNodes<Mesh::StructuredDeformableOfDimension<D>,BasisFunctionType> :
  public BasisOnMeshDofsNodesStructured<Mesh::StructuredDeformableOfDimension<D>,BasisFunctionType>,
  public std::enable_shared_from_this<BasisOnMeshDofsNodes<Mesh::StructuredDeformableOfDimension<D>,BasisFunctionType>>
{
public:
  //! constructor from python settings, it is possible to create a basisOnMesh object without geometry field, e.g. for the lower order mesh of a mixed formulation
  BasisOnMeshDofsNodes(std::shared_ptr<Partition::Manager> partitionManager, PyObject *specificSettings, bool noGeometryField=false);

  //! constructor from node positions
  BasisOnMeshDofsNodes(std::shared_ptr<Partition::Manager> partitionManager, const std::vector<Vec3> &nodePositions, const std::array<element_no_t,D> nElementsPerCoordinateDirection);

  typedef FieldVariable::FieldVariable<BasisOnMesh<Mesh::StructuredDeformableOfDimension<D>,BasisFunctionType>,3> GeometryFieldType;  ///< the class typename of the geometry field variable
  
  //! initialize geometry
  virtual void initialize();
  
protected:

  //! parse the node from python config into a vector
  void parseNodePositionsFromSettings(PyObject *specificSettings);

  //! set the values of the geometry field
  void setGeometryFieldValues();
  
  std::vector<double> localNodePositions_; //< Node positions to be inserted into geometry field, for local domain
 
};

/** Partial specialization for UnstructuredDeformable mesh
 */
template<int D,typename BasisFunctionType>
class BasisOnMeshDofsNodes<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType> :
  public BasisOnMeshGeometry<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>
{
public:
  //! inherit constructor
  using BasisOnMeshGeometry<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>::BasisOnMeshGeometry;

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

#include "basis_on_mesh/06_basis_on_mesh_dofs_nodes_regular_fixed.tpp"
#include "basis_on_mesh/06_basis_on_mesh_dofs_nodes_structured_deformable.tpp"
#include "basis_on_mesh/06_basis_on_mesh_dofs_nodes_unstructured_deformable.tpp"
