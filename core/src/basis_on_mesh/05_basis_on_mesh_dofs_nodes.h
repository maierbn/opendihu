#pragma once

#include <Python.h>  // has to be the first included header

#include <array>
#include "control/types.h"

#include "basis_on_mesh/04_basis_on_mesh_geometry.h"
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
};

/** Partial specialization for RegularFixed mesh
 */
template<int D,typename BasisFunctionType>
class BasisOnMeshDofsNodes<Mesh::StructuredRegularFixedOfDimension<D>,BasisFunctionType> :
  public BasisOnMeshGeometry<Mesh::StructuredRegularFixedOfDimension<D>,BasisFunctionType>
{
public:

  //! constructor from python settings
  BasisOnMeshDofsNodes(PyObject *specificSettings);

  //! construct from element numbers and physical extent
  BasisOnMeshDofsNodes(std::array<element_no_t, D> nElements, std::array<double, D> physicalExtent);

  typedef FieldVariable::FieldVariable<BasisOnMesh<Mesh::StructuredRegularFixedOfDimension<D>,BasisFunctionType>,3> GeometryFieldType;  ///< the class typename of the geometry field variable

  //! fill a vector with the node position entries, nodes will contain consecutively the (x,y,z) values of just all nodes, i.e. for Hermite not the derivatives
  void getNodePositions(std::vector<double> &nodes) const;

  //! get mesh width (=distance between nodes) of the given coordinate direction
  double meshWidth() const;

  //! return number of nodes
  node_no_t nNodes() const;

  //! return number of nodes in specified coordinate direction
  node_no_t nNodes(int dimension) const;

  //! return number of dofs
  dof_no_t nDofs() const;

protected:

  //! create the geometry field from meshWidth_
  void setupGeometryField();

  double meshWidth_;   ///< uniform mesh width, i.e. distance between nodes (not elements for quadratic element), this is a copy of the value which is stored in this->geometryField_
};

/** Partial specialization for StructuredDeformable mesh
 */
template<int D,typename BasisFunctionType>
class BasisOnMeshDofsNodes<Mesh::StructuredDeformableOfDimension<D>,BasisFunctionType> :
  public BasisOnMeshGeometry<Mesh::StructuredDeformableOfDimension<D>,BasisFunctionType>
{
public:
  //! constructor from python settings, it is possible to create a basisOnMesh object without geometry field, e.g. for the lower order mesh of a mixed formulation
  BasisOnMeshDofsNodes(PyObject *specificSettings, bool noGeometryField=false);

  //! constructor from node positions
  BasisOnMeshDofsNodes(const std::vector<Vec3> &nodePositions, const std::array<element_no_t,D> nElementsPerCoordinateDirection);

  typedef FieldVariable::FieldVariable<BasisOnMesh<Mesh::StructuredDeformableOfDimension<D>,BasisFunctionType>,3> GeometryFieldType;  ///< the class typename of the geometry field variable

  //! fill a vector with the node position entries, nodes will contain consecutively the (x,y,z) values of just all nodes, i.e. for Hermite not the derivatives
  void getNodePositions(std::vector<double> &nodes) const;

  //! return number of nodes
  node_no_t nNodes() const;

  //! return number of nodes in specified coordinate direction
  node_no_t nNodes(int dimension) const;

  //! return number of dofs
  dof_no_t nDofs() const;

  //!This creates geometry field
  void initialize() override;
  
protected:

  //! parse the node from python config into a vector
  void parseNodePositionsFromSettings(PyObject *specificSettings);

  //! set up the geometry field
  void initializeGeometryField();
  
  //! set the values of the geometry field
  void setGeometryFieldValues();
  
  
  std::vector<double> nodePositions_; //< Node positions to be inserted into geometry field
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

  //! return number of nodes
  node_no_t nNodes() const;

  //! return number of dofs
  dof_no_t nDofs() const;
};

}  // namespace

#include "basis_on_mesh/05_basis_on_mesh_dofs_nodes_regular.tpp"
#include "basis_on_mesh/05_basis_on_mesh_dofs_nodes_structured.tpp"
#include "basis_on_mesh/05_basis_on_mesh_dofs_nodes_unstructured.tpp"
