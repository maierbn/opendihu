#pragma once

#include <Python.h>  // has to be the first included header

#include <array>
#include "control/types.h"

#include "basis_on_mesh/03_basis_on_mesh_dofs.h"
#include "mesh/type_traits.h"

namespace BasisOnMesh
{
 
template<typename MeshType,typename BasisFunctionType>
class BasisOnMesh;

/** class to get node positions
 */
template<typename MeshType,typename BasisFunctionType>
class BasisOnMeshNodes
{
};

/** Partial specialization for RegularFixed mesh
 */
template<int D,typename BasisFunctionType>
class BasisOnMeshNodes<Mesh::StructuredRegularFixedOfDimension<D>,BasisFunctionType> :
  public BasisOnMeshDofs<Mesh::StructuredRegularFixedOfDimension<D>,BasisFunctionType>,
  public std::enable_shared_from_this<BasisOnMeshNodes<Mesh::StructuredRegularFixedOfDimension<D>,BasisFunctionType>>
{
public:
 
  //! inherit constructor
  //using BasisOnMeshDofs<Mesh::StructuredRegularFixedOfDimension<D>,BasisFunctionType>::BasisOnMeshDofs;
  //! constructor
  BasisOnMeshNodes(PyObject *specificSettings);
    
  //! construct from element numbers and physical extent
  BasisOnMeshNodes(std::array<element_no_t, D> nElements, std::array<double, D> physicalExtent);
  
  //! this assigns the geometry_ field variable's mesh pointer to this object, it is not possible from the constructor, therefore this extra method
  void initialize();
  
  //! fill a vector with the node position entries, nodes will contain consecutively the (x,y,z) values of just all nodes, i.e. for Hermite not the derivatives
  void getNodePositions(std::vector<double> &nodes) const;
  
  //! get mesh width of the given coordinate direction
  double meshWidth() const;
  
  //! return number of nodes
  node_no_t nNodes() const;
  
  //! return number of nodes in specified coordinate direction
  node_no_t nNodes(int dimension) const;
  
  //! return number of dofs
  node_no_t nDofs() const;
  
  //! return the geometry field entry (node position for Lagrange elements) of a specific dof
  Vec3 getGeometry(node_no_t dofGlobalNo) const;
  
  //! get all geometry entries for an element
  void getElementGeometry(element_no_t elementNo, std::array<Vec3, BasisOnMeshBaseDim<D,BasisFunctionType>::nDofsPerElement()> &values);
  
  typedef FieldVariable::FieldVariable<BasisOnMesh<Mesh::StructuredRegularFixedOfDimension<D>,BasisFunctionType>> FieldVariableType;  ///< the class typename of the geometry field variable
   
  //! return the internal geometry field variable
  FieldVariableType &geometryField();
  
  //! has no effect for structured meshes
  void addNonGeometryFieldVariables(std::vector<std::shared_ptr<FieldVariableType>> &fieldVariables){}

  //! if the geometry field is set
  bool hasGeometryField();
  
protected:
  
  double meshWidth_;   ///< uniform mesh width
 
  std::unique_ptr<FieldVariableType> geometry_;     ///< the geometry field variable
}; 

/** Partial specialization for StructuredDeformable mesh
 */
template<int D,typename BasisFunctionType>
class BasisOnMeshNodes<Mesh::StructuredDeformableOfDimension<D>,BasisFunctionType> :
  public BasisOnMeshDofs<Mesh::StructuredDeformableOfDimension<D>,BasisFunctionType>,
  public std::enable_shared_from_this<BasisOnMeshNodes<Mesh::StructuredDeformableOfDimension<D>,BasisFunctionType>>
{
public:
  //! constructor, it is possible to create a basisOnMesh object without geometry field, e.g. for the lower order mesh of a mixed formulation
  BasisOnMeshNodes(PyObject *specificSettings, bool noGeometryField=false);
  
  //! this assigns the geometry_ field variable's mesh pointer to this object, it is not possible from the constructor, therefore this extra method
  void initialize();
  
  //! fill a vector with the node position entries, nodes will contain consecutively the (x,y,z) values of just all nodes, i.e. for Hermite not the derivatives
  void getNodePositions(std::vector<double> &nodes) const;
  
  //! return number of nodes
  node_no_t nNodes() const;
  
  //! return number of nodes in specified coordinate direction
  node_no_t nNodes(int dimension) const;
  
  //! return number of dofs
  node_no_t nDofs() const;
  
  //! return the geometry field entry (node position for Lagrange elements) of a specific dof
  Vec3 getGeometry(node_no_t dofGlobalNo) const;
  
  //! get all geometry entries for an element
  void getElementGeometry(element_no_t elementNo, std::array<Vec3, BasisOnMeshBaseDim<D,BasisFunctionType>::nDofsPerElement()> &values);
 
  typedef FieldVariable::FieldVariable<BasisOnMesh<Mesh::StructuredDeformableOfDimension<D>,BasisFunctionType>> FieldVariableType;  ///< the class typename of the geometry field variable
 
  //! return the geometry field of this mesh
  FieldVariableType &geometryField();
  
  //! has no effect for structured meshes
  void addNonGeometryFieldVariables(std::vector<std::shared_ptr<FieldVariableType>> &fieldVariables){}
  
  //! if the geometry field is set
  bool hasGeometryField();
  
protected:
 
  //! parse the node from python config into a vector
  void parseNodePositionsFromSettings(PyObject *specificSettings, std::vector<double> &nodePositions);
  
  //! set the geometry field by the node positions in a vector
  void setGeometryField(std::vector<double> &nodePositions);
  
  std::unique_ptr<FieldVariableType> geometry_;     ///< the geometry field variable
  bool noGeometryField_;     ///< this is set if there is no geometry field stored. this is only needed for solid mechanics mixed formulation where the lower order basisOnMesh does not need its own geometry information
};

/** Partial specialization for UnstructuredDeformable mesh
 */
template<int D,typename BasisFunctionType>
class BasisOnMeshNodes<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType> :
  public BasisOnMeshDofs<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>
{
public:
  //! inherit constructor
  using BasisOnMeshDofs<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>::BasisOnMeshDofs;
  
  //! fill a vector with the node position entries, nodes will contain consecutively the (x,y,z) values of just all nodes, i.e. for Hermite not the derivatives
  void getNodePositions(std::vector<double> &nodes) const;
  
  //! return number of nodes
  node_no_t nNodes() const;
  
  //! return the geometry field entry (node position for Lagrange elements) of a specific dof
  Vec3 getGeometry(node_no_t dofGlobalNo) const;
  
  //! get all geometry entries for an element
  void getElementGeometry(element_no_t elementNo, std::array<Vec3, BasisOnMeshBaseDim<D,BasisFunctionType>::nDofsPerElement()> &values);
 
  typedef FieldVariable::FieldVariable<BasisOnMesh<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>> FieldVariableType;  ///< the type of a field variable on this mesh
 
  //! return a reference to this mesh' geometry field
  FieldVariableType &geometryField();
 
  //! if the geometry field is set
  bool hasGeometryField();
protected:
  // contains map<string, shared_ptr<FieldVariableType>> fieldVariable_
}; 

}  // namespace

#include "basis_on_mesh/04_basis_on_mesh_nodes_regular.tpp"
#include "basis_on_mesh/04_basis_on_mesh_nodes_structured.tpp"
#include "basis_on_mesh/04_basis_on_mesh_nodes_unstructured.tpp"
