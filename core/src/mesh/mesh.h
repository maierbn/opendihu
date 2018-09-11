#pragma once

#include <Python.h>  // has to be the first included header
#include <control/types.h>
#include <memory>
#include "partition/00_mesh_partition_base.h"

namespace Mesh
{

/**
 * A mesh object contains the topological information and size in terms of number of nodes.
 * Also the geometrical information for each element is stored.
 * The number of DOFs per node is not stored here, instead at data_management classes.
 * Only structured quadrilateral meshes are considered.
 */
class Mesh
{
public:
  //! construct mesh from python settings
  Mesh(PyObject *specificSettings);
  virtual ~Mesh() {}
  
  //! dimensionality of the mesh
  virtual int dimension() const = 0;
  
  //! number of nodes in the mesh stored in the local partition, this also includes ghost nodes
  virtual node_no_t nNodesLocalWithGhosts() const = 0;
  
  //! the number of non-ghost nodes stored in the local partition
  virtual node_no_t nNodesLocalWithoutGhosts() const = 0;
  
  //! number of elements in the mesh stored in the current partition
  virtual element_no_t nElementsLocal() const = 0;
  
  //! get the meshPartition of this mesh (defined in FunctionSpace)
  virtual std::shared_ptr<Partition::MeshPartitionBase> meshPartitionBase() = 0;
  
  //! initialize the mesh after creation
  virtual void initialize() = 0;
  
  //! set the name of the mesh
  void setMeshName(std::string meshName);
  
  //! get the name of the mesh
  std::string meshName();
protected:
  std::string meshName_;  ///< the name of this mesh, which can be given in the python config and is the key by which the mesh is stored in Mesh::Manager
  PyObject *specificSettings_;  ///< the python object of the settings for this mesh
};

/**
 * base class for a mesh with a dimension.
 */
template<int D>
class MeshOfDimension : public Mesh
{
public:
  //! construct mesh from python settings
  MeshOfDimension(PyObject *specificSettings);
  virtual ~MeshOfDimension() {}

  //! return the dimension/template argument D as method. This is also accessible from the base class.
  int dimension() const;

  //! return the template argument D (the dimension) as constexpr
  static constexpr int dim();

  //! number of nodes in the mesh stored in the local partition, this also includes ghost nodes
  virtual node_no_t nNodesLocalWithGhosts() const = 0;
  
  //! the number of non-ghost nodes stored in the local partition
  virtual node_no_t nNodesLocalWithoutGhosts() const = 0;
  
  //! number of elements in the mesh stored in the current partition
  virtual element_no_t nElementsLocal() const = 0;
  
};

}  // namespace

#include "mesh/mesh.tpp"
