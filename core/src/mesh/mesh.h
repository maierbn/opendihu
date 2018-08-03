#pragma once

#include <Python.h>  // has to be the first included header
#include <control/types.h>
#include <memory>

// forward declaration
namespace Partition
{
class MeshPartitionBase;
}

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
  
  //! number of nodes in the mesh
  virtual node_no_t nLocalNodes() const = 0;
  
  //! number of elements in the mesh
  virtual element_no_t nLocalElements() const = 0;
  
  //! get the meshPartition of this mesh (defined in BasisOnMesh)
  virtual std::shared_ptr<Partition::MeshPartitionBase> meshPartitionBase() = 0;
  
  //! initialize the mesh after creation
  virtual void initialize() = 0;
  
  //! set the name of the mesh
  void setMeshName(std::string meshName);
  
  //! get the name of the mesh
  std::string meshName();
protected:
  std::string meshName_;  ///< the name of this mesh, which can be given in the python config and is the key by which the mesh is stored in Mesh::Manager
};

/** dummy mesh to signal that no mesh was specified (meshManager will instead create a mesh with a single element)
 */
class None : public Mesh
{
public:
  using Mesh::Mesh;
  //! dimensionality of the mesh
  int dimension() const {return 0;}
  static constexpr int dim() {return 0;}
  
  //! number of nodes in the mesh
  node_no_t nLocalNodes() const {return 0;}
  
  //! number of elements in the mesh
  element_no_t nLocalElements() const {return 0;}
  
  //! initialization method
  void initialize(){}
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

  //! get the number of nodes of this mesh
  virtual node_no_t nLocalNodes() const = 0;
  
  //! number of elements in the mesh
  virtual element_no_t nLocalElements() const = 0;
};

}  // namespace

#include "mesh/mesh.tpp"
