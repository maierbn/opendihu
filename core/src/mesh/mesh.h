#pragma once

#include <Python.h>
#include <control/types.h>

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
  virtual int dimension() const = 0;
  virtual element_idx_t nNodes() const = 0;
private:
};

/**
 * base class for a mesh with a dimension.
 */
template<unsigned long D>
class MeshD : public Mesh
{
public:
  //! construct mesh from python settings
  MeshD(PyObject *specificSettings);
  virtual ~MeshD() {}
  
  //! return the dimension/template argument D as method. This is also accessible from the base class.
  int dimension() const; 
  
  //! return the template argument D (the dimension) as constexpr
  static constexpr int dim();
  
  //! get the number of nodes of this mesh
  virtual element_idx_t nNodes() const = 0;
private:
};

}  // namespace

#include "mesh/mesh.tpp"