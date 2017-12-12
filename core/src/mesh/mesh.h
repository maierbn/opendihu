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
  virtual int dimension() = 0;
  virtual element_idx_t nNodes() = 0;
private:
};

template<unsigned long D>
class MeshD : public Mesh
{
public:
  //! construct mesh from python settings
  MeshD(PyObject *specificSettings);
  virtual ~MeshD() {}
  virtual int dimension();
  virtual element_idx_t nNodes() = 0;
private:
};

}  // namespace

#include "mesh/mesh.tpp"