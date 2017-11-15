#pragma once

#include <Python.h>
#include <control/types.h>

namespace Mesh
{

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