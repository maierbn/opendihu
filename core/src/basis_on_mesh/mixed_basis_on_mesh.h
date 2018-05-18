#pragma once

#include <Python.h>  // has to be the first included header

#include <array>
#include "control/types.h"
#include "mesh/mesh.h"
#include "basis_function/mixed.h"

namespace BasisOnMesh
{

/** Mesh with two different node sets as overlay, used for mixed formulations. It inherits from Mesh to be able to store it as std::shared_ptr<Mesh::Mesh> in mesh manager.
 */
template<typename LowOrderBasisOnMeshType,typename HighOrderBasisOnMeshType>
class Mixed : public Mesh::Mesh
{
public:
  typedef typename HighOrderBasisOnMeshType::Mesh Mesh;   // the mesh must be the same for lowOrder and highOrder, only the BasisOnMesh is different
  typedef BasisFunction::Mixed<typename LowOrderBasisOnMeshType::BasisFunction, typename HighOrderBasisOnMeshType::BasisFunction> BasisFunction;
  typedef LowOrderBasisOnMeshType LowOrderBasisOnMesh;
  typedef HighOrderBasisOnMeshType HighOrderBasisOnMesh;

  // dimension, which is the same for the LowOrderBasisOnMeshType and the HighOrderBasisOnMeshType
  static constexpr int dim();

  //! contructor
  Mixed(PyObject *specificSettings);

  //! this assigns the geometry_ field variable's mesh pointer to this object, it is not possible from the constructor, therefore this extra method
  void initialize();

  //! the dimension of the mesh
  int dimension() const;

  //! number of nodes of the mesh (high order version)
  node_no_t nNodes() const;

  //! return the low order basisOnMesh object
  std::shared_ptr<LowOrderBasisOnMeshType> lowOrderBasisOnMesh();

  //! return the high order basisOnMesh object
  std::shared_ptr<HighOrderBasisOnMeshType> highOrderBasisOnMesh();

private:
  std::shared_ptr<LowOrderBasisOnMeshType> lowOrderBasisOnMesh_;     ///< the BasisOnMesh object for the low order basis
  std::shared_ptr<HighOrderBasisOnMeshType> highOrderBasisOnMesh_;   ///< the BasisOnMesh object for the high order basis
};

};  // namespace

#include "basis_on_mesh/mixed_basis_on_mesh.tpp"
