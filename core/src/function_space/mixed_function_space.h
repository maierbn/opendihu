#pragma once

#include <Python.h>  // has to be the first included header

#include <array>
#include "control/types.h"
#include "mesh/mesh.h"
#include "basis_function/mixed.h"

namespace FunctionSpace
{

/** Mesh with two different node sets as overlay, used for mixed formulations. It inherits from Mesh to be able to store it as std::shared_ptr<Mesh::Mesh> in mesh manager.
 */
template<typename LowOrderFunctionSpaceType,typename HighOrderFunctionSpaceType>
class Mixed : public Mesh::Mesh
{
public:
  typedef typename HighOrderFunctionSpaceType::Mesh Mesh;   // the mesh must be the same for lowOrder and highOrder, only the FunctionSpace is different
  typedef BasisFunction::Mixed<typename LowOrderFunctionSpaceType::BasisFunction, typename HighOrderFunctionSpaceType::BasisFunction> BasisFunction;
  typedef LowOrderFunctionSpaceType LowOrderFunctionSpace;
  typedef HighOrderFunctionSpaceType HighOrderFunctionSpace;

  // dimension, which is the same for the LowOrderFunctionSpaceType and the HighOrderFunctionSpaceType
  static constexpr int dim();

  //! contructor
  Mixed(PythonConfig specificSettings);

  //! this assigns the geometry_ field variable's mesh pointer to this object, it is not possible from the constructor, therefore this extra method
  void initialize();

  //! the dimension of the mesh
  int dimension() const;

  //! number of nodes of the mesh (high order version)
  node_no_t nNodesLocalWithGhosts() const;

  //! return the low order basisOnMesh object
  std::shared_ptr<LowOrderFunctionSpaceType> lowOrderFunctionSpace();

  //! return the high order basisOnMesh object
  std::shared_ptr<HighOrderFunctionSpaceType> highOrderFunctionSpace();

private:
  std::shared_ptr<LowOrderFunctionSpaceType> lowOrderFunctionSpace_;     //< the FunctionSpace object for the low order basis
  std::shared_ptr<HighOrderFunctionSpaceType> highOrderFunctionSpace_;   //< the FunctionSpace object for the high order basis
};

} // namespace

#include "function_space/mixed_function_space.tpp"
