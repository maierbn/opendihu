#pragma once

#include <Python.h>  // has to be the first included header

#include <array>
#include "control/types.h"

#include "basis_on_mesh/00_basis_on_mesh_base_dim.h"
#include "mesh/mesh.h"

namespace BasisOnMesh
{
template<typename MeshType,typename BasisFunctionType>
class BasisOnMeshFunction : 
  public BasisOnMeshBaseDim<MeshType::dim(),BasisFunctionType>,
  public MeshType
{
public:
  using MeshType::MeshType;
 
  //! evaluate the basis function corresponding to element-local dof dofIndex at xi, interval for xi is [0,1]
  static double phi(int dofIndex, std::array<double,MeshType::dim()> xi);
  
  //! evaluate the derivative of Phi(xi) w.r.t xi_i, where i is given by derivativeIdx
  static double dPhidxi(int dofIndex, int derivativeIdx, std::array<double,MeshType::dim()> xi);

  //! evaluate the first derivative of the basis function corresponding to element-local dof dofIndex at xi, interval for xi is [0,1]
  static std::array<double,MeshType::dim()> gradPhi(int dofIndex, std::array<double,MeshType::dim()> xi);

private:
 
  //! given an element-local dofIndex and dimension No (0 <= dimNo < D), return the basis function index in that direction
  static int getBasisFunctionIndex1D(int dofIndex, int dimNo);
};


}  // namespace

#include "basis_on_mesh/01_basis_on_mesh_function.tpp"
