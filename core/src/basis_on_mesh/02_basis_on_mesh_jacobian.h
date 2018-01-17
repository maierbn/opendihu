#pragma once

#include <array>
#include "control/types.h"

#include "basis_on_mesh/01_basis_on_mesh_function.h"
#include "basis_function/lagrange.h"
#include "mesh/type_traits.h"

namespace BasisOnMesh
{

// general algorithm to compute jacobian from basis functions
template<typename MeshType,typename BasisFunctionType,typename dummy = MeshType>
class BasisOnMeshJacobian :
  public BasisOnMeshFunction<MeshType,BasisFunctionType>
{
public:
  //! inherit constructor
  using BasisOnMeshFunction<MeshType,BasisFunctionType>::BasisOnMeshFunction;
  
  //! compute the jacobian matrix, geometryField is the node positions for Lagrange basis, node positions and derivatives for Hermite basis
  static std::array<Vec3,MeshType::dim()> computeJacobian(const std::array<Vec3,BasisOnMeshFunction<MeshType,BasisFunctionType>::nDofsPerElement()> &geometryField,
                                                          const std::array<double,MeshType::dim()> xi)
  {
    VLOG(3) << "computeJacobian generic for "<<xi;
    std::array<Vec3,MeshType::dim()> jacobian;
    for(int dimNo = 0; dimNo < MeshType::dim(); dimNo++)
    {
      jacobian[dimNo] = Vec3({0.0});
      for(int dofIndex = 0; dofIndex < BasisOnMeshFunction<MeshType,BasisFunctionType>::nDofsPerElement(); dofIndex++)
      {
        double coefficient = BasisOnMeshFunction<MeshType,BasisFunctionType>::dPhidxi(dofIndex, dimNo, xi);
        jacobian[dimNo] += coefficient * geometryField[dofIndex];
        VLOG(3) << "   col " << dimNo << " dof " << dofIndex << ", coeff: " << coefficient << ", node " << geometryField[dofIndex] 
         << " -> " << jacobian[dimNo];
      }
    }
    return jacobian;
  }

};

// partial specialization for linear Lagrange, D=1
template<typename MeshType>
class BasisOnMeshJacobian<MeshType,BasisFunction::Lagrange<1>,Mesh::isDim<1,MeshType>> :
  public BasisOnMeshFunction<MeshType,BasisFunction::Lagrange<1>>
{
public:
  //! inherit constructor
  using BasisOnMeshFunction<MeshType,BasisFunction::Lagrange<1>>::BasisOnMeshFunction;
  
  //! compute the jacobian matrix, geometryField is the node positions for Lagrange basis, node positions and derivatives for Hermite basis
  static std::array<Vec3,1> computeJacobian(const std::array<Vec3,BasisOnMeshFunction<MeshType,BasisFunction::Lagrange<1>>::nDofsPerElement()> &geometryField,
                                            const std::array<double,1> xi);
  
};

// partial specialization for linear Lagrange, D=2
template<typename MeshType>
class BasisOnMeshJacobian<MeshType,BasisFunction::Lagrange<1>,Mesh::isDim<2,MeshType>> :
  public BasisOnMeshFunction<MeshType,BasisFunction::Lagrange<1>>
{
public:
  //! inherit constructor
  using BasisOnMeshFunction<MeshType,BasisFunction::Lagrange<1>>::BasisOnMeshFunction;
  
  //! compute the jacobian matrix, geometryField is the node positions for Lagrange basis, node positions and derivatives for Hermite basis
  static std::array<Vec3,2> computeJacobian(const std::array<Vec3,BasisOnMeshFunction<MeshType,BasisFunction::Lagrange<1>>::nDofsPerElement()> &geometryField,
                                            const std::array<double,2> xi);
  
};

// partial specialization for linear Lagrange, D=3
template<typename MeshType>
class BasisOnMeshJacobian<MeshType,BasisFunction::Lagrange<1>,Mesh::isDim<3,MeshType>> :
  public BasisOnMeshFunction<MeshType,BasisFunction::Lagrange<1>>
{
public:
  //! inherit constructor
  using BasisOnMeshFunction<MeshType,BasisFunction::Lagrange<1>>::BasisOnMeshFunction;
  
  //! compute the jacobian matrix, geometryField is the node positions for Lagrange basis, node positions and derivatives for Hermite basis
  static std::array<Vec3,3> computeJacobian(const std::array<Vec3,BasisOnMeshFunction<MeshType,BasisFunction::Lagrange<1>>::nDofsPerElement()> &geometryField,
                                            const std::array<double,3> xi);
  
};

}  // namespace

#include "basis_on_mesh/02_basis_on_mesh_jacobian.tpp"