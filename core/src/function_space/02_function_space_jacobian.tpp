#include "function_space/02_function_space_jacobian.h"

#include "utility/math_utility.h"
#include "utility/vector_operators.h"
#include "easylogging++.h"

#include <cmath>
#include <array>

namespace FunctionSpace
{

//! compute the (geometry) jacobian matrix, geometryField is the node positions for Lagrange basis, node positions and derivatives for Hermite basis
template<typename MeshType,typename BasisFunctionType,typename dummy>
template<typename Vec3>
std::array<Vec3,MeshType::dim()> FunctionSpaceJacobian<MeshType,BasisFunctionType,dummy>::
computeJacobian(const std::array<Vec3,FunctionSpaceFunction<MeshType,BasisFunctionType>::nDofsPerElement()> &geometryField,
                                                        const std::array<double,MeshType::dim()> xi)
{
  std::array<Vec3,MeshType::dim()> jacobian;
  // loop over columns
  for (int dimNo = 0; dimNo < MeshType::dim(); dimNo++)
  {
    jacobian[dimNo] = Vec3({0.0});
    for (int dofIndex = 0; dofIndex < FunctionSpaceFunction<MeshType,BasisFunctionType>::nDofsPerElement(); dofIndex++)
    {
      double coefficient = FunctionSpaceFunction<MeshType,BasisFunctionType>::dphi_dxi(dofIndex, dimNo, xi);
      jacobian[dimNo] += coefficient * geometryField[dofIndex];

      if (VLOG_IS_ON(3))
        VLOG(3) << "   col " << dimNo << " dof " << dofIndex << ", coeff: " << coefficient << ", node " << geometryField[dofIndex]
          << " -> " << jacobian[dimNo];
    }
  }

// check for singularity
#ifndef NDEBUG

  // check if jacobian contains column with all zeros, then output a warning
  for (int dimNo = 0; dimNo < MeshType::dim(); dimNo++)
  {
    if (MathUtility::template normSquared<3>(jacobian[dimNo]) < 1e-12)
    {
      LOG(WARNING) << "Jacobian " << jacobian << " is singular (column " << dimNo << "), xi: " << xi << ", geometryField: " << geometryField;
      LOG(DEBUG) << "Enable debugging output with -vmodule=*jacobian*=3";
      if (std::is_same<BasisFunctionType,BasisFunction::Hermite>::value)
      {
        LOG(DEBUG) << "You are using Hermite polynomials, check if geometry is specified correctly using also the derivative dofs!";
      }
      break;
    }
  }
#endif

  return jacobian;
}

// specialization: Jacobian for 1D linear Lagrange basis
template<typename MeshType>
template<typename Vec3>
std::array<Vec3,1> FunctionSpaceJacobian<MeshType,BasisFunction::LagrangeOfOrder<1>,Mesh::isDim<1,MeshType>>::
computeJacobian(const std::array<Vec3,FunctionSpaceFunction<MeshType,BasisFunction::LagrangeOfOrder<1>>::nDofsPerElement()> &node, const std::array<double,1> xi)
{
  Vec3 jacobianColumn0 = (node[1]-node[0]);
  return std::array<Vec3,1>({jacobianColumn0});
}

// specialization: Jacobian for 2D linear Lagrange basis
template<typename MeshType>
template<typename Vec3>
std::array<Vec3,2> FunctionSpaceJacobian<MeshType,BasisFunction::LagrangeOfOrder<1>,Mesh::isDim<2,MeshType>>::
computeJacobian(const std::array<Vec3,FunctionSpaceFunction<MeshType,BasisFunction::LagrangeOfOrder<1>>::nDofsPerElement()> &node, const std::array<double,2> xi)
{
  double xi1 = xi[0];
  double xi2 = xi[1];
  Vec3 jacobianColumn0 = (1-xi2) * (node[1]-node[0]) + xi2 * (node[3]-node[2]);
  Vec3 jacobianColumn1 = (1-xi1) * (node[2]-node[0]) + xi1 * (node[3]-node[1]);
  VLOG(3) << "computeJacobian for (" <<xi1<< "," <<xi2<< ")";
  return std::array<Vec3,2>({jacobianColumn0, jacobianColumn1});
}

// specialization: Jacobian for 3D linear Lagrange basis
template<typename MeshType>
template<typename Vec3>
std::array<Vec3,3> FunctionSpaceJacobian<MeshType,BasisFunction::LagrangeOfOrder<1>,Mesh::isDim<3,MeshType>>::
computeJacobian(const std::array<Vec3,FunctionSpaceFunction<MeshType,BasisFunction::LagrangeOfOrder<1>>::nDofsPerElement()> &node, const std::array<double,3> xi)
{
  double xi1 = xi[0];
  double xi2 = xi[1];
  double xi3 = xi[2];

  Vec3 jacobianColumn0 =
    (1-xi2) * (1-xi3) * (node[1]-node[0])
    + xi2 * (1-xi3) * (node[3]-node[2])
    + (1-xi2) * xi3 * (node[5]-node[4])
    + xi2 * xi3 * (node[7]-node[6]);

  Vec3 jacobianColumn1 =
    (1-xi1) * (1-xi3) * (node[2]-node[0])
    + xi1 * (1-xi3) * (node[3]-node[1])
    + (1-xi1) * xi3 * (node[6]-node[4])
    + xi1 * xi3 * (node[7]-node[5]);

  Vec3 jacobianColumn2 =
    (1-xi1) * (1-xi2) * (node[4]-node[0])
    + xi1 * (1-xi2) * (node[5]-node[1])
    + (1-xi1) * xi2 * (node[6]-node[2])
    + xi1 * xi2 * (node[7]-node[3]);

  return std::array<Vec3,3>({jacobianColumn0, jacobianColumn1, jacobianColumn2});
}
/*
// general implementation of Jacobian
{
  VLOG(3) << "computeJacobian generic for " <<xi;
  std::array<Vec3,MeshType::dim()> jacobian;
  for (int dimNo = 0; dimNo < MeshType::dim(); dimNo++)
  {
    jacobian[dimNo] = Vec3({0.0});
    for (int dofIndex = 0; dofIndex < FunctionSpaceFunction<MeshType,BasisFunctionType>::nDofsPerElement(); dofIndex++)
    {
      double coefficient = FunctionSpaceFunction<MeshType,BasisFunctionType>::dphi_dxi(dofIndex, dimNo, xi);
      jacobian[dimNo] += coefficient * geometryField[dofIndex];
      VLOG(3) << "   col " << dimNo << " dof " << dofIndex << ", coeff: " << coefficient << ", node " << geometryField[dofIndex]
       << " -> " << jacobian[dimNo];
    }
  }
  return jacobian;
}
*/
} // namespace
