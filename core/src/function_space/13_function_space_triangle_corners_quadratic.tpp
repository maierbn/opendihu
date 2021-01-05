#include "function_space/13_function_space_triangle_corners.h"

#include "partition/mesh_partition/01_mesh_partition_structured.h"
#include "function_space/function_space.h"

namespace FunctionSpace
{

//--------------------------
// quadratic
template<typename Dummy>
double FunctionSpaceTriangleCorners<Mesh::StructuredDeformableOfDimension<3>, BasisFunction::LagrangeOfOrder<2>,Dummy>::
phi(int dofIndex, std::array<double,3> xi, element_no_t elementNoLocal) const
{
  using BasisFunctionType = BasisFunction::LagrangeOfOrder<2>;
  ::Mesh::face_or_edge_t edge;

  // if the current element is a corner triangle
  if (elementNoLocal >= 0 && this->hasTriangleCorners_ && this->meshPartition_->elementIsAtCorner(elementNoLocal, edge))
  {
    int basisFunctionIndex1D = this->getBasisFunctionIndex1D(dofIndex, 2);

    const double xi1 = xi[0];
    const double xi2 = xi[1];
    const double xi1m = 1 - xi1;
    const double xi2m = 1 - xi2;

    // depending on the orientation of the triangle
    switch(edge)
    {
    // 0-1-
    case ::Mesh::face_or_edge_t::edge0Minus1Minus:
      switch (dofIndex)
      {
      case 0:
      case 9:
      case 18:
        return 4*xi1m*xi2m * BasisFunctionType::phi(basisFunctionIndex1D, xi[2]);

      case 1:
      case 10:
      case 19:
        return 0.0;

      case 2:
      case 11:
      case 20:
        return xi2m*(2*xi2m - 1) * BasisFunctionType::phi(basisFunctionIndex1D, xi[2]);

      case 3:
      case 12:
      case 21:
        return 0.0;

      case 4:
      case 13:
      case 22:
        return 0.0;

      case 5:
      case 14:
      case 23:
        return 4*xi2m*(1 - xi1m - xi2m) * BasisFunctionType::phi(basisFunctionIndex1D, xi[2]);

      case 6:
      case 15:
      case 24:
        return xi1m*(2*xi1m - 1) * BasisFunctionType::phi(basisFunctionIndex1D, xi[2]);

      case 7:
      case 16:
      case 25:
        return 4*xi1m*(1 - xi1m - xi2m) * BasisFunctionType::phi(basisFunctionIndex1D, xi[2]);

      case 8:
      case 17:
      case 26:
        return (xi1m + xi2m - 1)*(2*xi1m + 2*xi2m - 1) * BasisFunctionType::phi(basisFunctionIndex1D, xi[2]);
      }
      break;

    // 0+1-
    case ::Mesh::face_or_edge_t::edge0Plus1Minus:
      switch (dofIndex)
      {
      case 0:
      case 9:
      case 18:
        return xi2m*(2*xi2m - 1) * BasisFunctionType::phi(basisFunctionIndex1D, xi[2]);

      case 1:
      case 10:
      case 19:
        return 0.0;

      case 2:
      case 11:
      case 20:
        return 4*xi1*xi2m * BasisFunctionType::phi(basisFunctionIndex1D, xi[2]);

      case 3:
      case 12:
      case 21:
        return 4*xi2m*(1 - xi1 - xi2m) * BasisFunctionType::phi(basisFunctionIndex1D, xi[2]);

      case 4:
      case 13:
      case 22:
        return 0.0;

      case 5:
      case 14:
      case 23:
        return 0.0;

      case 6:
      case 15:
      case 24:
        return (xi1 + xi2m - 1)*(2*xi1 + 2*xi2m - 1) * BasisFunctionType::phi(basisFunctionIndex1D, xi[2]);

      case 7:
      case 16:
      case 25:
        return 4*xi1*(1 - xi1 - xi2m) * BasisFunctionType::phi(basisFunctionIndex1D, xi[2]);

      case 8:
      case 17:
      case 26:
        return xi1*(2*xi1 - 1) * BasisFunctionType::phi(basisFunctionIndex1D, xi[2]);
      }
      break;

    // 0-1+
    case ::Mesh::face_or_edge_t::edge0Minus1Plus:
      switch (dofIndex)
      {
      case 0:
      case 9:
      case 18:
        return xi1m*(2*xi1m - 1) * BasisFunctionType::phi(basisFunctionIndex1D, xi[2]);

      case 1:
      case 10:
      case 19:
        return 4*xi1m*(1 - xi1m - xi2) * BasisFunctionType::phi(basisFunctionIndex1D, xi[2]);

      case 2:
      case 11:
      case 20:
        return (xi1m + xi2 - 1)*(2*xi1m + 2*xi2 - 1) * BasisFunctionType::phi(basisFunctionIndex1D, xi[2]);

      case 3:
      case 12:
      case 21:
        return 0.0;

      case 4:
      case 13:
      case 22:
        return 0.0;

      case 5:
      case 14:
      case 23:
        return 4*xi2*(1 - xi1m - xi2) * BasisFunctionType::phi(basisFunctionIndex1D, xi[2]);

      case 6:
      case 15:
      case 24:
        return 4*xi1m*xi2 * BasisFunctionType::phi(basisFunctionIndex1D, xi[2]);

      case 7:
      case 16:
      case 25:
        return 0.0;

      case 8:
      case 17:
      case 26:
        return xi2*(2*xi2 - 1) * BasisFunctionType::phi(basisFunctionIndex1D, xi[2]);
      }
      break;

    // 0+1+
    case ::Mesh::face_or_edge_t::edge0Plus1Plus:
      switch (dofIndex)
      {
      case 0:
      case 9:
      case 18:
        return (xi1 + xi2 -1)*(2*xi1 + 2*xi2 - 1) * BasisFunctionType::phi(basisFunctionIndex1D, xi[2]);

      case 1:
      case 10:
      case 19:
        return 4*xi1*(1 - xi1 - xi2) * BasisFunctionType::phi(basisFunctionIndex1D, xi[2]);

      case 2:
      case 11:
      case 20:
        return xi1*(2*xi1 - 1) * BasisFunctionType::phi(basisFunctionIndex1D, xi[2]);

      case 3:
      case 12:
      case 21:
        return 4*xi2*(1 - xi1 - xi2) * BasisFunctionType::phi(basisFunctionIndex1D, xi[2]);

      case 4:
      case 13:
      case 22:
        return 0.0;

      case 5:
      case 14:
      case 23:
        return 0.0;

      case 6:
      case 15:
      case 24:
        return xi2*(2*xi2 - 1) * BasisFunctionType::phi(basisFunctionIndex1D, xi[2]);

      case 7:
      case 16:
      case 25:
        return 0.0;

      case 8:
      case 17:
      case 26:
        return 4*xi1*xi2 * BasisFunctionType::phi(basisFunctionIndex1D, xi[2]);
      }
      break;

    default:
      break;
    }
  }

  return FunctionSpaceFunction<::Mesh::StructuredDeformableOfDimension<3>, BasisFunction::LagrangeOfOrder<2>>::
    phiHexahedralMesh(dofIndex, xi);
}

template<typename Dummy>
double FunctionSpaceTriangleCorners<Mesh::StructuredDeformableOfDimension<3>, BasisFunction::LagrangeOfOrder<2>,Dummy>::
dphi_dxi(int dofIndex, int derivativeIdx, std::array<double,3> xi, element_no_t elementNoLocal) const
{
  using BasisFunctionType = BasisFunction::LagrangeOfOrder<2>;
  ::Mesh::face_or_edge_t edge;
  assert(derivativeIdx >= 0);
  assert(derivativeIdx < 3);

  // if the current element is a corner triangle
  if (elementNoLocal >= 0 && this->hasTriangleCorners_ && this->meshPartition_->elementIsAtCorner(elementNoLocal, edge))
  {
    int basisFunctionIndex1D = this->getBasisFunctionIndex1D(dofIndex, 2);
    const double xi1 = xi[0];
    const double xi2 = xi[1];
    const double xi1m = 1 - xi1;
    const double xi2m = 1 - xi2;

    // depending on the orientation of the triangle
    switch(edge)
    {
    // 0-1-
    case ::Mesh::face_or_edge_t::edge0Minus1Minus:
      switch (dofIndex)
      {
      case 0:
      case 9:
      case 18:
        if (derivativeIdx == 0)
        {
          return -4*xi2m * BasisFunctionType::phi(basisFunctionIndex1D, xi[2]);
        }
        else if (derivativeIdx == 1)
        {
          return -4*xi1m * BasisFunctionType::phi(basisFunctionIndex1D, xi[2]);
        }
        else if (derivativeIdx == 2)
        {
          return 4*xi1m*xi2m * BasisFunctionType::dphi_dxi(basisFunctionIndex1D, xi[2]);
        }

      case 1:
      case 10:
      case 19:
        return 0.0;

      case 2:
      case 11:
      case 20:
        if (derivativeIdx == 0)
        {
          return 0;
        }
        else if (derivativeIdx == 1)
        {
          return -(4*xi2m - 1) * BasisFunctionType::phi(basisFunctionIndex1D, xi[2]);  // 2*xi2m^2 - xi2m  ->  4*xi2m - 1
        }
        else if (derivativeIdx == 2)
        {
          return xi2m*(2*xi2m - 1) * BasisFunctionType::dphi_dxi(basisFunctionIndex1D, xi[2]);
        }

      case 3:
      case 12:
      case 21:
        return 0.0;

      case 4:
      case 13:
      case 22:
        return 0.0;

      case 5:
      case 14:
      case 23:
        if (derivativeIdx == 0)
        {
          return 4*xi2m * BasisFunctionType::phi(basisFunctionIndex1D, xi[2]);
        }
        else if (derivativeIdx == 1)
        {
          return -(4 - 4*xi1m - 8*xi2m) * BasisFunctionType::phi(basisFunctionIndex1D, xi[2]);   // 4*xi2m - 4*xi2m*xi1m - 4*xi2m*xi2m  ->  4 - 4*xi1m - 8*xi2m
        }
        else if (derivativeIdx == 2)
        {
          return 4*xi2m*(1 - xi1m - xi2m) * BasisFunctionType::dphi_dxi(basisFunctionIndex1D, xi[2]);
        }

      case 6:
      case 15:
      case 24:
        if (derivativeIdx == 0)
        {
          return -(4*xi1m - 1) * BasisFunctionType::phi(basisFunctionIndex1D, xi[2]);   // 2*xi1m^2 - xi1m  -> 4*xi1m - 1
        }
        else if (derivativeIdx == 1)
        {
          return 0;
        }
        else if (derivativeIdx == 2)
        {
          return xi1m*(2*xi1m - 1) * BasisFunctionType::dphi_dxi(basisFunctionIndex1D, xi[2]);
        }

      case 7:
      case 16:
      case 25:
        if (derivativeIdx == 0)
        {
          return -(4 - 8*xi1m - 4*xi2m) * BasisFunctionType::phi(basisFunctionIndex1D, xi[2]);  // 4*xi1m - 4*xi1m^2 - 4*xi1m*xi2m  -> 4 - 8*xi1m - 4*xi2m
        }
        else if (derivativeIdx == 1)
        {
          return 4*xi1m * BasisFunctionType::phi(basisFunctionIndex1D, xi[2]);  // 4*xi1m - 4*xi1m^2 - 4*xi1m*xi2m  -> -4*xi1m
        }
        else if (derivativeIdx == 2)
        {
          return 4*xi1m*(1 - xi1m - xi2m) * BasisFunctionType::dphi_dxi(basisFunctionIndex1D, xi[2]);
        }

      case 8:
      case 17:
      case 26:
        if (derivativeIdx == 0)
        {
          return -(4*xi1m + 4*xi2m - 3) * BasisFunctionType::phi(basisFunctionIndex1D, xi[2]);    // (2*xi1m^2 + 4*xi2m*xi1m + 2*xi2m^2 - 3*xi1m - 3*xi2m + 1) -> 4*xi1m + 4*xi2m - 3
        }
        else if (derivativeIdx == 1)
        {
          return -(4*xi2m + 4*xi1m - 3) * BasisFunctionType::phi(basisFunctionIndex1D, xi[2]);    // (2*xi1m^2 + 4*xi2m*xi1m + 2*xi2m^2 - 3*xi1m - 3*xi2m + 1) -> 4*xi1m + 4*xi2m - 3
        }
        else if (derivativeIdx == 2)
        {
          return (xi1m + xi2m - 1)*(2*xi1m + 2*xi2m - 1) * BasisFunctionType::dphi_dxi(basisFunctionIndex1D, xi[2]);
        }
      }
      break;

    // 0+1-
    case ::Mesh::face_or_edge_t::edge0Plus1Minus:
      switch (dofIndex)
      {
      case 0:
      case 9:
      case 18:
        if (derivativeIdx == 0)
        {
          return 0;
        }
        else if (derivativeIdx == 1)
        {
          return -(4*xi2m - 1) * BasisFunctionType::phi(basisFunctionIndex1D, xi[2]);      // 2*xi2m^2 - xi2m -> 4*xi2m - 1
        }
        else if (derivativeIdx == 2)
        {
          return xi2m*(2*xi2m - 1) * BasisFunctionType::dphi_dxi(basisFunctionIndex1D, xi[2]);
        }

      case 1:
      case 10:
      case 19:
        return 0.0;

      case 2:
      case 11:
      case 20:
        if (derivativeIdx == 0)
        {
          return 4*xi2m * BasisFunctionType::phi(basisFunctionIndex1D, xi[2]);
        }
        else if (derivativeIdx == 1)
        {
          return -4*xi1 * BasisFunctionType::phi(basisFunctionIndex1D, xi[2]);
        }
        else if (derivativeIdx == 2)
        {
          return 4*xi1*xi2m * BasisFunctionType::dphi_dxi(basisFunctionIndex1D, xi[2]);
        }

      case 3:
      case 12:
      case 21:
        if (derivativeIdx == 0)
        {
          return -4*xi2m * BasisFunctionType::phi(basisFunctionIndex1D, xi[2]);
        }
        else if (derivativeIdx == 1)
        {
          return -(4 - 4*xi1 - 8*xi2m) * BasisFunctionType::phi(basisFunctionIndex1D, xi[2]);    // 4*xi2m - 4*xi2m*xi1 - 4*xi2m^2   -> 4 - 4*xi1 - 8*xi2m
        }
        else if (derivativeIdx == 2)
        {
          return 4*xi2m*(1 - xi1 - xi2m) * BasisFunctionType::dphi_dxi(basisFunctionIndex1D, xi[2]);
        }

      case 4:
      case 13:
      case 22:
        return 0.0;

      case 5:
      case 14:
      case 23:
        return 0.0;

      case 6:
      case 15:
      case 24:
        if (derivativeIdx == 0)
        {
          return (4*xi1 + 4*xi2m - 3) * BasisFunctionType::phi(basisFunctionIndex1D, xi[2]);      // (xi1 + xi2m - 1)*(2*xi1 + 2*xi2m - 1)  -> 4*xi1 + 4*xi2m - 3
        }
        else if (derivativeIdx == 1)
        {
          return -(4*xi1 + 4*xi2m - 3) * BasisFunctionType::phi(basisFunctionIndex1D, xi[2]);
        }
        else if (derivativeIdx == 2)
        {
          return (xi1 + xi2m - 1)*(2*xi1 + 2*xi2m - 1) * BasisFunctionType::dphi_dxi(basisFunctionIndex1D, xi[2]);
        }

      case 7:
      case 16:
      case 25:
        if (derivativeIdx == 0)
        {
          return (4 - 8*xi1 - 4*xi2m) * BasisFunctionType::phi(basisFunctionIndex1D, xi[2]);     // 4*xi1 - 4*xi1^2 - 4*xi1*xi2m  -> 4 - 8*xi1 - 4*xi2m
        }
        else if (derivativeIdx == 1)
        {
          return 4*xi1 * BasisFunctionType::phi(basisFunctionIndex1D, xi[2]);     // 4*xi1 - 4*xi1^2 - 4*xi1*xi2m  -> -4*xi1
        }
        else if (derivativeIdx == 2)
        {
          return 4*xi1*(1 - xi1 - xi2m) * BasisFunctionType::dphi_dxi(basisFunctionIndex1D, xi[2]);
        }

      case 8:
      case 17:
      case 26:
        if (derivativeIdx == 0)
        {
          return (4*xi1 - 1) * BasisFunctionType::phi(basisFunctionIndex1D, xi[2]);    // 2*xi1^2 - xi1  -> 4*xi1 - 1
        }
        else if (derivativeIdx == 1)
        {
          return 0;
        }
        else if (derivativeIdx == 2)
        {
          return xi1*(2*xi1 - 1) * BasisFunctionType::dphi_dxi(basisFunctionIndex1D, xi[2]);
        }
      }
      break;

    // 0-1+
    case ::Mesh::face_or_edge_t::edge0Minus1Plus:
      switch (dofIndex)
      {
      case 0:
      case 9:
      case 18:
        if (derivativeIdx == 0)
        {
          return -(4*xi1m - 1) * BasisFunctionType::phi(basisFunctionIndex1D, xi[2]);    // (2*xi1m^2 - xi1m) -> 4*xi1m - 1
        }
        else if (derivativeIdx == 1)
        {
          return 0;
        }
        else if (derivativeIdx == 2)
        {
          return xi1m*(2*xi1m - 1) * BasisFunctionType::dphi_dxi(basisFunctionIndex1D, xi[2]);
        }

      case 1:
      case 10:
      case 19:
        if (derivativeIdx == 0)
        {
          return -(4 - 8*xi1m - 4*xi2) * BasisFunctionType::phi(basisFunctionIndex1D, xi[2]);    //  4*xi1m - 4*xi1m^2 - 4*xi1m*xi2 -> 4 - 8*xi1m - 4*xi2
        }
        else if (derivativeIdx == 1)
        {
          return -4*xi1m * BasisFunctionType::phi(basisFunctionIndex1D, xi[2]);   //  4*xi1m - 4*xi1m^2 - 4*xi1m*xi2 -> - 4*xi1m
        }
        else if (derivativeIdx == 2)
        {
          return 4*xi1m*(1 - xi1m - xi2) * BasisFunctionType::dphi_dxi(basisFunctionIndex1D, xi[2]);
        }

      case 2:
      case 11:
      case 20:
        if (derivativeIdx == 0)
        {
          return -(4*xi1m + 4*xi2 - 3) * BasisFunctionType::phi(basisFunctionIndex1D, xi[2]);    // 4*xi1m + 4*xi2 - 3
        }
        else if (derivativeIdx == 1)
        {
          return (4*xi1m + 4*xi2 - 3) * BasisFunctionType::phi(basisFunctionIndex1D, xi[2]);
        }
        else if (derivativeIdx == 2)
        {
          return (xi1m + xi2 - 1)*(2*xi1m + 2*xi2 - 1) * BasisFunctionType::dphi_dxi(basisFunctionIndex1D, xi[2]);
        }

      case 3:
      case 12:
      case 21:
        return 0.0;

      case 4:
      case 13:
      case 22:
        return 0.0;

      case 5:
      case 14:
      case 23:
        if (derivativeIdx == 0)
        {
          return 4*xi2 * BasisFunctionType::phi(basisFunctionIndex1D, xi[2]);   // 4*xi2 - 4*xi2*xi1m - 4*xi2^2 ->
        }
        else if (derivativeIdx == 1)
        {
          return (4 - 4*xi1m - 8*xi2) * BasisFunctionType::phi(basisFunctionIndex1D, xi[2]);    // 4*xi2 - 4*xi2*xi1m - 4*xi2^2 -> 4 - 4*xi1m - 8*xi2
        }
        else if (derivativeIdx == 2)
        {
          return 4*xi2*(1 - xi1m - xi2) * BasisFunctionType::dphi_dxi(basisFunctionIndex1D, xi[2]);
        }

      case 6:
      case 15:
      case 24:
        if (derivativeIdx == 0)
        {
          return -4*xi2 * BasisFunctionType::phi(basisFunctionIndex1D, xi[2]);
        }
        else if (derivativeIdx == 1)
        {
          return 4*xi1m * BasisFunctionType::phi(basisFunctionIndex1D, xi[2]);
        }
        else if (derivativeIdx == 2)
        {
          return 4*xi1m*xi2 * BasisFunctionType::dphi_dxi(basisFunctionIndex1D, xi[2]);
        }

      case 7:
      case 16:
      case 25:
        return 0.0;

      case 8:
      case 17:
      case 26:
        if (derivativeIdx == 0)
        {
          return 0;      // 2*xi2^2 - xi2 -> 4*xi2 - 1
        }
        else if (derivativeIdx == 1)
        {
          return (4*xi2 - 1) * BasisFunctionType::phi(basisFunctionIndex1D, xi[2]);
        }
        else if (derivativeIdx == 2)
        {
          return xi2*(2*xi2 - 1) * BasisFunctionType::dphi_dxi(basisFunctionIndex1D, xi[2]);
        }
      }
      break;

    // 0+1+
    case ::Mesh::face_or_edge_t::edge0Plus1Plus:
      switch (dofIndex)
      {
      case 0:
      case 9:
      case 18:
        if (derivativeIdx == 0)
        {
          return (4*xi1 + 4*xi2 - 3) * BasisFunctionType::phi(basisFunctionIndex1D, xi[2]);   // (xi1 + xi2 -1)*(2*xi1 + 2*xi2 - 1)
        }
        else if (derivativeIdx == 1)
        {
          return (4*xi1 + 4*xi2 - 3) * BasisFunctionType::phi(basisFunctionIndex1D, xi[2]);
        }
        else if (derivativeIdx == 2)
        {
          return (xi1 + xi2 -1)*(2*xi1 + 2*xi2 - 1) * BasisFunctionType::dphi_dxi(basisFunctionIndex1D, xi[2]);
        }

      case 1:
      case 10:
      case 19:
        if (derivativeIdx == 0)
        {
          return (4 - 8*xi1 - 4*xi2) * BasisFunctionType::phi(basisFunctionIndex1D, xi[2]);    // 4*xi1 - 4*xi1^2 - 4*xi1*xi2 -> 4 - 8*xi1 - 4*xi2
        }
        else if (derivativeIdx == 1)
        {
          return -4*xi1 * BasisFunctionType::phi(basisFunctionIndex1D, xi[2]);    // 4*xi1 - 4*xi1^2 - 4*xi1*xi2 -> -4*xi1
        }
        else if (derivativeIdx == 2)
        {
          return 4*xi1*(1 - xi1 - xi2) * BasisFunctionType::dphi_dxi(basisFunctionIndex1D, xi[2]);
        }

      case 2:
      case 11:
      case 20:
        if (derivativeIdx == 0)
        {
          return (4*xi1 - 1) * BasisFunctionType::phi(basisFunctionIndex1D, xi[2]);    // (2*xi1^2 - xi1) -> 4*xi1 - 1
        }
        else if (derivativeIdx == 1)
        {
          return 0;
        }
        else if (derivativeIdx == 2)
        {
          return xi1*(2*xi1 - 1) * BasisFunctionType::dphi_dxi(basisFunctionIndex1D, xi[2]);
        }

      case 3:
      case 12:
      case 21:
        if (derivativeIdx == 0)
        {
          return -4*xi2 * BasisFunctionType::phi(basisFunctionIndex1D, xi[2]);    // 4*xi2 - 4*xi2*xi1 - 4*xi2^2 -> -4*xi2
        }
        else if (derivativeIdx == 1)
        {
          return (4 - 4*xi1 - 8*xi2) * BasisFunctionType::phi(basisFunctionIndex1D, xi[2]);     // 4*xi2 - 4*xi2*xi1 - 4*xi2^2 -> 4 - 4*xi1 - 8*xi2
        }
        else if (derivativeIdx == 2)
        {
          return 4*xi2*(1 - xi1 - xi2) * BasisFunctionType::dphi_dxi(basisFunctionIndex1D, xi[2]);
        }

      case 4:
      case 13:
      case 22:
        return 0.0;

      case 5:
      case 14:
      case 23:
        return 0.0;

      case 6:
      case 15:
      case 24:
        if (derivativeIdx == 0)
        {
          return 0;
        }
        else if (derivativeIdx == 1)
        {
          return (4*xi2 - 1) * BasisFunctionType::phi(basisFunctionIndex1D, xi[2]);    // 2*xi2^2 - xi2  -> 4*xi2 - 1
        }
        else if (derivativeIdx == 2)
        {
          return xi2*(2*xi2 - 1) * BasisFunctionType::dphi_dxi(basisFunctionIndex1D, xi[2]);
        }

      case 7:
      case 16:
      case 25:
        return 0.0;

      case 8:
      case 17:
      case 26:
        if (derivativeIdx == 0)
        {
          return 4*xi2 * BasisFunctionType::phi(basisFunctionIndex1D, xi[2]);
        }
        else if (derivativeIdx == 1)
        {
          return 4*xi1 * BasisFunctionType::phi(basisFunctionIndex1D, xi[2]);
        }
        else if (derivativeIdx == 2)
        {
          return 4*xi1*xi2 * BasisFunctionType::dphi_dxi(basisFunctionIndex1D, xi[2]);
        }
      }
      break;

    default:
      break;
    }
  }

  return FunctionSpaceFunction<::Mesh::StructuredDeformableOfDimension<3>, BasisFunction::LagrangeOfOrder<2>>::
    dphi_dxiHexahedralMesh(dofIndex, derivativeIdx, xi);
}

//! set the dependent dofs in the given field variable by interpolating the independent dofs of the triangle basis
template<typename Dummy>
void FunctionSpaceTriangleCorners<Mesh::StructuredDeformableOfDimension<3>, BasisFunction::LagrangeOfOrder<2>,Dummy>::
interpolateNonDofValuesInFieldVariable(std::shared_ptr<FieldVariable::FieldVariableBaseFunctionSpace<FunctionSpace<::Mesh::StructuredDeformableOfDimension<3>,BasisFunction::LagrangeOfOrder<2>>>> fieldVariable, int componentNo) const
{
  if (!fieldVariable || !this->hasTriangleCorners_)
    return;

  const int nDofsPerElement = 27;
  ::Mesh::face_or_edge_t edge;

  // get local values
  std::vector<double> valuesLocal;
  fieldVariable->getValuesWithoutGhosts(componentNo, valuesLocal);

  std::vector<double> valuesLocalNew = valuesLocal;

  std::shared_ptr<FunctionSpacePointInElement<::Mesh::StructuredDeformableOfDimension<3>,BasisFunction::LagrangeOfOrder<2>>> functionSpace
    = fieldVariable->functionSpace();
  dof_no_t nDofsLocalWithoutGhosts = functionSpace->nDofsLocalWithoutGhosts();

  // iterate over elements
  for (element_no_t elementNoLocal = 0; elementNoLocal < functionSpace->nElementsLocal(); elementNoLocal++)
  {
    // if the element is a triangle at the corner
    if (this->meshPartition_->elementIsAtCorner(elementNoLocal, edge))
    {
      std::array<dof_no_t,nDofsPerElement> dofNosLocal = functionSpace->getElementDofNosLocal(elementNoLocal);

      if (VLOG_IS_ON(1))
        VLOG(1) << "quadratic element " << elementNoLocal << ", dofNosLocal: " << dofNosLocal;

      int maxZLevel = 3;
      if (dofNosLocal[18] >= nDofsLocalWithoutGhosts)
        maxZLevel = 2;

      // iterate over z levels
      for (int i = 0; i < maxZLevel; i++)
      {
        std::stringstream message;
        if (VLOG_IS_ON(1))
        {
          message << "el " << elementNoLocal << " " << ::Mesh::getString(edge) << " "
            << "(" << valuesLocal[dofNosLocal[i*9 + 1]] << "," << valuesLocal[dofNosLocal[i*9 + 3]]
            << "," << valuesLocal[dofNosLocal[i*9 + 4]] << "," << valuesLocal[dofNosLocal[i*9 + 5]]
            << "," << valuesLocal[dofNosLocal[i*9 + 7]] << ")";
        }

        // depending on the orientation of the triangle
        switch(edge)
        {
        // 0-1-
        case ::Mesh::face_or_edge_t::edge0Minus1Minus:

          //valuesLocalNew[dofNosLocal[i*9 + 0]] = 0.5*(valuesLocal[dofNosLocal[i*9 + 2]] + valuesLocal[dofNosLocal[i*9 + 6]]);

          valuesLocalNew[dofNosLocal[i*9 + 1]] = 0.5*(valuesLocal[dofNosLocal[i*9 + 0]] + valuesLocal[dofNosLocal[i*9 + 2]]);
          valuesLocalNew[dofNosLocal[i*9 + 3]] = 0.5*(valuesLocal[dofNosLocal[i*9 + 0]] + valuesLocal[dofNosLocal[i*9 + 6]]);
          valuesLocalNew[dofNosLocal[i*9 + 4]] = 0.5*(valuesLocal[dofNosLocal[i*9 + 0]] + valuesLocal[dofNosLocal[i*9 + 8]]);
          break;

        // 0+1-
        case ::Mesh::face_or_edge_t::edge0Plus1Minus:
          //valuesLocalNew[dofNosLocal[i*9 + 2]] = 0.5*(valuesLocal[dofNosLocal[i*9 + 0]] + valuesLocal[dofNosLocal[i*9 + 8]]);

          valuesLocalNew[dofNosLocal[i*9 + 1]] = 0.5*(valuesLocal[dofNosLocal[i*9 + 0]] + valuesLocal[dofNosLocal[i*9 + 2]]);
          valuesLocalNew[dofNosLocal[i*9 + 4]] = 0.5*(valuesLocal[dofNosLocal[i*9 + 2]] + valuesLocal[dofNosLocal[i*9 + 6]]);
          valuesLocalNew[dofNosLocal[i*9 + 5]] = 0.5*(valuesLocal[dofNosLocal[i*9 + 2]] + valuesLocal[dofNosLocal[i*9 + 8]]);
          break;

        // 0-1+
        case ::Mesh::face_or_edge_t::edge0Minus1Plus:
          //valuesLocalNew[dofNosLocal[i*9 + 6]] = 0.5*(valuesLocal[dofNosLocal[i*9 + 0]] + valuesLocal[dofNosLocal[i*9 + 8]]);

          valuesLocalNew[dofNosLocal[i*9 + 3]] = 0.5*(valuesLocal[dofNosLocal[i*9 + 0]] + valuesLocal[dofNosLocal[i*9 + 6]]);
          valuesLocalNew[dofNosLocal[i*9 + 4]] = 0.5*(valuesLocal[dofNosLocal[i*9 + 2]] + valuesLocal[dofNosLocal[i*9 + 6]]);
          valuesLocalNew[dofNosLocal[i*9 + 7]] = 0.5*(valuesLocal[dofNosLocal[i*9 + 6]] + valuesLocal[dofNosLocal[i*9 + 8]]);
          break;

        // 0+1+
        case ::Mesh::face_or_edge_t::edge0Plus1Plus:
          //valuesLocalNew[dofNosLocal[i*9 + 8]] = 0.5*(valuesLocal[dofNosLocal[i*9 + 2]] + valuesLocal[dofNosLocal[i*9 + 6]]);

          valuesLocalNew[dofNosLocal[i*9 + 4]] = 0.5*(valuesLocal[dofNosLocal[i*9 + 0]] + valuesLocal[dofNosLocal[i*9 + 8]]);
          valuesLocalNew[dofNosLocal[i*9 + 5]] = 0.5*(valuesLocal[dofNosLocal[i*9 + 2]] + valuesLocal[dofNosLocal[i*9 + 8]]);
          valuesLocalNew[dofNosLocal[i*9 + 7]] = 0.5*(valuesLocal[dofNosLocal[i*9 + 6]] + valuesLocal[dofNosLocal[i*9 + 8]]);
          break;

        default:
          break;
        }

        if (VLOG_IS_ON(1))
        {
          message << " -> "
            << "(" << valuesLocalNew[dofNosLocal[i*9 + 1]] << "," << valuesLocalNew[dofNosLocal[i*9 + 3]]
            << "," << valuesLocalNew[dofNosLocal[i*9 + 4]] << "," << valuesLocalNew[dofNosLocal[i*9 + 5]]
            << "," << valuesLocalNew[dofNosLocal[i*9 + 7]] << ")";
          VLOG(1) << message.str();
        }
      }
    }
  }

  fieldVariable->setValuesWithoutGhosts(componentNo, valuesLocalNew);
}

template<typename Dummy>
template<int nComponents>
void FunctionSpaceTriangleCorners<Mesh::StructuredDeformableOfDimension<3>, BasisFunction::LagrangeOfOrder<2>,Dummy>::
interpolateNonDofValuesInFieldVariable(std::shared_ptr<FieldVariable::FieldVariable<FunctionSpace<::Mesh::StructuredDeformableOfDimension<3>,BasisFunction::LagrangeOfOrder<2>>,nComponents>> fieldVariable) const
{
  for (int componentNo = 0; componentNo < nComponents; componentNo++)
  {
    interpolateNonDofValuesInFieldVariable(fieldVariable, componentNo);
  }
}


} // namespace
