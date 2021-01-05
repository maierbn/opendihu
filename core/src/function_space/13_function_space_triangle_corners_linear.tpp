#include "function_space/13_function_space_triangle_corners.h"

#include "partition/mesh_partition/01_mesh_partition_structured.h"

namespace FunctionSpace
{
//--------------------------
// linear

template<typename Dummy>
double FunctionSpaceTriangleCorners<Mesh::StructuredDeformableOfDimension<3>,BasisFunction::LagrangeOfOrder<1>,Dummy>::
phi(int dofIndex, std::array<double,3> xi, element_no_t elementNoLocal) const
{
  using BasisFunctionType = BasisFunction::LagrangeOfOrder<1>;
  ::Mesh::face_or_edge_t edge;

  // if the current element is a corner triangle
  if (elementNoLocal >= 0 && this->hasTriangleCorners_ && this->meshPartition_->elementIsAtCorner(elementNoLocal, edge))
  {
    int basisFunctionIndex1D = this->getBasisFunctionIndex1D(dofIndex, 2);
    const double xi1 = xi[0];
    const double xi2 = xi[1];

    // depending on the orientation of the triangle
    switch(edge)
    {
    // 0-1-
    case ::Mesh::face_or_edge_t::edge0Minus1Minus:
      switch (dofIndex)
      {
      case 0:
      case 4:
        return 0.0;

      case 1:
      case 5:
        return (1-xi2) * BasisFunctionType::phi(basisFunctionIndex1D, xi[2]);

      case 2:
      case 6:
        return (1-xi1) * BasisFunctionType::phi(basisFunctionIndex1D, xi[2]);

      case 3:
      case 7:
        return (-1+xi1+xi2) * BasisFunctionType::phi(basisFunctionIndex1D, xi[2]);  // (1-(1-xi1)-(1-xi2)) = (1-1+xi1-1+xi2) = -1+xi1+xi2
      }
      break;

    // 0+1-
    case ::Mesh::face_or_edge_t::edge0Plus1Minus:
      switch (dofIndex)
      {
      case 0:
      case 4:
        return (1-xi2) * BasisFunctionType::phi(basisFunctionIndex1D, xi[2]);

      case 1:
      case 5:
        return 0.0;

      case 2:
      case 6:
        return (-xi1+xi2) * BasisFunctionType::phi(basisFunctionIndex1D, xi[2]);  // (1-xi1-(1-xi2)) = (-xi1+xi2)

      case 3:
      case 7:
        return xi1 * BasisFunctionType::phi(basisFunctionIndex1D, xi[2]);
      }
      break;

    // 0-1+
    case ::Mesh::face_or_edge_t::edge0Minus1Plus:
      switch (dofIndex)
      {
      case 0:
      case 4:
        return (1-xi1) * BasisFunctionType::phi(basisFunctionIndex1D, xi[2]);

      case 1:
      case 5:
        return (xi1-xi2) * BasisFunctionType::phi(basisFunctionIndex1D, xi[2]);   // (1-(1-xi1)-xi2) = (xi1-xi2)

      case 2:
      case 6:
        return 0.0;

      case 3:
      case 7:
        return xi2 * BasisFunctionType::phi(basisFunctionIndex1D, xi[2]);
      }
      break;

    // 0+1+
    case ::Mesh::face_or_edge_t::edge0Plus1Plus:
      switch (dofIndex)
      {
      case 0:
      case 4:
        return (1-xi1-xi2) * BasisFunctionType::phi(basisFunctionIndex1D, xi[2]);

      case 1:
      case 5:
        return xi1 * BasisFunctionType::phi(basisFunctionIndex1D, xi[2]);

      case 2:
      case 6:
        return xi2 * BasisFunctionType::phi(basisFunctionIndex1D, xi[2]);

      case 3:
      case 7:
        return 0.0;
      }
      break;

    default:
      break;
    }
  }

  // return normal function value
  return FunctionSpaceFunction<::Mesh::StructuredDeformableOfDimension<3>, BasisFunction::LagrangeOfOrder<1>>::
    phiHexahedralMesh(dofIndex, xi);
}

template<typename Dummy>
double FunctionSpaceTriangleCorners<Mesh::StructuredDeformableOfDimension<3>,BasisFunction::LagrangeOfOrder<1>,Dummy>::
dphi_dxi(int dofIndex, int derivativeIdx, std::array<double,3> xi, element_no_t elementNoLocal) const
{
  using BasisFunctionType = BasisFunction::LagrangeOfOrder<1>;
  ::Mesh::face_or_edge_t edge;
  assert(derivativeIdx >= 0);
  assert(derivativeIdx < 3);

  // if the current element is a corner triangle
  if (elementNoLocal >= 0 && this->hasTriangleCorners_ && this->meshPartition_->elementIsAtCorner(elementNoLocal, edge))
  {
    int basisFunctionIndex1D = this->getBasisFunctionIndex1D(dofIndex, 2);
    const double xi1 = xi[0];
    const double xi2 = xi[1];

    // depending on the orientation of the triangle
    switch(edge)
    {
    // 0-1-
    case ::Mesh::face_or_edge_t::edge0Minus1Minus:
      switch (dofIndex)
      {
      case 0:
      case 4:
        return 0.0;

      case 1:
      case 5:
        if (derivativeIdx == 0)
        {
          return 0;
        }
        else if (derivativeIdx == 1)
        {
          return -BasisFunctionType::phi(basisFunctionIndex1D, xi[2]);      // (1-xi3) or xi3
        }
        else if (derivativeIdx == 2)
        {
          return (1-xi2) * BasisFunctionType::dphi_dxi(basisFunctionIndex1D, xi[2]);    // -1 or 1
        }

      case 2:
      case 6:
        if (derivativeIdx == 0)
        {
          return -BasisFunctionType::phi(basisFunctionIndex1D, xi[2]);
        }
        else if (derivativeIdx == 1)
        {
          return 0;
        }
        else if (derivativeIdx == 2)
        {
          return (1-xi1) * BasisFunctionType::dphi_dxi(basisFunctionIndex1D, xi[2]);
        }

      case 3:
      case 7:
        if (derivativeIdx == 0)
        {
          return BasisFunctionType::phi(basisFunctionIndex1D, xi[2]);
        }
        else if (derivativeIdx == 1)
        {
          return BasisFunctionType::phi(basisFunctionIndex1D, xi[2]);
        }
        else if (derivativeIdx == 2)
        {
          return (-1+xi1+xi2) * BasisFunctionType::dphi_dxi(basisFunctionIndex1D, xi[2]);
        }
      }
      break;

    // 0+1-
    case ::Mesh::face_or_edge_t::edge0Plus1Minus:
      switch (dofIndex)
      {
      case 0:
      case 4:
        if (derivativeIdx == 0)
        {
          return 0;
        }
        else if (derivativeIdx == 1)
        {
          return -BasisFunctionType::phi(basisFunctionIndex1D, xi[2]);
        }
        else if (derivativeIdx == 2)
        {
          return (1-xi2) * BasisFunctionType::dphi_dxi(basisFunctionIndex1D, xi[2]);
        }

      case 1:
      case 5:
        return 0.0;

      case 2:
      case 6:
        if (derivativeIdx == 0)
        {
          return -BasisFunctionType::phi(basisFunctionIndex1D, xi[2]);
        }
        else if (derivativeIdx == 1)
        {
          return BasisFunctionType::phi(basisFunctionIndex1D, xi[2]);
        }
        else if (derivativeIdx == 2)
        {
          return (-xi1+xi2) * BasisFunctionType::dphi_dxi(basisFunctionIndex1D, xi[2]);
        }

      case 3:
      case 7:
        if (derivativeIdx == 0)
        {
          return BasisFunctionType::phi(basisFunctionIndex1D, xi[2]);
        }
        else if (derivativeIdx == 1)
        {
          return 0;
        }
        else if (derivativeIdx == 2)
        {
          return xi1 * BasisFunctionType::dphi_dxi(basisFunctionIndex1D, xi[2]);
        }
      }
      break;

    // 0-1+
    case ::Mesh::face_or_edge_t::edge0Minus1Plus:
      switch (dofIndex)
      {
      case 0:
      case 4:
        if (derivativeIdx == 0)
        {
          return -BasisFunctionType::phi(basisFunctionIndex1D, xi[2]);
        }
        else if (derivativeIdx == 1)
        {
          return 0;
        }
        else if (derivativeIdx == 2)
        {
          return (1-xi1) * BasisFunctionType::dphi_dxi(basisFunctionIndex1D, xi[2]);
        }

      case 1:
      case 5:
        if (derivativeIdx == 0)
        {
          return BasisFunctionType::phi(basisFunctionIndex1D, xi[2]);
        }
        else if (derivativeIdx == 1)
        {
          return -BasisFunctionType::phi(basisFunctionIndex1D, xi[2]);
        }
        else if (derivativeIdx == 2)
        {
          return (xi1-xi2) * BasisFunctionType::dphi_dxi(basisFunctionIndex1D, xi[2]);
        }

      case 2:
      case 6:
        return 0.0;

      case 3:
      case 7:
        if (derivativeIdx == 0)
        {
          return 0;
        }
        else if (derivativeIdx == 1)
        {
          return BasisFunctionType::phi(basisFunctionIndex1D, xi[2]);
        }
        else if (derivativeIdx == 2)
        {
          return xi2 * BasisFunctionType::dphi_dxi(basisFunctionIndex1D, xi[2]);
        }
      }
      break;

    // 0+1+
    case ::Mesh::face_or_edge_t::edge0Plus1Plus:
      switch (dofIndex)
      {
      case 0:
      case 4:
        if (derivativeIdx == 0)
        {
          return -BasisFunctionType::phi(basisFunctionIndex1D, xi[2]);
        }
        else if (derivativeIdx == 1)
        {
          return -BasisFunctionType::phi(basisFunctionIndex1D, xi[2]);
        }
        else if (derivativeIdx == 2)
        {
          return (1-xi1-xi2) * BasisFunctionType::dphi_dxi(basisFunctionIndex1D, xi[2]);
        }
        break;

      case 1:
      case 5:
        if (derivativeIdx == 0)
        {
          return BasisFunctionType::phi(basisFunctionIndex1D, xi[2]);
        }
        else if (derivativeIdx == 1)
        {
          return 0;
        }
        else if (derivativeIdx == 2)
        {
          return xi1 * BasisFunctionType::dphi_dxi(basisFunctionIndex1D, xi[2]);
        }
        break;

      case 2:
      case 6:
        if (derivativeIdx == 0)
        {
          return 0;
        }
        else if (derivativeIdx == 1)
        {
          return BasisFunctionType::phi(basisFunctionIndex1D, xi[2]);
        }
        else if (derivativeIdx == 2)
        {
          return xi2 * BasisFunctionType::dphi_dxi(basisFunctionIndex1D, xi[2]);
        }

      case 3:
      case 7:
        return 0.0;
      }
      break;

    default:
      break;
    }
  }

  return FunctionSpaceFunction<::Mesh::StructuredDeformableOfDimension<3>, BasisFunction::LagrangeOfOrder<1>>::
    dphi_dxiHexahedralMesh(dofIndex, derivativeIdx, xi);
}

//! set the dependent dofs in the given field variable by interpolating the independent dofs of the triangle basis
template<typename Dummy>
void FunctionSpaceTriangleCorners<Mesh::StructuredDeformableOfDimension<3>, BasisFunction::LagrangeOfOrder<1>,Dummy>::
interpolateNonDofValuesInFieldVariable(std::shared_ptr<FieldVariable::FieldVariableBaseFunctionSpace<FunctionSpace<::Mesh::StructuredDeformableOfDimension<3>,BasisFunction::LagrangeOfOrder<1>>>> fieldVariable, int componentNo) const
{
  if (!fieldVariable || !this->hasTriangleCorners_)
    return;

  const int nDofsPerElement = 8;
  ::Mesh::face_or_edge_t edge;

  // get local values
  std::vector<double> valuesLocal;
  fieldVariable->getValuesWithoutGhosts(componentNo, valuesLocal);


  std::shared_ptr<FunctionSpacePointInElement<::Mesh::StructuredDeformableOfDimension<3>,BasisFunction::LagrangeOfOrder<1>>> functionSpace
    = fieldVariable->functionSpace();
  dof_no_t nDofsLocalWithoutGhosts = functionSpace->nDofsLocalWithoutGhosts();

  std::vector<double> valuesLocalNew = valuesLocal;
  //LOG(DEBUG) << "size valuesLocal: " << valuesLocal.size();

  // iterate over elements
  for (element_no_t elementNoLocal = 0; elementNoLocal < functionSpace->nElementsLocal(); elementNoLocal++)
  {
    // if the element is a triangle at the corner
    if (this->meshPartition_->elementIsAtCorner(elementNoLocal, edge))
    {
      std::array<dof_no_t,nDofsPerElement> dofNosLocal = functionSpace->getElementDofNosLocal(elementNoLocal);

      //LOG(DEBUG) << "linear element " << elementNoLocal << ", dofNosLocal: " << dofNosLocal;

      // depending on the orientation of the triangle
      switch(edge)
      {
      // 0-1-
      case ::Mesh::face_or_edge_t::edge0Minus1Minus:
        valuesLocalNew[dofNosLocal[0]] = 0.5*(valuesLocal[dofNosLocal[1]] + valuesLocal[dofNosLocal[2]]);
        if (dofNosLocal[4] < nDofsLocalWithoutGhosts)
          valuesLocalNew[dofNosLocal[4]] = 0.5*(valuesLocal[dofNosLocal[5]] + valuesLocal[dofNosLocal[6]]);
        break;

      // 0+1-
      case ::Mesh::face_or_edge_t::edge0Plus1Minus:
        valuesLocalNew[dofNosLocal[1]] = 0.5*(valuesLocal[dofNosLocal[0]] + valuesLocal[dofNosLocal[3]]);
        if (dofNosLocal[4] < nDofsLocalWithoutGhosts)
          valuesLocalNew[dofNosLocal[5]] = 0.5*(valuesLocal[dofNosLocal[4]] + valuesLocal[dofNosLocal[7]]);
        break;

      // 0-1+
      case ::Mesh::face_or_edge_t::edge0Minus1Plus:
        valuesLocalNew[dofNosLocal[2]] = 0.5*(valuesLocal[dofNosLocal[0]] + valuesLocal[dofNosLocal[3]]);
        if (dofNosLocal[4] < nDofsLocalWithoutGhosts)
          valuesLocalNew[dofNosLocal[6]] = 0.5*(valuesLocal[dofNosLocal[4]] + valuesLocal[dofNosLocal[7]]);
        break;

      // 0+1+
      case ::Mesh::face_or_edge_t::edge0Plus1Plus:
        valuesLocalNew[dofNosLocal[3]] = 0.5*(valuesLocal[dofNosLocal[1]] + valuesLocal[dofNosLocal[2]]);
        if (dofNosLocal[4] < nDofsLocalWithoutGhosts)
          valuesLocalNew[dofNosLocal[7]] = 0.5*(valuesLocal[dofNosLocal[5]] + valuesLocal[dofNosLocal[6]]);
        break;

      default:
        break;
      }
    }
  }

  fieldVariable->setValuesWithoutGhosts(componentNo, valuesLocalNew);
}

template<typename Dummy>
template<int nComponents>
void FunctionSpaceTriangleCorners<Mesh::StructuredDeformableOfDimension<3>, BasisFunction::LagrangeOfOrder<1>,Dummy>::
interpolateNonDofValuesInFieldVariable(std::shared_ptr<FieldVariable::FieldVariable<FunctionSpace<::Mesh::StructuredDeformableOfDimension<3>,BasisFunction::LagrangeOfOrder<1>>,nComponents>> fieldVariable) const
{
  for (int componentNo = 0; componentNo < nComponents; componentNo++)
  {
    interpolateNonDofValuesInFieldVariable(fieldVariable, componentNo);
  }
}

} // namespace
