#include "quadrature/triangular_prism.h"

#include "easylogging++.h"
#include <set>

namespace Quadrature
{

// stubs for D != 3
template<int D, typename GaussQuadrature>
std::array<std::array<double,D>,TriangularPrismBase<GaussQuadrature>::numberEvaluations()> TriangularPrism<D,GaussQuadrature>::
samplingPoints()
{
  return std::array<std::array<double,D>,TriangularPrismBase<GaussQuadrature>::numberEvaluations()>{};
}

template<int D, typename GaussQuadrature>
const std::array<double,TriangularPrismBase<GaussQuadrature>::numberEvaluations()> TriangularPrism<D,GaussQuadrature>::
quadratureWeights()
{
  return std::array<double,TriangularPrismBase<GaussQuadrature>::numberEvaluations()>{};
}

template<int D, typename GaussQuadrature>
template<typename ValueType, long unsigned int nEntriesEvaluationArray>
ValueType TriangularPrism<D,GaussQuadrature>::
computeIntegral(const typename std::array<ValueType,nEntriesEvaluationArray> &evaluations)
{
  ValueType result{};
  return result;
}

template<int D, typename GaussQuadrature>
template<typename ValueType, long unsigned int nEntriesEvaluationArray>
void TriangularPrism<D,GaussQuadrature>::
adjustEntriesforPrism(typename std::array<ValueType,nEntriesEvaluationArray> &evaluations,
                      Mesh::face_or_edge_t edge, bool usesDofPairs, bool isQuadraticElement0, int nEntriesPerDof0,
                      bool isQuadraticElement1, int nEntriesPerDof1)
{
}

// actual implementations

// 3D tensor product integration
template<typename GaussQuadrature>
template<typename ValueType, long unsigned int nEntriesEvaluationArray>
ValueType TriangularPrism<3,GaussQuadrature>::
computeIntegral(const typename std::array<ValueType,nEntriesEvaluationArray> &evaluations)
{
  const int numberEvaluations = TriangularPrism<3,GaussQuadrature>::numberEvaluations();
  assert(numberEvaluations <= nEntriesEvaluationArray);
  const std::array<double,numberEvaluations> weights = TriangularPrism<3,GaussQuadrature>::quadratureWeights();

  ValueType result{};
  for (int i = 0; i < numberEvaluations; i++)
  {
    result += weights[i] * evaluations[i];
  }
  //LOG(DEBUG) << "quadrature " << weights << "*" << evaluations << " -> " << result;

  return result;
}

// reorder entries
template<typename GaussQuadrature>
template<typename ValueType, long unsigned int nEntriesEvaluationArray>
void TriangularPrism<3,GaussQuadrature>::
adjustEntriesforPrism(typename std::array<ValueType,nEntriesEvaluationArray> &evaluations,
                      Mesh::face_or_edge_t edge, bool usesDofPairs, bool isQuadraticElement0, int nEntriesPerDof0,
                      bool isQuadraticElement1, int nEntriesPerDof1)
{
  //typename std::array<ValueType,TriangularPrism<3,GaussQuadrature>::numberEvaluations()> &evaluations
  int interpolatedDofLinear = -1;
  std::set<int> interpolatedDofsQuadratic;

  if (edge == Mesh::face_or_edge_t::edge0Minus1Minus)
  {
    interpolatedDofLinear = 0;
    interpolatedDofsQuadratic = std::set<int>{1,3,4};
  }
  else if (edge == Mesh::face_or_edge_t::edge0Plus1Minus)
  {
    interpolatedDofLinear = 1;
    interpolatedDofsQuadratic = std::set<int>{1,4,5};
  }
  else if (edge == Mesh::face_or_edge_t::edge0Minus1Plus)
  {
    interpolatedDofLinear = 2;
    interpolatedDofsQuadratic = std::set<int>{3,4,7};
  }
  else if (edge == Mesh::face_or_edge_t::edge0Plus1Plus)
  {
    interpolatedDofLinear = 3;
    interpolatedDofsQuadratic = std::set<int>{4,5,7};
  }

  int nDofsXY0 = 4;
  int nDofsZ0 = 2;
  if (isQuadraticElement0)
  {
    nDofsXY0 = 9;
    nDofsZ0 = 3;
  }

  int nDofsXY1 = 4;
  int nDofsZ1 = 2;
  if (isQuadraticElement1)
  {
    nDofsXY1 = 9;
    nDofsZ1 = 3;
  }

  if (!usesDofPairs)
  {
    nDofsXY1 = 1;
    nDofsZ1 = 1;
  }

  int nDofs0 = nDofsXY0 * nDofsZ0;
  int nDofs1 = nDofsXY1 * nDofsZ1;

  VLOG(1) << "edge " << Mesh::getString(edge) << ", isQuadratic: " << isQuadraticElement0 << "," << isQuadraticElement1
    << ", nEntriesPerDof: " << nEntriesPerDof0 << "," << nEntriesPerDof1 << ", nDofs: " << nDofs0 << "," << nDofs1;

  // iterate over pairs of dofs in the hex element
  // iterate over first dof
  for (int dofNoZ1 = 0; dofNoZ1 < nDofsZ1; dofNoZ1++)
  {
    for (int dofNoXY1 = 0; dofNoXY1 < nDofsXY1; dofNoXY1++)
    {
      int dofNo1 = dofNoZ1 * nDofsXY1 + dofNoXY1;

      bool dofIsInterpolated1 = false;

      // if the first dof is an interpolated dof in the prism
      if (usesDofPairs &&
          ((isQuadraticElement1 && interpolatedDofsQuadratic.find(dofNoXY1) != interpolatedDofsQuadratic.end())
          || (!isQuadraticElement1 && dofNoXY1 == interpolatedDofLinear)))
      {
        dofIsInterpolated1 = true;
      }

      // iterate over second dof
      for (int dofNoZ0 = 0; dofNoZ0 < nDofsZ0; dofNoZ0++)
      {
        for (int dofNoXY0 = 0; dofNoXY0 < nDofsXY0; dofNoXY0++)
        {
          int dofNo0 = dofNoZ0 * nDofsXY0 + dofNoXY0;

          bool dofIsInterpolated0 = false;

          // if the second dof is an interpolated dof in the prism
          if ((isQuadraticElement0 && interpolatedDofsQuadratic.find(dofNoXY0) != interpolatedDofsQuadratic.end())
            || (!isQuadraticElement0 && dofNoXY0 == interpolatedDofLinear))
          {
            dofIsInterpolated0 = true;
          }

          // iterate over the entries per dofs, for both dofs
          for (int entryNoPerDof1 = 0; entryNoPerDof1 < nEntriesPerDof1; entryNoPerDof1++)
          {
            for (int entryNoPerDof0 = 0; entryNoPerDof0 < nEntriesPerDof0; entryNoPerDof0++)
            {
              // calculate index in evaluations array
              int index = (dofNo1*nEntriesPerDof1 + entryNoPerDof1)*nDofs0*nEntriesPerDof0 + dofNo0*nEntriesPerDof0 + entryNoPerDof0;
              assert(index < evaluations.size());

              // set entry to 0
              if (dofIsInterpolated1 || dofIsInterpolated0)
              {
                if (dofNo0 == dofNo1 && entryNoPerDof0 == entryNoPerDof1)
                {
                  evaluations[index] = 1;
                }
                else
                {
                  evaluations[index] = 0;
                }
              }
            }
          }
        }
      }
    }
  }
}

}  // namespace
