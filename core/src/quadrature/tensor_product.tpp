#include "quadrature/tensor_product.h"

#include "utility/math_utility.h"

namespace Quadrature
{

template<unsigned int D, typename Quadrature>
constexpr int TensorProductBase<D,Quadrature>::
numberEvaluations()
{
  return MathUtility::pow(Quadrature::numberEvaluations(),D);
}

// 1D sampling points
template<typename Quadrature>
std::array<std::array<double,1>,TensorProductBase<1,Quadrature>::numberEvaluations()> TensorProduct<1,Quadrature>::
samplingPoints()
{
  std::array<std::array<double,1>,TensorProductBase<1,Quadrature>::numberEvaluations()> samplingPoints;
  std::array<double, Quadrature::numberEvaluations()> samplingPoints1D = Quadrature::samplingPoints();
  for(int x=0; x<Quadrature::numberEvaluations(); x++)
  {
    samplingPoints[x][0] = samplingPoints1D[x];
  }
  return samplingPoints;
}

// 2D sampling points
template<typename Quadrature>
std::array<std::array<double,2>,TensorProductBase<2,Quadrature>::numberEvaluations()> TensorProduct<2,Quadrature>::
samplingPoints()
{
  std::array<std::array<double,2>,TensorProductBase<2,Quadrature>::numberEvaluations()> samplingPoints;
  std::array<double, Quadrature::numberEvaluations()> samplingPoints1D = Quadrature::samplingPoints();

  int samplingPointNo = 0;
  for(int y=0; y<Quadrature::numberEvaluations(); y++)
  {
    for(int x=0; x<Quadrature::numberEvaluations(); x++, samplingPointNo++)
    {
      samplingPoints[samplingPointNo] = {samplingPoints1D[x], samplingPoints1D[y]};
    }
  }
  return samplingPoints;
}

// 3D sampling points
template<typename Quadrature>
std::array<std::array<double,3>,TensorProductBase<3,Quadrature>::numberEvaluations()> TensorProduct<3,Quadrature>::
samplingPoints()
{
  std::array<std::array<double,3>,TensorProductBase<3,Quadrature>::numberEvaluations()> samplingPoints;
  std::array<double, Quadrature::numberEvaluations()> samplingPoints1D = Quadrature::samplingPoints();

  int samplingPointNo = 0;
  for(int z=0; z<Quadrature::numberEvaluations(); z++)
  {
    for(int y=0; y<Quadrature::numberEvaluations(); y++)
    {
      for(int x=0; x<Quadrature::numberEvaluations(); x++, samplingPointNo++)
      {
        samplingPoints[samplingPointNo] = {samplingPoints1D[x], samplingPoints1D[y], samplingPoints1D[z]};
      }
    }
  }
  return samplingPoints;
}

// 1D integration
template<typename Quadrature>
template<typename ValueType>
ValueType TensorProduct<1,Quadrature>::
computeIntegral(const std::array<ValueType, TensorProductBase<1,Quadrature>::numberEvaluations()> &evaluations)
{
  return Quadrature::computeIntegral(evaluations);
}

// 2D tensor product integration
template<typename Quadrature>
template<typename ValueType>
ValueType TensorProduct<2,Quadrature>::
computeIntegral(const std::array<ValueType, TensorProductBase<2,Quadrature>::numberEvaluations()> &evaluations)
{
  // integrate by calling Quadrature in each direction
  std::array<ValueType, Quadrature::numberEvaluations()> evaluationsY;
  for(int y = 0; y < Quadrature::numberEvaluations(); y++)
  {
    // index of first evaluation that belongs to the list for the current y
    size_t offset = y*Quadrature::numberEvaluations();

    evaluationsY[y] = Quadrature::template computeIntegral<ValueType>(evaluations.begin()+offset);
  }
  return Quadrature::computeIntegral(evaluationsY);
}

// 3D tensor product integration
template<typename Quadrature>
template<typename ValueType>
ValueType TensorProduct<3,Quadrature>::
computeIntegral(const std::array<ValueType, TensorProductBase<3,Quadrature>::numberEvaluations()> &evaluations)
{
  // integrate by calling Quadrature in each direction
  std::array<ValueType, Quadrature::numberEvaluations()> evaluationsZ;
  for(int z=0; z<Quadrature::numberEvaluations(); z++)
  {
    std::array<ValueType, Quadrature::numberEvaluations()> evaluationsY;
    for(int y=0; y<Quadrature::numberEvaluations(); y++)
    {
      // index of first evaluation that belongs to the list for the current y
      size_t offset = y*Quadrature::numberEvaluations() + z*Quadrature::numberEvaluations()*Quadrature::numberEvaluations();

      evaluationsY[y] = Quadrature::template computeIntegral<ValueType>(evaluations.begin()+offset);
    }
    evaluationsZ[z] = Quadrature::computeIntegral(evaluationsY);
  }
  return Quadrature::computeIntegral(evaluationsZ);
}

}
