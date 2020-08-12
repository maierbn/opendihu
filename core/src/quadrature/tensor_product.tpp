#include "quadrature/tensor_product.h"

#include "utility/math_utility.h"

namespace Quadrature
{

template<unsigned int D, typename Quadrature>
constexpr int TensorProductBase<D,Quadrature>::
numberEvaluations()
{
  // compile-time power function
  return MathUtility::powConst(Quadrature::numberEvaluations(),D);
}

// 1D sampling points
template<typename Quadrature>
std::array<std::array<double,1>,TensorProductBase<1,Quadrature>::numberEvaluations()> TensorProduct<1,Quadrature>::
samplingPoints()
{
  std::array<std::array<double,1>,TensorProductBase<1,Quadrature>::numberEvaluations()> samplingPoints;
  std::array<double, Quadrature::numberEvaluations()> samplingPoints1D = Quadrature::samplingPoints();
  for (int x=0; x<Quadrature::numberEvaluations(); x++)
  {
    samplingPoints[x][0] = samplingPoints1D[x];
  }
  return samplingPoints;
}

// 2D sampling points
template<typename Quadrature>
std::array<Vec2,TensorProductBase<2,Quadrature>::numberEvaluations()> TensorProduct<2,Quadrature>::
samplingPoints()
{
  std::array<Vec2,TensorProductBase<2,Quadrature>::numberEvaluations()> samplingPoints;
  std::array<double, Quadrature::numberEvaluations()> samplingPoints1D = Quadrature::samplingPoints();

  int samplingPointNo = 0;
  for (int j = 0; j < Quadrature::numberEvaluations(); j++)
  {
    for (int i = 0; i < Quadrature::numberEvaluations(); i++, samplingPointNo++)
    {
      samplingPoints[samplingPointNo] = {samplingPoints1D[i], samplingPoints1D[j]};   // x y
    }
  }
  return samplingPoints;
}

// 3D sampling points
template<typename Quadrature>
std::array<Vec3,TensorProductBase<3,Quadrature>::numberEvaluations()> TensorProduct<3,Quadrature>::
samplingPoints()
{
  std::array<Vec3,TensorProductBase<3,Quadrature>::numberEvaluations()> samplingPoints;
  std::array<double, Quadrature::numberEvaluations()> samplingPoints1D = Quadrature::samplingPoints();

  int samplingPointNo = 0;
  for (int k = 0; k < Quadrature::numberEvaluations(); k++)
  {
    for (int j = 0; j < Quadrature::numberEvaluations(); j++)
    {
      for (int i = 0; i < Quadrature::numberEvaluations(); i++, samplingPointNo++)
      {
        samplingPoints[samplingPointNo] = {samplingPoints1D[i], samplingPoints1D[j], samplingPoints1D[k]};  // x y z
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
  const std::array<double,Quadrature::numberEvaluations()> weights = Quadrature::quadratureWeights();

  ValueType result{};
  for (int j = 0; j < Quadrature::numberEvaluations(); j++)
  {
    for (int i = 0; i < Quadrature::numberEvaluations(); i++)
    {
      const int index = j*Quadrature::numberEvaluations() + i;
      result += weights[j]*weights[i] * evaluations[index];
    }
  }
  return result;
  /*
  // integrate by calling Quadrature in each direction
  std::array<ValueType, Quadrature::numberEvaluations()> evaluationsY;
  for (int y = 0; y < Quadrature::numberEvaluations(); y++)
  {
    // index of first evaluation that belongs to the list for the current y
    size_t offset = y*Quadrature::numberEvaluations();

    evaluationsY[y] = Quadrature::template computeIntegral<ValueType>(evaluations.begin()+offset);
  }
  return Quadrature::computeIntegral(evaluationsY);
  */
}

// 3D tensor product integration
template<typename Quadrature>
template<typename ValueType>
ValueType TensorProduct<3,Quadrature>::
computeIntegral(const std::array<ValueType, TensorProductBase<3,Quadrature>::numberEvaluations()> &evaluations)
{
  const std::array<double,Quadrature::numberEvaluations()> weights = Quadrature::quadratureWeights();

  ValueType result{};
  for (int k = 0; k < Quadrature::numberEvaluations(); k++)
  {
    for (int j = 0; j < Quadrature::numberEvaluations(); j++)
    {
      for (int i = 0; i < Quadrature::numberEvaluations(); i++)
      {
        const int index = k*Quadrature::numberEvaluations()*Quadrature::numberEvaluations() + j*Quadrature::numberEvaluations() + i;
        result += weights[k]*weights[j]*weights[i] * evaluations[index];
      }
    }
  }
  return result;
/*
  // integrate by calling Quadrature in each direction
  std::array<ValueType, Quadrature::numberEvaluations()> evaluationsZ;
  for (int z=0; z<Quadrature::numberEvaluations(); z++)
  {
    std::array<ValueType, Quadrature::numberEvaluations()> evaluationsY;
    for (int y=0; y<Quadrature::numberEvaluations(); y++)
    {
      // index of first evaluation that belongs to the list for the current y
      size_t offset = y*Quadrature::numberEvaluations() + z*Quadrature::numberEvaluations()*Quadrature::numberEvaluations();

      evaluationsY[y] = Quadrature::template computeIntegral<ValueType>(evaluations.begin()+offset);
    }
    evaluationsZ[z] = Quadrature::computeIntegral(evaluationsY);
  }
  return Quadrature::computeIntegral(evaluationsZ);
  */
}

}  // namespace
