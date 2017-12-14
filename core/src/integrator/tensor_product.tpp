
namespace Integrator 
{

template <unsigned int D, typename Integrator> 
constexpr int TensorProductBase<D,Integrator>::
numberEvaluations()
{
  return pow(Integrator::numberEvaluations(),D);
}

template <unsigned int D, typename Integrator> 
constexpr int TensorProductBase<D,Integrator>::
samplingArraySize()
{
  return numberEvaluations()*D;
}

// 1D sampling points
template <typename Integrator>
std::array<double, TensorProductBase<1,Integrator>::samplingArraySize()> TensorProduct<1,Integrator>::
samplingPoints()
{
  return Integrator::samplingPoints();
}

// 2D sampling points
template <typename Integrator>
std::array<double, TensorProductBase<2,Integrator>::samplingArraySize()> TensorProduct<2,Integrator>::
samplingPoints()
{
  std::array<double, TensorProductBase<2,Integrator>::samplingArraySize()> samplingPoints;
  std::array<double, Integrator::numberEvaluations()> samplingPoints1D = Integrator::samplingPoints();
  
  int samplingPointNo = 0;
  for(int y=0; y<Integrator::numberEvaluations(); y++)
  {
    for(int x=0; x<Integrator::numberEvaluations(); x++, samplingPointNo++)
    {
      samplingPoints[samplingPointNo*2 + 0] = samplingPoints1D[x];
      samplingPoints[samplingPointNo*2 + 1] = samplingPoints1D[y];
    }
  }
  return samplingPoints;
}

// 3D sampling points
template <typename Integrator>
std::array<double, TensorProductBase<3,Integrator>::samplingArraySize()> TensorProduct<3,Integrator>::
samplingPoints()
{
  std::array<double, TensorProductBase<3,Integrator>::samplingArraySize()> samplingPoints;
  std::array<double, Integrator::numberEvaluations()> samplingPoints1D = Integrator::samplingPoints();
    
  int samplingPointNo = 0;
  for(int z=0; z<Integrator::numberEvaluations(); z++)
  {
    for(int y=0; y<Integrator::numberEvaluations(); y++)
    {
      for(int x=0; x<Integrator::numberEvaluations(); x++, samplingPointNo++)
      {
        samplingPoints[samplingPointNo*3 + 0] = samplingPoints1D[x];
        samplingPoints[samplingPointNo*3 + 1] = samplingPoints1D[y];
        samplingPoints[samplingPointNo*3 + 2] = samplingPoints1D[z];
      }
    }
  }
  return samplingPoints;
}

// 1D integration
template <typename Integrator>
double TensorProduct<1,Integrator>::
integrate(std::array<double, TensorProductBase<1,Integrator>::numberEvaluations()> &evaluations)
{
  return Integrator::integrate(evaluations);
}

// 2D tensor product integration
template <typename Integrator>
double TensorProduct<2,Integrator>::
integrate(std::array<double, TensorProductBase<2,Integrator>::numberEvaluations()> &evaluations)
{
  // integrate by calling Integrator in each direction
  std::array<double, Integrator::numberEvaluations()> evaluationsY;
  for(int y=0; y<Integrator::numberEvaluations(); y++)
  {
    size_t offset = y*Integrator::numberEvaluations();
    evaluationsY[y] = Integrator::integrate(evaluations.data()+offset);
  }
  return Integrator::integrate(evaluationsY);
}

// 3D tensor product integration
template <typename Integrator>
double TensorProduct<3,Integrator>::
integrate(std::array<double, TensorProductBase<3,Integrator>::numberEvaluations()> &evaluations)
{
  // integrate by calling Integrator in each direction
  std::array<double, Integrator::numberEvaluations()> evaluationsZ;
  for(int z=0; z<Integrator::numberEvaluations(); z++)
  {
    std::array<double, Integrator::numberEvaluations()> evaluationsY;
    for(int y=0; y<Integrator::numberEvaluations(); y++)
    {
      size_t offset = y*Integrator::numberEvaluations() + z*Integrator::numberEvaluations()*Integrator::numberEvaluations();
      evaluationsY[y] = Integrator::integrate(evaluations.data()+offset);
    }
    evaluationsZ[z] = Integrator::integrate(evaluationsY);
  }
  return Integrator::integrate(evaluationsZ);
}

}