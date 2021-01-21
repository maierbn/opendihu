
#include <Vc/Vc>
#include <iostream>
#include <cstdlib>
#include <vector>
#include <algorithm>
#include <array>
#include <cmath>
#include <chrono>

#include "easylogging++.h"
#include "semt/Semt.h"

#include "opendihu.h"

using namespace std;
/*
struct FiberPointData
{
  Vc::double_v values[8];
};
Vc_DECLARE_ALLOCATOR(FiberPointData)

int f(int a)
{
  double v = 0;
  for(long long int j = 0; j < 1e8; j++)
  {
    if (j % 17 == 5 || j % 13 == 2 || j%7 == 6)
      v += 2*a;
    if (j%2 == 0)
      v *= 2.;
    if (j%3 == 2)
      v /= 3.21;
  }
  return v;
}

class A
{
public:

  void init()
  {
    v_.resize(5);
    std::iota(v_.begin(), v_.end(), 0);
  }
  
  void print()
  {
    LOG(DEBUG) << "v: ";
    for (int i = 0; i < v_.size(); i++)
      LOG(DEBUG) << v_[i];
  }

  void getV(std::vector<int> &&v)
  {
    v = v_;
  }
protected:

  std::vector<int> v_;
};
*/
double gen()
{
  double d = rand()*1e-5;
  if (rand()%2 == 0)
  {
    return -d;
  }
  return d;
}
int main(int argc, char *argv[])
{
  // test Vc
  
  Vc::array<double,14> data;
  std::iota(data.begin(), data.end(), 10);
  
  LOG(INFO) << "data: " << data;
  //auto indexes = Vc::double_v::IndexType::IndexesFromZero();
  auto indexes = Vc::double_v([](int n) { return n; });
  Vc::double_v gathered = data[indexes];  // gathered == [0, 1, 2, ...]
  LOG(INFO) << "gathered: " << gathered;
  
  //Vc::double_v gathered2 = data[Vc::IndexesFromZero];  // gathered == [0, 1, 2, ...]
  //LOG(DEBUG) << "gathered2: " << gathered2;
  
  auto indexes2 = Vc::double_v([](int n) { return n+1; });
  LOG(INFO) << StringUtility::demangle(typeid(indexes2).name());
  Vc::double_v gathered3 = data[indexes2];  // gathered == [0, 1, 2, ...]
  LOG(INFO) << "indexes2: " << indexes2 << ", gathered3: " << gathered3 << ",";

  for (int i = 0; i < int(ceil((float)data.size() / Vc::double_v::size())); i++)
  {
    auto indexes = Vc::double_v([&i](int n) { return i*Vc::double_v::size() + n; });
    Vc::double_v v = data[indexes];
    LOG(INFO) << "indexes: " << indexes << ", v: " << v;
  }*/

  VecD<3,Vc::double_v> a{Vc::double_v(Vc::Zero),Vc::double_v(Vc::Zero),Vc::double_v(Vc::One)};
  a[0][0] = 0.1;
  a[0][1] = 0.2;
  a[0][2] = 0.3;
  a[0][3] = 0.4;
  Vc::double_v inverseNorm = 1./MathUtility::norm<3,Vc::double_v>(a);
  
  LOG(INFO) << a;
  LOG(INFO) << "inverseNorm: " << inverseNorm;
  LOG(INFO) << a*inverseNorm;
  VecD<3,Vc::double_v> normalizedA = a*inverseNorm;
  LOG(INFO) << MathUtility::norm<3,Vc::double_v>(normalizedA);
  
  MathUtility::normalize<3,Vc::double_v>(a);
  LOG(INFO) << a;
  exit(0);
#if 0  
  Vc::double_v v;
  for (int i = 0; i < Vc::double_v::size(); i++)
  {
    v[i] = i * (i%2==0? 1 : -1);
  }
  
  Vc::double_v epsilon = 0;
  LOG(INFO) << "epsilon: " << epsilon << ", v: " << v << ", abs: " << Vc::abs(v);
  Vc::where(Vc::abs(v) < 2.5) | epsilon = 1;
  LOG(INFO) << "epsilon: " << epsilon;

#if 0
  // measurement
  const long long N = 1e5;
  Vc::array<double,N> valuesA;
  Vc::array<double,N> valuesB;
  std::array<double,N> valuesArrayA;
  std::array<double,N> valuesArrayB;
  
  for (long long i = 0; i < N; i++)
  {
    valuesA[i] = (double)(gen());
    valuesB[i] = (double)(gen());
    valuesArrayA[i] = valuesA[i];
    valuesArrayB[i] = valuesB[i];
  }
  
  Vc::double_v cc = 0;
  
	auto start0 = chrono::steady_clock::now();

  const long long nEntries = int(ceil((float)valuesA.size() / Vc::double_v::size()));
  for (long long i = 0; i < nEntries; i++)
  {
    auto indexes = Vc::double_v([&i](int n) { return i*Vc::double_v::size() + n; });
    Vc::double_v a = valuesA[indexes];
    Vc::double_v b = valuesB[indexes];
    cc += a*b;
    //LOG(INFO) << "indexes: " << indexes << ", v: " << v;
  }
  
  double result = 0;
  for (long long i = 0; i < Vc::double_v::size(); i++)
  {
    result += cc[i];
  }
  
	auto end0 = chrono::steady_clock::now();
  
  LOG(INFO) << "resulta: " << result;
	auto start1 = chrono::steady_clock::now();
  
  result = 0;
  for (long long i = 0; i < valuesA.size(); i++)
  {
    result += valuesArrayA[i] * valuesArrayB[i];
  }
  LOG(INFO) << "resultb: " << result;
  
	auto end1 = chrono::steady_clock::now();

    
  LOG(INFO) << "duration Vc : " << chrono::duration_cast<chrono::microseconds>(end0 - start0).count();
  LOG(INFO) << "duration std: " << chrono::duration_cast<chrono::microseconds>(end1 - start1).count();

  
  // loop over elements
  for (int elementNoLocal = 0; elementNoLocal < nElementsLocal; elementNoLocal += Vc::double_v::size())
  {
    Vc::int_v elementNoLocalv([elementNoLocal, nElementsLocal](int i)
    {
      return (i >= Vc::double_v::size() || elementNoLocal+i >= nElementsLocal? -1: elementNoLocal+i); 
    });
  
    // get indices of element-local dofs
    
    std::array<Vc::int_v,nDisplacementsDofsPerElement> dofNosLocalv = displacementsFunctionSpace->getElementDofNosLocal(elementNoLocalv);
    
    LOG(INFO) << "element " << elementNoLocalv << ", type: " << StringUtility::demangle(typeid(decltype(elementNoLocalv)).name());
    LOG(INFO) << "dofNosLocal: " << dofNosLocalv;
  }

  
  int nFiberPoints = 100;
  int nVcVectors = (nFiberPoints + Vc::double_v::Size - 1) / Vc::double_v::Size;

  std::vector<FiberPointData> values(nVcVectors);
  
  for (int i = 0; i < nVcVectors; i++)
  {
    values[i].values[0] = 8.0;
  }
  

  return EXIT_SUCCESS;
}
