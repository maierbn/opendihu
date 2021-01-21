// g++ -std=c++14 std-simd.cpp -I $OPENDIHU_HOME/dependencies/vc/install/include/ -L$OPENDIHU_HOME/dependencies/vc/install/lib -lVc -I/store/software/std-simd && ./a.out
// g++ -std=c++17 std-simd.cpp -I $OPENDIHU_HOME/dependencies/vc/install/include/ -L$OPENDIHU_HOME/dependencies/vc/install/lib -lVc -I/store/software/std-simd && ./a.out

#include <iostream>
#include <cstdlib>

#if __cplusplus >= 201703L
#define HAVE_STDSIMD
#endif

#ifdef HAVE_STDSIMD

#include <experimental/simd>
// define compatibility layer such that std::simd appears like Vc
namespace Vc
{
  // define types
  template<typename T,int n>
  using array = std::experimental::fixed_size_simd<T,n>;
  using double_v = std::experimental::native_simd<double>;
  using int_v = std::experimental::native_simd<int>;
  
  // import math functions
  using std::experimental::abs;
  using std::experimental::exp;
  using std::experimental::log;
  
  // define constants
  double_v One = 1;
  double_v Zero = 0;
  
  // define inline if functions
  template<typename T>
  constexpr T iif(const typename T::mask_type& mask, const T& trueValue, const T& falseValue)
  {
    T result(falseValue);
    where(mask, result) = trueValue;
    return result;
  }
  
  template<typename T> 
  constexpr T iif (bool condition, const T &trueValue, const T &falseValue)
  {
    return condition ? trueValue : falseValue;
  }
  template<typename T> 
  constexpr T isnegative (bool condition, const T &trueValue, const T &falseValue)
  {
    return condition ? trueValue : falseValue;
  }
  
  template<typename T>
  typename T::mask_type isnegative(T x)
  {
    return x < 0;
  }
  
  template<typename Mask>
  int count(const Mask &mask)
  {
    int result(0);
    for (int i = 0; i < Mask::size(); i++)
      result += (int)(mask[i]);
    return result;
  }
  
  template<typename T>
  typename T::value_type min(const T &x)
  {
    typename T::value_type result(x[0]);
    for (int i = 1; i < T::size(); i++)
      result = std::min(result, x[i]);
    return result;
  }
  
  template<typename T>
  typename T::value_type max(const T &x)
  {
    typename T::value_type result(x[0]);
    for (int i = 1; i < T::size(); i++)
      result = std::max(result, x[i]);
    return result;
  }
}
#else

#include <Vc/Vc>

namespace Vc_1
{
  template<typename Mask>
  int count(const Mask &mask)
  {
    return mask.count();
  }
  
  template<typename T>
  typename T::value_type min(const T &x)
  {
    return x.min();
  }
  
  template<typename T>
  typename T::value_type max(const T &x)
  {
    return x.max();
  }
}

#endif

// output operator
std::ostream &operator<<(std::ostream &stream, const Vc::double_v &v)
{
  stream << "(";
  for (int i = 0; i < Vc::double_v::size(); i++)
  {
    if (i != 0)
      stream << ",";
    stream << v[i];
  }
  stream << ")";
  return stream;
}

int main(int argc, char *argv[])
{
  Vc::double_v a,b,c,d,f,g;
  a[0] = -1;
  a[1] = 2;
  b[0] = 10;
  b[1] = 1;
  c[0] = 0;
  c[1] = 0.1;
  d[0] = 10;
  d[1] = 10;
  f = Vc::double_v(Vc::Zero);
  g = Vc::double_v(Vc::One);
  
  std::cout << "C++: " << __cplusplus << std::endl;
#ifdef HAVE_STDSIMD
  std::cout << "have std::simd" << std::endl;
#endif
  std::cout << "double size: " << Vc::double_v::size() << ", int size: " << Vc::int_v::size() << std::endl;
  std::cout << "a=" << a << ", b=" << b << ", c=" << c << ", d=" << d << std::endl;
  std::cout << "a+b=" << a << ", log(a)=" << log(a) << ",abs(a)=" << abs(a) << ",exp(a)=" << exp(a) << std::endl;
  std::cout << "a+b=" << a << ", log(a)=" << Vc::log(a) << ",abs(a)=" << Vc::abs(a) << ",exp(a)=" << Vc::exp(a) << std::endl;
  std::cout << (double)a[0] << ",f=" << f << ",g=" << g << std::endl;
  where(a > 0, a) += 2;
  std::cout << "a=" << a << std::endl;
  Vc::double_v h = Vc::iif(b>5,Vc::double_v(Vc::Zero),a);
  Vc::double_v j = Vc::iif(true,Vc::double_v(Vc::Zero),a);
  std::cout << "h=" << h << ", j=" << j << std::endl;
  
  Vc::double_v e = log(c);
  if (any_of(isfinite(e)))
    std::cout << "isfinite" << std::endl;
  else
    std::cout << "not isfinite" << std::endl;
  
  Vc::double_v k;
  k[0] = -1;
  k[1] = 1;
  std::cout << "k=" << k << ", isnegative.count(k): " << Vc::count(Vc::isnegative(k)) << std::endl;
  //std::cout << "pow(b,2)=" << pow(b,2) << std::endl;
  std::cout << "b=" << b << ", min(b)=" << Vc::min(b) << ", max(b)=" << Vc::max(b) << std::endl;
  
  return EXIT_SUCCESS;
}
