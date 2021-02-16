import sys, os, multiprocessing
from .Package import Package
import subprocess

##
##  Handles std_simd library
##
class std_simd(Package):

  def __init__(self, **kwargs):
    defaults = {
        'download_url': 'https://github.com/VcDevel/std-simd/archive/master.zip',
    }
    defaults.update(kwargs)
    super(std_simd, self).__init__(**defaults)
    self.ext = '.cpp'
    #self.sub_dirs = [
    #    ('include/mysql', 'lib'),
    #    ('include/mysql', 'lib64'),
    #]
    self.headers = ['experimental/simd']
    self.libs = []
    #self.extra_libs = ['lapack', 'blas']
    #self.set_rpath = False    # do not use dynamic linkage
  
    self.check_text = r'''
// g++ -std=c++14 std-simd.cpp -I $OPENDIHU_HOME/dependencies/vc/install/include/ -L$OPENDIHU_HOME/dependencies/vc/install/lib -lVc -I/store/software/std-simd && ./a.out
// g++ -std=c++17 std-simd.cpp -I $OPENDIHU_HOME/dependencies/vc/install/include/ -L$OPENDIHU_HOME/dependencies/vc/install/lib -lVc -I/store/software/std-simd && ./a.out
#include <iostream>
#include <cstdlib>
#include <vc_or_std_simd.h>

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
  
  return EXIT_SUCCESS;
}

 
'''

  def check(self, ctx):
      
    # Setup the build handler.
    self.set_build_handler([
      "mkdir -p ${PREFIX}",
      "cd ${PREFIX} && ln -s ${SOURCE_DIR} include",
      "cd ${PREFIX}/include && wget https://github.com/maierbn/std_simd_vc_wrapper/archive/main.zip && unzip main.zip && mv std_simd_vc_wrapper-main/* ."
    ])
    self.number_output_lines = 10
    
    env = ctx.env
    ctx.Message('Checking for std_simd ...      ')
    self.check_options(env)

    res = super(std_simd, self).check(ctx)

    self.check_required(res[0], ctx)
    ctx.Result(res[0])
    return res[0]
