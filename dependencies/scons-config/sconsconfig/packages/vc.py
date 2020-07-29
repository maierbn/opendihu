import sys, os, multiprocessing
from .Package import Package
import subprocess

##
##  Handles Vc library
##
class Vc(Package):

  def __init__(self, **kwargs):
    defaults = {
        'download_url': 'https://github.com/VcDevel/Vc/archive/1.4.1.tar.gz',
    }
    defaults.update(kwargs)
    super(Vc, self).__init__(**defaults)
    self.ext = '.cpp'
    #self.sub_dirs = [
    #    ('include/mysql', 'lib'),
    #    ('include/mysql', 'lib64'),
    #]
    self.headers = ['Vc/Vc']
    self.libs = ['Vc']
    #self.extra_libs = ['lapack', 'blas']
    #self.set_rpath = False    # do not use dynamic linkage
    self.check_text = r'''
#include <Vc/Vc>
#include <iostream>
#include <iomanip>
using Vc::float_v;
int Vc_CDECL main()
{
    // allocate memory for our initial x and y coordinates. Note that you can also put it into a
    // normal float C-array but that you then must ensure alignment to Vc::VectorAlignment!
    Vc::Memory<float_v, 1000> x_mem;
    Vc::Memory<float_v, 1000> y_mem;
    Vc::Memory<float_v, 1000> r_mem;
    Vc::Memory<float_v, 1000> phi_mem;
    
        // fill the memory with values from -1.f to 1.f
    for (size_t i = 0; i < x_mem.vectorsCount(); ++i) {
        x_mem.vector(i) = float_v::Random() * 2.f - 1.f;
        y_mem.vector(i) = float_v::Random() * 2.f - 1.f;
    }
    
    // calculate the polar coordinates for all coordinates and overwrite the euclidian coordinates
    // with the result
    for (size_t i = 0; i < x_mem.vectorsCount(); ++i) {
        const float_v x = x_mem.vector(i);
        const float_v y = y_mem.vector(i);
        r_mem.vector(i) = Vc::sqrt(x * x + y * y);
        float_v phi = Vc::atan2(y, x) * 57.295780181884765625f; // 180/pi
        phi(phi < 0.f) += 360.f;
        phi_mem.vector(i) = phi;
    }
    
    // print the results
    for (size_t i = 0; i < x_mem.entriesCount(); ++i) {
        std::cout << std::setw(3) << i << ": ";
        std::cout << std::setw(10) << x_mem[i] << ", " << std::setw(10) << y_mem[i] << " -> ";
        std::cout << std::setw(10) << r_mem[i] << ", " << std::setw(10) << phi_mem[i] << '\n';
    }
    return 0;
}
'''

  def check(self, ctx):
      
    # get number of available processors
    p = multiprocessing.cpu_count()

    # Setup the build handler.
    self.set_build_handler([
      "mkdir -p build; cd build; cmake -DCMAKE_INSTALL_PREFIX=${PREFIX} -DBUILD_TESTING=OFF -DBUILD_SHARED_LIBS=ON -DCMAKE_BUILD_TYPE=Release ..",
      "cd build; make -j{}".format(p),
      "cd build; make install",
    ])
    self.number_output_lines = 166
    
    env = ctx.env
    ctx.Message('Checking for Vc ...            ')
    self.check_options(env)

    res = super(Vc, self).check(ctx)

    self.check_required(res[0], ctx)
    ctx.Result(res[0])
    return res[0]
