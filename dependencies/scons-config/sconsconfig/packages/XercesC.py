import sys, os, multiprocessing
from .Package import Package


class XercesC(Package):

  def __init__(self, **kwargs):
    defaults = {
      'download_url': 'https://archive.apache.org/dist/xerces/c/3/sources/xerces-c-3.1.4.tar.gz',
    }
    defaults.update(kwargs)
    super(XercesC, self).__init__(**defaults)
    self.libs=['xerces-c']
    self.extra_libs=[]
    self.ext = '.cpp'
    self.check_text = r'''
#include <xercesc/util/PlatformUtils.hpp>
// Other include files, declarations, and non-Xerces-C++ initializations.

using namespace xercesc;

int main(int argc, char* argv[])
{
  try {
    XMLPlatformUtils::Initialize();
  }
  catch (const XMLException& toCatch) {
    // Do your failure processing here
    return 1;
  }

  // Do your actual work with Xerces-C++ here.

  XMLPlatformUtils::Terminate();

  // Other terminations and cleanup.
  return 0;
}
'''
    # get number of available processors
    p = multiprocessing.cpu_count()

    # Setup the build handler. I'm going to assume this will work for all architectures.
    self.set_build_handler([
        './configure --prefix=${PREFIX}',
        'make -j{}'.format(p),
        'make install'
    ])

  def check(self, ctx):
    env = ctx.env
    ctx.Message('Checking for XercesC ...       ')
    self.check_options(env)

    res = super(XercesC, self).check(ctx)

    self.check_required(res[0], ctx)
    ctx.Result(res[0])
    return res[0]
