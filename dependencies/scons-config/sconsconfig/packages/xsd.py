import sys, os, multiprocessing
from .Package import Package


class xsd(Package):

  def __init__(self, **kwargs):
    defaults = {
      'download_url': 'http://www.codesynthesis.com/download/xsd/4.0/linux-gnu/x86_64/xsd-4.0.0-x86_64-linux-gnu.tar.bz2',
    }
    defaults.update(kwargs)
    super(xsd, self).__init__(**defaults)
    self.libs=[]
    self.extra_libs=[]
    self.ext = '.cpp'
    self.check_text = r'''
int main(int argc, char* argv[])
{
  return 0;
}
'''
    # Setup the build handler. I'm going to assume this will work for all architectures.
    self.set_build_handler([
        'rm -rf ${PREFIX}; ln -s ${SOURCE_DIR} ${PREFIX}'
    ])

  def check(self, ctx):
    env = ctx.env
    ctx.Message('Checking for xsd ...           ')
    self.check_options(env)

    res = super(xsd, self).check(ctx)

    self.check_required(res[0], ctx)
    ctx.Result(res[0])
    return res[0]
