import sys, os, multiprocessing, subprocess
from .Package import Package

class OpenCOR(Package):

  def __init__(self, **kwargs):
    defaults = {
      'download_url': 'https://github.com/opencor/opencor/releases/download/v0.6/OpenCOR-0-6-Linux.tar.gz'
      
    }
    defaults.update(kwargs)
    super(OpenCOR, self).__init__(**defaults)
    self.ext = '.cpp'
    self.sub_dirs = [
      ('include', 'lib'),
    ]
    
    self.check_text = r'''
      #include <stdlib.h>
      #include <stdio.h>
      #include <iostream>
      #include <sstream>
      #include "opencor.h"
      
      int main(int argc, char* argv[])
      {
        std::stringstream command;
        command << OPENCOR_BINARY << " --version";
        int ret = system(command.str().c_str());
        if (ret == 0)
        {
          std::cout << "opencor found.";
          return EXIT_SUCCESS;
        }
        else
        {
          std::cout << "opencor not found";
          return EXIT_FAILURE;
        }
      }
    '''
  
    # Setup the build handler.
    self.libs = []
    self.headers = ['opencor.h']
    
    self.set_build_handler([
      'mkdir -p ${PREFIX}/bin ${PREFIX}/include',
      'ln -s ${SOURCE_DIR}/formats ${PREFIX}/include/formats',
      'echo "#define OPENCOR_BINARY \\\"${PREFIX}/bin/opencor\\\"" > ${PREFIX}/include/opencor.h',
      'echo "#define OPENCOR_FORMATS_DIRECTORY \\\"${PREFIX}/include/formats\\\"" >> ${PREFIX}/include/opencor.h',
      '$echo "#!/bin/bash" > ${PREFIX}/bin/opencor && echo \'${SOURCE_DIR}/bin/OpenCOR $*\' >> ${PREFIX}/bin/opencor && chmod +x ${PREFIX}/bin/opencor'
    ])

  def check(self, ctx):
    env = ctx.env
    ctx.Message('Checking for OpenCOR ...       ')
    self.check_options(env)

    res = super(OpenCOR, self).check(ctx)

    self.check_required(res[0], ctx)
    ctx.Result(res[0])
    return res[0]
