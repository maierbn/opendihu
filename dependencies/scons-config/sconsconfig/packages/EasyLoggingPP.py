import sys, os
from distutils import sysconfig
from Package import Package
from Package import Package

check_text = r'''
#include <stdlib.h>
#include <stdio.h>
#include "easylogging++.h"

INITIALIZE_EASYLOGGINGPP

int main(int argc, char* argv[])
{

   START_EASYLOGGINGPP(argc, argv);

   LOG(INFO) << "hi";
   return EXIT_SUCCESS;
}
'''

class EasyLoggingPP(Package):
  
    def __init__(self, **kwargs):
        defaults = {
            'download_url': 'https://github.com/muflihun/easyloggingpp/releases/download/v9.95.3/easyloggingpp_v9.95.3.tar.gz'
        }
        defaults.update(kwargs)
        super(EasyLoggingPP, self).__init__(**defaults)
        self.ext = '.cpp'
        #self.sub_dirs = [
        #    ('include/mysql', 'lib'),
        #    ('include/mysql', 'lib64'),
        #]
        self.headers = ['easylogging++.h']
        self.sources = ['easylogging++.cc']
        self.libs = []
        self.extra_libs = []
        self.check_text = check_text
        self.static = False
        
        # Setup the build handler.
        self.set_build_handler([
            'rm README.txt',
            'mkdir -p  ${PREFIX}/include',
            'mkdir -p  ${PREFIX}/src',
            'mv ${SOURCE_DIR}/easylogging++.h ${PREFIX}/include',
            'mv ${SOURCE_DIR}/easylogging++.cc ${PREFIX}/src',
            'g++ -c ${PREFIX}/src/easylogging++.cc -I${PREFIX}/include -std=c++11 -DELPP_FEATURE_CRASH_LOG -o ${PREFIX}/src/easylogging++.o',
        ])
        self.build_flags = '-std=c++14'
        
    def check(self, ctx):
        env = ctx.env
        ctx.Message('Checking for EasyLoggingPP ... ')
        self.check_options(env)

        res = super(EasyLoggingPP, self).check(ctx)

        self.check_required(res[0], ctx)
        ctx.Result(res[0])
        return res[0]
