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
#VERSION='https://github.com/zuhd-org/easyloggingpp/archive/v9.96.7.tar.gz'
#VERSION = 'https://github.com/muflihun/easyloggingpp/releases/download/v9.95.3/easyloggingpp_v9.95.3.tar.gz'
VERSION='https://github.com/maierbn/easyloggingpp/archive/master.zip'

class EasyLoggingPP(Package):
  
    def __init__(self, **kwargs):
        defaults = {
            'download_url': VERSION
        }
        defaults.update(kwargs)
        super(EasyLoggingPP, self).__init__(**defaults)
        self.ext = '.cpp'
        #self.sub_dirs = [
        #    ('include/mysql', 'lib'),
        #    ('include/mysql', 'lib64'),
        #]
        self.headers = ['easylogging++.h']
        self.sources = ['easylogging++.cpp']
        self.libs = []
        self.extra_libs = []
        self.check_text = check_text
        self.static = False
        
        #self.build_flags = '-std=c++14'
        
    def check(self, ctx):
        env = ctx.env
        ctx.Message('Checking for EasyLogging++ ... ')

        # Setup the build handler.
        self.set_build_handler([
            'mkdir -p  ${PREFIX}/include',
            'mkdir -p  ${PREFIX}/src',
            'cp ${SOURCE_DIR}/src/easylogging++.h ${PREFIX}/include',
            'cp ${SOURCE_DIR}/src/easylogging++.cc ${PREFIX}/src/easylogging++.cpp',
            env["CXX"]+' -std=c++11 -c ${PREFIX}/src/easylogging++.cpp -I${PREFIX}/include -DELPP_FEATURE_CRASH_LOG -DELPP_NO_DEFAULT_LOG_FILE -o ${PREFIX}/src/easylogging++.o',
            env["CXX"]+' -std=c++11 -c ${PREFIX}/src/easylogging++.cpp -I${PREFIX}/include -DELPP_FEATURE_CRASH_LOG -DELPP_NO_DEFAULT_LOG_FILE \
              -DELPP_DISABLE_DEBUG_LOGS -DELPP_DISABLE_VERBOSE_LOGS -DELPP_DISABLE_TRACE_LOGS -DELPP_DISABLE_LOGGING_FLAGS_FROM_ARG -DNDEBUG -o ${PREFIX}/src/easylogging++.release.o',
        ])
        
        self.check_options(env)

        res = super(EasyLoggingPP, self).check(ctx)

        self.check_required(res[0], ctx)
        ctx.Result(res[0])
        return res[0]
