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
VERSION='https://github.com/zuhd-org/easyloggingpp/archive/v9.96.7.tar.gz'
#VERSION = 'https://github.com/muflihun/easyloggingpp/releases/download/v9.95.3/easyloggingpp_v9.95.3.tar.gz'

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
        self.sources = ['easylogging++.cc']
        self.libs = []
        self.extra_libs = []
        self.check_text = check_text
        self.static = False
        
        #self.build_flags = '-std=c++14'
        
    def check(self, ctx):
        env = ctx.env
        ctx.Message('Checking for EasyLoggingPP ... ')
        
        flag = ""
        if env["CXX"] == "g++":
          flag = " -std=c++11"
        elif "pgc++" in env["CXX"]:
          flag = " -std=c++11"
        if "v9.96.7" in VERSION:
          sub_dir = '/src'
          linenumber1 = '1984'
          linenumber2 = '1985'
          ln3 = '396'
          ln4 = '397'
          ln5 = '398'
          ln6 = '399'
          ln7 = '400'
          ln8 = '401'
        else:
          sub_dir = ''
          linenumber1 = '1906'
          linenumber2 = '1907'
          ln3 = '389'
          ln4 = '390'
          ln5 = '391'
          ln6 = '392'
          ln7 = '393'
          ln8 = '394'



        # Setup the build handler.
        self.set_build_handler([
            'rm -f README.txt',
            'mkdir -p  ${PREFIX}/include',
            'mkdir -p  ${PREFIX}/src',
            'cp ${SOURCE_DIR}' + sub_dir + '/easylogging++.h ${PREFIX}/include',
            'cp ${SOURCE_DIR}' + sub_dir + '/easylogging++.cc ${PREFIX}/src/easylogging++.cpp',
            'sed -i \'' + linenumber1 + 'i      m_modules.insert(std::make_pair(ss.str(), level));\' ${PREFIX}/src/easylogging++.cpp',
            'sed -i \'' + linenumber2 + 'i      addSuffix(ss, ".tpp", ".hh");\' ${PREFIX}/src/easylogging++.cpp',
            'sed -i \'' + ln3 + 'i #if __cplusplus >= 201103L\' ${PREFIX}/include/easylogging++.h',
            'sed -i \'' + ln4 + 'i __attribute__((weak))\' ${PREFIX}/include/easylogging++.h',
            'sed -i \'' + ln5 + 'i void operator delete(void * ptr, std::size_t){ ::operator delete(ptr);}\' ${PREFIX}/include/easylogging++.h',
            'sed -i \'' + ln6 + 'i __attribute__((weak))\' ${PREFIX}/include/easylogging++.h',
            'sed -i \'' + ln7 + 'i void operator delete[](void * ptr, std::size_t){ ::operator delete(ptr);}\' ${PREFIX}/include/easylogging++.h',
            'sed -i \'' + ln8 + 'i #endif  // __cplusplus >= 201103L\' ${PREFIX}/include/easylogging++.h',
            env["CXX"]+flag+' -c ${PREFIX}/src/easylogging++.cpp -I${PREFIX}/include -DELPP_FEATURE_CRASH_LOG -DELPP_NO_DEFAULT_LOG_FILE -o ${PREFIX}/src/easylogging++.o',
        ])
        
        self.check_options(env)

        res = super(EasyLoggingPP, self).check(ctx)

        self.check_required(res[0], ctx)
        ctx.Result(res[0])
        return res[0]
