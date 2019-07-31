import sys, os
from Package import Package

class boost(Package):

    def __init__(self, **kwargs):
        defaults = {
            'download_url': 'https://dl.bintray.com/boostorg/release/1.70.0/source/boost_1_70_0.tar.bz2',
        }
        defaults.update(kwargs)
        super(boost, self).__init__(**defaults)
        self.ext = '.cpp'
        self.sub_dirs = [
            ('include', 'lib'),
        ]
        self.headers = ['boost/optional.hpp']
        self.libs = [['boost_wserialization', 'boost_system', 'boost_program_options', 'boost_filesystem', 'boost_serialization']]
        self.check_text = r'''
#include <iostream>
#include <boost/optional.hpp>

#include <boost/filesystem.hpp>

using namespace boost::filesystem;


int main(int argc, char* argv[]) {
   
  path p{"/tmp"};

  std::cout << p.native() << '\n';
  std::cout << p.string() << '\n';
  std::wcout << p.wstring() << '\n';
  std::cout << p.generic_string() << '\n';
  std::wcout << p.generic_wstring() << '\n';
}
'''
        self.number_output_lines = 14292
        
    def check(self, ctx):
        env = ctx.env
        
        self.set_build_handler([
            './bootstrap.sh',
            #'./b2 link=static variant=release runtime-link=static cxxflags=-std=c++14 install --prefix=${PREFIX}  --with-system --with-filesystem --with-serialization --with-program_options'
            './b2 link=static variant=release runtime-link=static cxxstd=14 cxxflags=-std=c++14 install --prefix=${PREFIX}  --with-system --with-filesystem --with-serialization --with-program_options'
        ])
        
        ctx.Message('Checking for boost ...         ')
        self.check_options(env)

        res = super(boost, self).check(ctx)

        self.check_required(res[0], ctx)
        ctx.Result(res[0])
        return res[0]
