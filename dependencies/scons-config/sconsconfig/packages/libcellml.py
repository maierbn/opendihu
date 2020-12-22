import sys, os, multiprocessing
from .Package import Package

class libcellml(Package):

    def __init__(self, **kwargs):
        defaults = {
            'download_url': 'https://github.com/cellml/libcellml/archive/develop.zip',
        }
        defaults.update(kwargs)
        super(libcellml, self).__init__(**defaults)
        self.ext = '.cpp'
        #self.sub_dirs = [
        #    ('include/mysql', 'lib'),
        #    ('include/mysql', 'lib64'),
        #]
        #self.headers = ['mysql.h']
        #self.libs = ['mysqlclient']
        #self.extra_libs = ['lapack', 'blas']
        self.check_text = r'''
#include <iostream>
#include <cstdlib>
#include <libcellml>
#include <memory>
int main(int argc, char* argv[]) {
  const std::string e =
          "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
          "<model xmlns=\"http://www.cellml.org/cellml/2.0#\">"
              "<units name=\"valid_name\"/>"
          "</model>";

  libcellml::Model m;

  libcellml::UnitsPtr u = std::make_shared<libcellml::Units>();;
  u->setName("valid_name");

  m.addUnits(u);

  libcellml::Printer printer;
  const std::string a = printer.printModel(m);
  return EXIT_SUCCESS;
}
'''

        # get number of available processors
        p = multiprocessing.cpu_count()

        # Setup the build handler.
        self.set_build_handler([
          "mkdir -p ${PREFIX}/../build ",
          "cd ${PREFIX}/../build && cmake \
            -DBUILD_TYPE=Release \
            -DINSTALL_PREFIX=${PREFIX} \
            -DLIBXML2_LIBRARIES=../../libxml2/install/lib/libxml2.a \
            -DLIBXML2_INCLUDE_DIR=../../libxml2/install/include \
            -DLIBCELLML_BUILD_SHARED=Off \
            -DLIBCELLML_COVERAGE=Off \
            -DLIBCELLML_MEMCHECK=Off \
            -DLIBCELLML_UNIT_TESTS=Off \
            ${SOURCE_DIR} && make && make install",
          #"ln -s ${PREFIX}/include/libcellml/module/libcellml ${PREFIX}/include/libcellml"
        ])
        self.number_output_lines = 84
          
        self.libs = ["cellml"]
        self.headers = ["libcellml/model.h", "libcellml"]

    def check(self, ctx):
        env = ctx.env
        ctx.Message('Checking for libcellml ... ')
        self.check_options(env)

        res = super(libcellml, self).check(ctx)

        self.check_required(res[0], ctx)
        ctx.Result(res[0])
        return res[0]
