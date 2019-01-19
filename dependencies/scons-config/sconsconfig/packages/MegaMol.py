import sys, os, multiprocessing, subprocess
from Package import Package

class MegaMol(Package):

  def __init__(self, **kwargs):
    defaults = {
      'download_url': 'https://github.com/UniStuttgart-VISUS/megamol/archive/master.zip'
    }
    defaults.update(kwargs)
    super(MegaMol, self).__init__(**defaults)
    self.ext = '.cpp'
    self.sub_dirs = [
        ('include', 'lib'),
    ]
    self.headers = ['Console.h']
    #self.libs = ['mysqlclient']
    #self.extra_libs = ['lapack', 'blas']
    self.check_text = r'''
      #include <stdlib.h>
      #include <stdio.h>
      #include <Console.h>
      int main(int argc, char* argv[])
      {
        megamol_main(argc,argv);
        return EXIT_SUCCESS;
      }
'''

    # Setup the build handler.
    #self.libs = ["MegaMol3.6m"]
    
  def check(self, ctx):
    env = ctx.env
    ctx.Message('Checking for MegaMol ...       ')
    self.check_options(env)

    # download and so on
    self.set_build_handler([
      'mkdir -p ${PREFIX}',
      "sed -i 's|int megamol_main(int argc, char\* argv\[\]) {|int main(int argc, char\* argv\[\]) {|g' ${SOURCE_DIR}/console/src/Console.cpp",  # reverse eventual change to Console.cpp
      'cd ${SOURCE_DIR} && mkdir -p build && cd build && \
      '+ctx.env["cmake"]+' -DCMAKE_INSTALL_PREFIX=${PREFIX} -DUSE_MPI=ON -DMPI_GUESS_LIBRARY_NAME=  .. \
      && make && make install',
      "cp ${SOURCE_DIR}/console/src/Console.cpp ${SOURCE_DIR}/console/src/Console.cpp.backup && \
      sed -i 's|int main(int argc, char\* argv\[\]) {|int megamol_main(int argc, char\* argv\[\]) {|g' ${SOURCE_DIR}/console/src/Console.cpp",  # replace main() by megamol_main in Console.cpp
      'mkdir -p ${PREFIX}/include && printf "#pragma once\nint megamol_main(int argc, char* argv[]);" > ${PREFIX}/include/Console.h',
      '!cd ${SOURCE_DIR}/build/console && echo "linking mmconsole will fail, but that is expected and okay!" && make'   # build Console.cpp.o with megamol_main() instead of main(), this fails when linking mmconsole, but that is okay
    ])
    self.number_output_lines = 12050

    res = super(MegaMol, self).check(ctx)
    self.check_required(res[0], ctx)
    
    # define custom linker flags, this needs the src tree present
    cwd = os.getcwd()
    opendihu_home = os.getcwd()
    ctx.Log("opendihu_home: {}\n".format(opendihu_home))
    
    path = os.path.join(opendihu_home, "dependencies/megamol/src")
    if os.path.exists(path):
      entries = os.listdir(path)
      if len(entries) == 1:
        if (os.path.isdir(os.path.join(path, entries[0]))):
          path = os.path.join(path, entries[0])
      
      megamol_link_path = os.path.join(path, "build/console/CMakeFiles/mmconsole.dir/link.txt")
      ctx.Log("megamol_link_path: {}\n".format(megamol_link_path))
      with open(megamol_link_path, "r") as f:
        linker_flags = f.readline()
        linker_flags = linker_flags[linker_flags.find(" "):]
        #print(linker_flags)

        new_linker_flags = []
        for item in linker_flags.split():
          if item == "../mmconsole" or item == "-o":
            continue
            
          if item[0] != "-" and item[0] != "/":
            item = os.path.join(os.path.join(opendihu_home, "dependencies/megamol/src/megamol-master/build/console/"), item)
          new_linker_flags.append(item)
        
        ctx.Log("new_linker_flags: {}\n".format(new_linker_flags))
        self.link_flags = new_linker_flags
        
      # run again with new linker flags and modified Console.cpp
      self.set_build_handler(None)

      res = super(MegaMol, self).check(ctx)
      self.check_required(res[0], ctx)
    else:
      ctx.Log("path {} does not exist\n".format(path))
      
    ctx.Result(res[0])
    return res[0]
