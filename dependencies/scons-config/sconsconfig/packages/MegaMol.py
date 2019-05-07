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
    self.headers = ['Console.h', 'zmq.h']
    self.libs = ['zmq']
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

    # get number of available processors
    p = multiprocessing.cpu_count()

    # download and so on
    self.set_build_handler([
      'mkdir -p ${PREFIX} && mkdir -p ${DEPENDENCIES_DIR}/megamol/ospray',  # create folders
      'cd ${DEPENDENCIES_DIR}/megamol && wget https://datapacket.dl.sourceforge.net/project/ispcmirror/v1.11.0/ispc-v1.11.0-linux.tar.gz && tar xf ispc-v1.11.0-linux.tar.gz',
      'export PATH=$PATH:${DEPENDENCIES_DIR}/megamol/ispc-v1.11.0-linux/bin && cd ${DEPENDENCIES_DIR}/megamol && git clone https://github.com/embree/embree.git embree; cd embree && git checkout v3.5.2 && mkdir -p build && cd build && cmake -DCMAKE_C_COMPILER='+ctx.env["CC"]+' -DCMAKE_CXX_COMPILER='+ctx.env["CXX"]+' -DEMBREE_TUTORIALS=OFF -DCMAKE_INSTALL_PREFIX=${PREFIX} .. && make -j '+str(p)+' && make install', # build embree
      'export PATH=$PATH:${DEPENDENCIES_DIR}/megamol/ispc-v1.11.0-linux/bin && cd ${DEPENDENCIES_DIR}/megamol && git clone https://github.com/ospray/ospray.git ospray; cd ospray && git checkout v1.8.5 && mkdir -p build && cd build && cmake -DCMAKE_C_COMPILER='+ctx.env["CC"]+' -DCMAKE_CXX_COMPILER='+ctx.env["CXX"]+' -DEMBREE_DIR=${PREFIX} -DCMAKE_INSTALL_PREFIX=${PREFIX} .. && make -j '+str(p)+' && make install', # build ospray
      "sed -i 's|int megamol_main(int argc, char\* argv\[\]) {|int main(int argc, char\* argv\[\]) {|g' ${SOURCE_DIR}/console/src/Console.cpp",  # reverse eventual change to Console.cpp
      'cd ${SOURCE_DIR} && mkdir -p build && cd build && \
      '+ctx.env["cmake"]+' -DCMAKE_INSTALL_PREFIX=${PREFIX} -DUSE_MPI=ON -DMPI_GUESS_LIBRARY_NAME= -DADIOS2_DIR=${DEPENDENCIES_DIR}/adios/install/lib/cmake/adios2 \
      -DBUILD_CINEEMATICCAMERA_PLUGIN=OFF -DBUILD_GEOMETRY_CALLS_PLUGIN=ON -DBUILD_GUI_PLUGIN=OFF -DBUILD_IMAGEVIEWER2_PLUGIN=OFF -DBUILD_IMAGE_CALLS_PLUGIN=OFF -DBUILD_INFOVIS_PLUGIN=OFF -DBUILD_MDAO2_PLUGIN=OFF \
      -DBUILD_MMSTD_MOLDYN_PLUGIN=OFF -DBUILD_MMSTD_TRISOUP_PLUGIN=OFF -DBUILD_MMSTD_VOLUME_PLUGIN=OFF -DBUILD_PROTEIN_PLUGIN=OFF -DBUILD_PWDEMOS_PLUGIN=OFF -DBUILD_MWDEMOS_PLUGIN=OFF -DBUILD_PROTEIN_CALLS_PLUGIN=OFF -DUSE_LIBUNWIND=OFF \
      -DUSE_EXTERNAL_ADIOS=ON \
      -DBUILD_OSPRAY_PLUGIN=ON -DBUILD_PBS_PLUGIN=ON -Dospray_DIR=${PREFIX} \
      CXXFLAGS="-Wno-int-in-bool-context -Wno-endif-labels -Wno-reorder -Wno-pedantic" .. \
      && make -j '+str(p)+' && make install',
      "cp ${SOURCE_DIR}/console/src/Console.cpp ${SOURCE_DIR}/console/src/Console.cpp.backup && \
      sed -i 's|int main(int argc, char\* argv\[\]) {|int megamol_main(int argc, char\* argv\[\]) {|g' ${SOURCE_DIR}/console/src/Console.cpp",  # replace main() by megamol_main in Console.cpp
      'mkdir -p ${PREFIX}/include && printf "#pragma once\nint megamol_main(int argc, char* argv[]);" > ${PREFIX}/include/Console.h',
      '!cd ${SOURCE_DIR}/build/console && echo "linking mmconsole will fail, but that is expected and okay!" && make'   # build Console.cpp.o with megamol_main() instead of main(), this fails when linking mmconsole, but that is okay
    ])
    self.number_output_lines = 12050

    res = super(MegaMol, self).check(ctx)
    self.check_required(res[0], ctx)
    ctx.Log("Note, the previous check failed, because mmconsole cannot be linked. This is expected and no error (Just repeating in case you didn't read the full log.)\n")
      
    # define custom linker flags, this needs the src tree present
    cwd = os.getcwd()
    opendihu_home = os.getcwd()
    ctx.Log("opendihu_home: {}\n".format(opendihu_home))
    
    megamol_src_path = os.path.join(opendihu_home, "dependencies/megamol/src")
    if os.path.exists(megamol_src_path):
      entries = os.listdir(megamol_src_path)
      if len(entries) == 1:
        if (os.path.isdir(os.path.join(megamol_src_path, entries[0]))):
          megamol_src_path = os.path.join(megamol_src_path, entries[0])
      else:
        ctx.Message("\nError! There are multiple folders in {}: \n{}\nPlease remove all but the one containing MegaMol!\n".format(megamol_src_path,entries))
        megamol_src_path = os.path.join(megamol_src_path, entries[0])
      
      ctx.Log("megamol_src_path: {}\n".format(megamol_src_path))
      
      megamol_link_path = os.path.join(megamol_src_path, "build/console/CMakeFiles/mmconsole.dir/link.txt")
      ctx.Log("megamol_link_path: {}\n".format(megamol_link_path))
      try:
        with open(megamol_link_path, "r") as f:
          linker_flags = f.readline()
          linker_flags = linker_flags[linker_flags.find(" "):]
          #print(linker_flags)

          new_linker_flags = []
          for item in linker_flags.split():
            if item == "../mmconsole" or item == "-o":
              continue
              
            if item[0] != "-" and item[0] != "/":
              item = os.path.join(os.path.join(megamol_src_path, "build/console/"), item)
            new_linker_flags.append(item)
          
          ctx.Log("new_linker_flags: {}\n".format(new_linker_flags))
          self.link_flags = new_linker_flags
          
        ctx.Log("The next check should succeed if MegaMol was installed correctly.\n")
      except:
        ctx.Log("Could not open {}".format(megamol_link_path))
      
      # run again with new linker flags and modified Console.cpp
      #self.set_build_handler(None)

      res = super(MegaMol, self).check(ctx)
      self.check_required(res[0], ctx)
    else:
      ctx.Log("megamol src path {} does not exist\n".format(megamol_src_path))
      
    ctx.Result(res[0])
    return res[0]
