import sys, os, multiprocessing, subprocess, socket
from .Package import Package

#
# Python requires bzip2 (for matplotlib to work)
#
class Python(Package):

  def __init__(self, **kwargs):
    defaults = {
      #'download_url': 'https://www.python.org/ftp/python/2.7.14/Python-2.7.14.tgz',
      #'download_url': 'https://www.python.org/ftp/python/2.7.12/Python-2.7.12.tgz',
      #'download_url': 'https://www.python.org/ftp/python/3.4.6/Python-3.4.6.tgz'
      #'download_url': 'https://www.python.org/ftp/python/3.5.3/Python-3.5.3.tgz'
      'download_url': 'https://www.python.org/ftp/python/3.6.5/Python-3.6.5.tgz'
      #'download_url': 'https://github.com/stackless-dev/stackless/archive/3.6-slp.zip'  # stackless-python
        
    }
    defaults.update(kwargs)
    super(Python, self).__init__(**defaults)
    self.ext = '.cpp'
    self.sub_dirs = [
        ('include', 'lib'),
    #    ('include/mysql', 'lib64'),
    ]
    #self.headers = ['mysql.h']
    #self.libs = ['mysqlclient']
    #self.extra_libs = ['lapack', 'blas']
    self.check_text = r'''
      #include <stdlib.h>
      #include <stdio.h>
      #include <Python.h>
      int main(int argc, char* argv[])
      {
        Py_Initialize();
        return EXIT_SUCCESS;
      }
'''

  def check(self, ctx):
    env = ctx.env
    
    if socket.gethostname() != 'cmcs09':
      # on normal host
        
      self.libs = ["python3.6m"]
      self.headers = ["Python.h"]
  
      # check configuration of gcc
      #gcc_config = subprocess.check_output(["gcc", "-v"], stderr=subprocess.STDOUT).decode("utf-8")
     
      # extract and output gcc version
      #pos1 = gcc_config.find("gcc version")
      #pos2 = pos1+11+gcc_config[pos1+11:].find(" ")
      #pos3 = pos2+gcc_config[pos2+1:].find(" ")
      #print("GCC version: {}".format(gcc_config[pos2+1:pos3+1]))
      gcc_config = ""
   
      # if gcc was compiled such that -fuse-linker-plugin is available, compile with optimizations, (now disabled because it takes long and sometimes fails)
      if "--enable-plugin" in gcc_config and False:
        #print("gcc has --enable-plugin, compile python with optimizations")
        self.set_build_handler([
          'mkdir -p ${PREFIX}',
          'cd ${SOURCE_DIR} && chmod +x ./configure && CC="'+env['CC']+'" CXX="'+env['CXX']+'" ./configure --enable-shared --enable-optimizations --prefix=${PREFIX} \
            LDFLAGS="-Wl,--rpath=${PREFIX}/lib -L${DEPENDENCIES_DIR}/bzip2/install/lib" \
            CFLAGS="-I${DEPENDENCIES_DIR}/bzip2/install/include" \
            && make && make install',
          '$export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${PREFIX}/lib',
          'cd ${PREFIX}/include && echo "#define PYTHON_HOME_DIRECTORY \\"${PREFIX}\\"\n" > python_home.h',
          "sed -i 's/#define clock_t long/\/\/#define clock_t long/g' ${PREFIX}/include/python3.6m/pyconfig.h",   # this is needed for pgi on argon
          "sed -i 's/#define gid_t int/\/\/#define gid_t int/g' ${PREFIX}/include/python3.6m/pyconfig.h",
          "sed -i 's/#define uid_t int/\/\/#define uid_t int/g' ${PREFIX}/include/python3.6m/pyconfig.h",
        ])
        self.number_output_lines = 9823
      else:       
        #if env["CC"] == "gcc":
        #  print("gcc has no --enable-plugin, compile python without optimizations")
        self.set_build_handler([
          'mkdir -p ${PREFIX}',
          'cd ${SOURCE_DIR} && chmod +x ./configure && CC="'+env['CC']+'" CXX="'+env['CXX']+'" ./configure --enable-shared --prefix=${PREFIX} \
            LDFLAGS="-Wl,--rpath=${PREFIX}/lib -L${DEPENDENCIES_DIR}/bzip2/install/lib" \
            CFLAGS="-I${DEPENDENCIES_DIR}/bzip2/install/include" \
            && make && make install',
          '$export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${PREFIX}/lib',
          'cd ${PREFIX}/include && echo "#define PYTHON_HOME_DIRECTORY \\"${PREFIX}\\"\n" > python_home.h',
          "sed -i 's/#define clock_t long/\/\/#define clock_t long/g' ${PREFIX}/include/python3.6m/pyconfig.h",  # this is needed for pgi on argon
          "sed -i 's/#define gid_t int/\/\/#define gid_t int/g' ${PREFIX}/include/python3.6m/pyconfig.h",
          "sed -i 's/#define uid_t int/\/\/#define uid_t int/g' ${PREFIX}/include/python3.6m/pyconfig.h",
        ])
    
        self.number_output_lines = 7082
        
      # build stackless-python
      if False:
        self.set_build_handler([
          'mkdir -p ${PREFIX}',
          'cd ${SOURCE_DIR} && chmod +x ./configure && ./configure --enable-shared --prefix=${PREFIX} \
            LDFLAGS="-Wl,--rpath=${PREFIX}/lib -L${DEPENDENCIES_DIR}/bzip2/install/lib" \
            CFLAGS="-I${DEPENDENCIES_DIR}/bzip2/install/include" \
            \
            && export l1=$(grep "CC=" Makefile -n | head -n 1 | cut -d: -f1) \
            && export l2=$(grep "CXX=" Makefile -n | head -n 1 | cut -d: -f1) \
            && sed -i "$l1s/.*/CC=		gcc -pthread -fPIC/" Makefile \
            && sed -i "$l2s/.*/CXX=		g++ -pthread -fPIC/" Makefile \
            && make && make install',
          '$export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${PREFIX}/lib',
          'cd ${PREFIX}/include && echo "#define PYTHON_HOME_DIRECTORY \\"${PREFIX}\\"\n" > python_home.h',
        ])

    
    # on cmcs09, use PGI compiler and hardcoded paths
    if socket.gethostname() == 'cmcs09':
      # pgi compiler

      # Setup the build handler.
      self.libs = ["python3.6m"]
      self.link_flags = ["-L/lib/x86_64-linux-gnu","-ldl","-lutil"]
      self.headers = ["Python.h"]

      #print("ATTENTION! Manipulating ENVIRONMENT VARIABLES to install Python 'PGI-unaware'!")
      print("cmcs09: Using pre-compiled Python module.")
      # get all environment variables which are needed to be reset for a moment.
      # store them, set them as needed for now, do what is necessary, reset them later to continue with PGI
#          env_path = os.getenv("PATH","none")
#          if env_path =="none":
#            print("ERROR in scons package Python.py: PATH does not exist")
#          else:
#            previousPATH=env_path
#            os.environ["PATH"] = "/home/kraemer/perl5/perlbrew/bin:/home/kraemer/perl5/perlbrew/perls/perl-5.10.1/bin:/afs/.mathe/home/cmcs/share/environment-modules/Packages/gcc/7.2.0/bin:/afs/.mathematik.uni-stuttgart.de/home/cmcs/share/environment-modules/Packages/cmake/3.5.0/bin:/afs/.mathematik.uni-stuttgart.de/home/cmcs/share/environment-modules/Modules/3.2.10/bin:/usr/local/bin:/usr/bin:/bin:/usr/local/games:/usr/games:/usr/local/cuda-9.0/bin:/opt/thinlinc/bin"
#            print("{}".format(previousPATH))
#            print("{}".format(env_path))
#          env_libPath = os.getenv("LIBRARY_PATH","none")
#          if env_libPath=="none":
#            print("ERROR in scons package Python.py: LIBRARY_PATH does not exist")
#          else:
#            previousLIBRARY_PATH=env_libPath
#            os.environ["LIBRARY_PATH"] = "/afs/.mathe/home/cmcs/share/environment-modules/Packages/gcc/7.2.0/lib64"
#            print("{}".format(previousLIBRARY_PATH))
#            print("{}".format(env_libPath))
#          env_ldLibPath = os.getenv("LD_LIBRARY_PATH","none")
#          if env_ldLibPath=="none":
#            print("ERROR in scons package Python.py: LD_LIBRARY_PATH does not exist")
#          else:
#            previousLD_LIBRARY_PATH=env_ldLibPath
#            os.environ["LIBRARY_PATH"] = "/afs/.mathe/home/cmcs/share/environment-modules/Packages/gcc/7.2.0/lib64:/usr/local/cuda-10.0/lib64"
#            print("{}".format(previousLD_LIBRARY_PATH))
#            print("{}".format(env_ldLibPath))


      #-fprofile-arcs -ftest-coverage -fprofile-generate 
      self.set_build_handler([
        'mkdir -p ${PREFIX}',
        'cd ${SOURCE_DIR} && chmod +x ./configure && ./configure --enable-shared --prefix=${PREFIX} \
          LDFLAGS="-Wl,--rpath=${PREFIX}/lib -L${DEPENDENCIES_DIR}/bzip2/install/lib" \
          CFLAGS="-I${DEPENDENCIES_DIR}/bzip2/install/include " CXX=g++ CC=gcc \
          && make && make install',
        '$export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${PREFIX}/lib',
        'cd ${PREFIX}/include && echo "#define PYTHON_HOME_DIRECTORY \\"${PREFIX}\\"\n" > python_home.h',
      ])
      self.number_output_lines = 7082
      # reset the previously changed environment variables to their original PGI-case values.

    ctx.Message('Checking for Python ...        ')
    self.check_options(env)

    res = super(Python, self).check(ctx)

    self.check_required(res[0], ctx)
    ctx.Result(res[0])
    return res[0]
