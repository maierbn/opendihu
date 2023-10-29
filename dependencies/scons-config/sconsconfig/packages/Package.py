import os, sys, copy, shutil, subprocess, shlex
import time
import datetime
import socket
from threading import Thread
import sconsconfig.utils as utils
from sconsconfig.utils import conv
from SCons.Variables import BoolVariable

## Backup a list of environment variables.
# @param[in] env An SCons environment.
# @param[in] names A list of environment names to backup.
# @returns A dictionary of backed up names and values.
def env_backup(env, names):
  names = conv.to_iter(names)
  bkp = {}
  for n in names:
    if n not in env:
      bkp[n] = None
    else:
      bkp[n] = env[n]
  return bkp

## Backup the existing environment and update with provided keywords.
# @see env_backup
# @param[in] env An SCons environment.
# @returns A dictionary of backed up names and values.
def env_setup(env, **kw):
  bkp = env_backup(env, list(kw.keys()))
  env.Replace(**kw)
  return bkp

## Restore a set of previously backed up environment macros.
# @see env_backup
# @param[in] env An SCons environment.
# @param[in] bkp A dictionary of backed up key/values returned by env_backup.
def env_restore(env, bkp):
  '''
  '''
  for n, v in bkp.items():
    if v is None:
      del env[n]
    else:
      env[n] = v

## Base class for all configuration packages.
class Package(object):

  ## Default include/library sub-directories to search.
  DEFAULT_SUB_DIRS = [('include', 'lib'), ('include', 'lib64')]
  one_shot_options = []  # options that should only be applied once per call to scons, e.g. REDOWNLOAD should only happen for the first target, not again for further targets
  compilers_checked = False  # if the C and CXX compilers have been checked

  ##
  # @param[in] required Boolean indicating whether the configuration should fail if this package
  #           cannot be found.
  # @param[in] download_url Location to download the package if required.
  def __init__(self, required=True, download_url=''):
    self.name = self.__class__.__name__
    self.required = required
    self.found = False
    self.headers = []         # header files of this package, that need to be included by the main program, have to be in directory "include"
    self.libs = []            # libraries of this package, that need to be linked against by the main program, have to be in directory "lib" or "lib64"
    self.extra_libs = []      # extra i.e. extern libraries that are required for the package
    self.sources = []         # sources that need to be compiled with every program, have to be in subdirectory "src", note, these are not the sources needed to compile the library
    self.sub_dirs = self.DEFAULT_SUB_DIRS
    self.check_text = ''      # a test program that checks if the package is included successfully
    self.ext = '.c'           # extension of the check program, this determines if it is compiled with g++ or gcc
    self.options = []         # scons options that should be available to the user via the config script (default includes *_DIR, *_DOWNLOAD etc.)
    self.test_names = ['Check' + self.name]
    self.custom_tests = {'Check' + self.name: self.check}
    self.auto_add_libs=True
    self.run=True
    self.download_url = download_url    # the url from where to download the package
    self.build_handlers = {}            # the build handlers, use set_build_handler to set
    self.build_flags = ''               # additional C flags to use for compiling the test program
    self.number_output_lines = False    # number of output lines in typical compilation output (False to disable), used for monitoring compilation progress on stdout
    self.static = False                 # if the compiled test program is a static library
    self.set_rpath = True               # if the rpath in the linker should also be set (dynamic linkage)
    self.link_flags = None              # additional linker flags that can directly be set by the derived class
    self.last_build_log = ""            # after performing a build this contains the console output of the build process
    
    self.base_dir = None                # will be set to the base directory that contains "include" and "lib"
    self._used_inc_dirs = None
    self._used_libs = None


  def check_compilers(self, ctx):
    
    if Package.compilers_checked:
      return True
    
    ctx.Log("CC: {}, CXX: {}\n".format(ctx.env["CC"], ctx.env["CXX"]))
    for compiler in [ctx.env["CC"], ctx.env["CXX"]]:
      cmd = "{} --version".format(compiler)
      
      ctx.Log("Checking compiler {}\n".format(compiler))

      # Make a file to log stdout from the commands.
      stdout_log = open('stdout.log', 'w')
      try:
        subprocess.check_call(cmd, stdout=stdout_log, stderr=subprocess.STDOUT, shell=True)
        
        # get output
        with open('stdout.log') as f:
          output = f.read()
        ctx.Log("$"+cmd+"\n")
        ctx.Log(output+"\n")
      except:
        self.command_running = False
        stdout_log.close()
        with open('stdout.log') as f:
          output = f.read()
        ctx.Log("Command failed: \n"+output)
        
        # try again with "-V" instead of "--version" (for cray compiler)
        cmd = "{} -V".format(compiler)
        
        # Make a file to log stdout from the commands.
        stdout_log = open('stdout.log', 'w')
        try:
          subprocess.check_call(cmd, stdout=stdout_log, stderr=subprocess.STDOUT, shell=True)
          
          # get output
          with open('stdout.log') as f:
            output = f.read()
          ctx.Log("$"+cmd+"\n")
          ctx.Log(output+"\n")
        except:
          self.command_running = False
          stdout_log.close()
          with open('stdout.log') as f:
            output = f.read()
          ctx.Log("Command failed: \n"+output)
          
          sys.stdout.write('\nError: Compiler \"{}\" not found. Set the cc and CC variables appropriately.'.format(compiler))
          ctx.Log('Compiler \"{}\" not found.\n'.format(compiler))
          return False
      
      if output[0:3] == "gcc":
        major_no = None
        try:
          gcc_version = output[0:output.find("\n")]
          gcc_version = gcc_version[gcc_version.rfind(" ")+1:]
          major_no = (int)(gcc_version[0:gcc_version.find(".")])
          ctx.Log("GCC major version no: {}\n".format(major_no))
        except:
          pass
        if major_no and major_no < 7:
          print(("\n\nError! You are using GCC version {}, but at least version 7 is required.\nInstall a newer version of GCC or adjust the cc and CC variables in user-variables.scons.py to specify a newer compiler.\n".format(gcc_version)));
          quit()
    Package.compilers_checked = True
    return True
    
  ##
  # TODO: Make more general
  def include_directories(self):
    return os.path.join(self.base_dir, 'include')

  ##
  # TODO: This needs to be more general.
  def libraries(self):
    if self._used_libs:
      return os.path.join(self.base_dir, 'lib', 'lib' + self._used_libs[0] + '.so')
    else:
      return ''
      
  def determine_dependencies_dir(self, ctx):
    
    dependencies_dir_candidates = [
        os.path.join(os.getcwd(),'dependencies'), os.path.join(os.getcwd(),'../dependencies'), os.path.join(os.getcwd(),'../../dependencies')]
    
    # add directory from environment variable "OPENDIHU_HOME"
    if os.environ.get('OPENDIHU_HOME') is not None:
      if os.environ.get('OPENDIHU_HOME') != "":
        dependencies_dir_candidates = [os.path.join(os.environ.get('OPENDIHU_HOME'),'dependencies')] + dependencies_dir_candidates
    
    dependencies_dir = os.path.join(os.getcwd(),'dependencies')
    for directory in dependencies_dir_candidates:
      if directory is not None:
        if os.path.exists(directory):
          dependencies_dir = directory
          break
    
    self.check_text = self.check_text.replace('${DEPENDENCIES_DIR}', dependencies_dir)
    return dependencies_dir

  ## Run the configuration checks for this package.
  # @param[in,out] ctx The configuration context, retrieved from SCons.
  def check(self, ctx, **kwargs):
    ctx.Log('\n====================================\n')
    ctx.Log('Beginning check for %s\n'%self.name)
    env = ctx.env
    name = self.name
    libs = self.libs
    extra_libs = self.extra_libs
    sub_dirs = self.sub_dirs

    upp = name.upper()

    # check if compiler works
    if not self.check_compilers(ctx):
      return (False, 0)
    
    dependencies_dir = self.determine_dependencies_dir(ctx)
    
    #ctx.Log("ctx.env.items: "+str(ctx.env.items())+"\n")
    
    # Check if the user requested to download this package.
    if self.have_any_options(env, upp + '_DOWNLOAD', 'DOWNLOAD_ALL', upp + '_REBUILD'):

      # Perform the auto-management of this package.
      res = self.download_and_build_package(ctx)

      # Here, the package location is in <NAME>_DIR, either by option or by download and build
      if res[0]:
        value = self.get_option(env, upp + '_DIR')
        res = self.try_location(ctx, value, **kwargs)
        
        if res[0]:
          # If it all seemed to work, write a dummy file to indicate this package has been built.
          ctx.Log("Build was successful, touch file \"scons_build_success\".\n")
          
          success_file = open(os.path.join(env['current_source_dir'],'scons_build_success'), 'w')
          success_file.write('  ')
          success_file.close()  
        
        else:
          self._msg = '\n\nUnable to validate a %s installation at:\n %s\nInspect "config.log" to see what went wrong.\n'%(name, value)
          # ctx.Log(msg)
          # print(msg)
          # env.Exit(1)

    elif self.have_option(env, upp + '_DIR'):
      value = self.get_option(env, upp + '_DIR')
      ctx.Log('Found option %s = %s\n'%(upp + '_DIR', value))
      res = self.try_location(ctx, value, **kwargs)
      if not res[0]:
        self._msg = '\n\nUnable to validate a %s installation at:\n %s\nInspect "config.log" to see what went wrong.'%(name, value)
        # ctx.Log(msg)
        # print(msg)
        # env.Exit(1)

    elif self.have_any_options(env, upp + '_INC_DIR', upp + '_LIB_DIR', upp + '_LIBS'):
      inc_dirs = self.get_option(env, upp + '_INC_DIR').split(';')
      lib_dirs = self.get_option(env, upp + '_LIB_DIR').split(';')
      if self.have_option(env, upp + '_LIBS'):
        cur_libs = [list(map(env.File, self.get_option(env, upp + '_LIBS').split(';')))]
        cur_extra_libs = []
      else:
        cur_libs = libs
        cur_extra_libs = extra_libs
      ctx.Log('Found options:\n')
      if inc_dirs:
        kwargs['inc_dirs'] = inc_dirs
        ctx.Log(' %s = %s\n'%(upp + '_INC_DIR', str(inc_dirs)))
      if lib_dirs:
        kwargs['lib_dirs'] = lib_dirs
        ctx.Log(' %s = %s\n'%(upp + '_LIB_DIR', str(lib_dirs)))
      if self.have_option(env, upp + '_LIBS'):
        kwargs['libs'] = cur_libs
        kwargs['extra_libs'] = []
        ctx.Log(' %s = %s\n'%(upp + '_LIBS', str([l.path for l in cur_libs])))

      res = (True, 0)
      if inc_dirs:
        if not self.try_headers(ctx, inc_dirs):
          res = (False, 0)
      
      if self.set_rpath:
        bkp = env_setup(ctx.env,
                CPPPATH=ctx.env.get('CPPPATH', []) + inc_dirs,
                LIBPATH=ctx.env.get('LIBPATH', []) + lib_dirs,
                RPATH=ctx.env.get('RPATH', []) + lib_dirs)
      else:
        bkp = env_setup(ctx.env,
                CPPPATH=ctx.env.get('CPPPATH', []) + inc_dirs,
                LIBPATH=ctx.env.get('LIBPATH', []) + lib_dirs)
        
      if res[0]:
        res = self.try_libs(ctx, libs, **kwargs)
      if not res[0]:
        env_restore(ctx.env, bkp)
        self._msg = '\n\nUnable to validate a %s installation using:\n'%name
        if self.have_option(env, upp + '_INC_DIR'):
          self._msg += '  Header directories: %s\n'%str(inc_dirs)
        if self.have_option(env, upp + '_LIB_DIR'):
          self._msg += '  Library directories: %s\n'%str(lib_dirs)
        if self.have_option(env, upp + '_LIBS'):
          self._msg += '  Libraries: %s\n'%str([l.path for l in cur_libs])
        self._msg += 'Inspect "config.log" to see what went wrong.\n'
        # ctx.Log(msg)
        # print(msg)
        # env.Exit(1)

    else:
      ctx.Log('None of the options \"' + upp + '_DOWNLOAD\", \"DOWNLOAD_ALL\", \"' + upp + '_REBUILD\", \"' + upp + '_DIR\", \"' + upp + '_INC_DIR\", \"' + upp + '_LIB_DIR\", \"' + upp + '_LIBS\" found. Search in current location \".\".\n')
      self.base_dir = "."
      res = self.try_libs(ctx, libs, extra_libs, **kwargs)

      if not res[0]:
        common_dirs = [os.path.join(os.getcwd(),'dependencies'), os.path.join(os.getcwd(),'../dependencies'), os.path.join(os.getcwd(),'../../dependencies'),
                       '/usr', '/usr/lib/openmpi', '/usr/local', os.environ['HOME'], os.path.join(os.environ['HOME'], 'soft'), '/sw']
        ctx.Log('\n  Trying common locations: %s\n'%str(common_dirs))
        res = (0, '')
        
        # loop over common directories
        for cd in common_dirs:
          res = self.try_location_subtree(ctx, cd, **kwargs)
          if res[0]:
            break
            
    disable_checks = False
    if self.have_option(env, "DISABLE_CHECKS"):
      if ctx.env.get('DISABLE_CHECKS', []):
        ctx.Log('Disable checks because DISABLE_CHECKS is set')
        disable_checks = True
    if self.have_option(env, upp + '_DISABLE_CHECKS'):
      if ctx.env.get(upp + '_DISABLE_CHECKS', []):
        ctx.Log('Disable checks because '+upp + '_DISABLE_CHECKS is set')
        disable_checks = True
    if disable_checks:
      res = (True,'')
    return res

  def try_location_subtree(self, ctx, base_dirs, **kwargs):
    cd = base_dirs
    name = self.name
    upp = name.upper()
    env = ctx.env
    res = (False, None)
    
    ctx.Log('  ------------------------------\n')
    
    # check if path exists
    if not os.path.exists(cd):
      ctx.Log('%s does not exist.\n'%cd)
      return res
  
    ctx.Log('  Looking in %s\n'%cd)
    
    # try if this is the correct location
    res = self.try_location(ctx, cd, **kwargs)
    if res[0]:
      env[upp + '_DIR'] = cd
      return res
    
    # loop over subdirectories  
    walk_directories = False            
    for d in os.listdir(cd):
      ctx.Log("    Directory %s\n"%d)
      
      # descend if subdirectory is "install" or "build" or contains name of package
      if (d.lower() == "install" or d.lower() == "build" or d.lower().find(name.lower()) != -1) and os.path.isdir(os.path.join(cd, d)):
        d = os.path.join(cd, d)
        ctx.Log('  Descend to %s\n'%d)
        
        res = self.try_location_subtree(ctx, d, **kwargs)
        if res[0]:
          return res
    return res

  def check_required(self, result, ctx=None):
    name = self.name
    upp = name.upper()

    # Stash the result.
    self.found = bool(result)

    # Failed specified location?
    if not result and hasattr(self, '_msg') and self.required:
      ctx.Log(self._msg)
      print((self._msg))
      sys.exit(1)

    # General failure?
    if not result and self.required:
      print('\n')
      print(('Unable to locate required package %s. You can do the following:'%name))
      print(('  - Define %s_DIR with a directory containing "include" and "lib" or "lib64" subdirectories.'%upp))
      print(('  - Define %s_INC_DIR, %s_LIB_DIR and %s_LIBS'%(upp, upp, upp)))
      print(('  - Set %s_DOWNLOAD to True to automatically download and install the package.'%upp))
      print(('    If you additionally set %s_REDOWNLOAD to True it forces a fresh download '%upp))
      print('    if it was already done earlier.')
      print('If you already did that, something went wrong. Inspect "config.log" to see what went wrong.\n\nSometimes it helps to delete the cache and try again:\nrm -rf .scon*')
      sys.exit(1)

    # If the package is not required but was found anyway, add a preprocessor
    # flag indicating such.
    elif result:
      if ctx:
        ctx.env.AppendUnique(CPPDEFINES=['HAVE_' + upp])

  def add_options(self, vars):
    name = self.name
    upp = name.upper()
    vars.Add(BoolVariable('DISABLE_CHECKS', help='Allow failure of checks for all packages.', default=False))
    vars.Add(BoolVariable('DISABLE_RUN', help='Do not attempt to run executables for checks, only compile, for all packages.', default=False))
    vars.Add(upp + '_DIR', help='Location of %s.'%name)
    vars.Add(upp + '_INC_DIR', help='Location of %s header files.'%name)
    vars.Add(upp + '_LIB_DIR', help='Location of %s libraries.'%name)
    vars.Add(upp + '_LIBS', help='%s libraries.'%name)
    vars.Add(upp + '_DISABLE_CHECKS', help='Allow failure of check for %s.'%name)
    vars.Add(upp + '_DISABLE_RUN', help='Do not attempt to run executables for checks for %s.'%name)
    if self.download_url:
      vars.Add(BoolVariable(upp + '_DOWNLOAD', help='Download and use a local copy of %s.'%name, default=False))
      vars.Add(BoolVariable(upp + '_REDOWNLOAD', help='Force update of previously downloaded copy of %s.'%name, default=False))
    vars.Add(BoolVariable(upp + '_REBUILD', help='Force new build of previously downloaded copy of %s, even if it was installed successfully.'%name, default=False))
    vars.Add(BoolVariable('MPI_IGNORE_MPICC', help='Disable retrieving mpi compile information from mpicc --showme.', default=False))
    vars.Add(BoolVariable('MPI_DEBUG', help='Build MPI with debugging support for memcheck.', default=False))
    vars.Add(BoolVariable('PETSC_DEBUG', help='Build Petsc with debugging support.', default=False))
    
    self.options.extend([upp + '_DIR', upp + '_INC_DIR', upp + '_LIB_DIR', upp + '_LIBS', upp + '_DOWNLOAD'])

  ## Set the build handler for an architecture and operating system. Pass in None to the handler
  # to clear the current handler.
  # @param[in] handler The handler string to use.
  # @param[in] sys_id The kind of system to set for.
  def set_build_handler(self, handler, sys_id=None):
    if handler is None and sys_id in self.build_handlers:
      del self.build_handlers[sys_id]
    else:
      self.build_handlers[sys_id] = handler

  def get_build_handler(self):
    import platform
    os = platform.system()
    arch = platform.architecture()[0]
    sys_id = os + '_' + arch
    if sys_id in self.build_handlers:
      return self.build_handlers[sys_id]
    if os in self.build_handlers:
      return self.build_handlers[os]
    if arch in self.build_handlers:
      return self.build_handlers[arch]
    return self.build_handlers.get(None, None)

  def download_and_build_package(self, ctx):
    #sys.stdout.write('\n')
    
    # Are we forcing this?
    upp = self.name.upper()
    force_redownload = self.have_option(ctx.env, upp + '_REDOWNLOAD')
    force_rebuild = self.have_option(ctx.env, upp + '_REBUILD')
    
    if force_redownload:
      ctx.Log("Option force redownload is set.\n")
      ctx.Log("Already performed one shot options that are prevented a second time: "+str(Package.one_shot_options)+"\n")
      if upp + '_REDOWNLOAD' in Package.one_shot_options:
        ctx.Log("Redownload is a one shot option and was performed earlier, disable it this time.\n")
        Package.one_shot_options.remove(upp + '_REDOWNLOAD')
        force_redownload = False
      else:
        Package.one_shot_options.append(upp + '_REDOWNLOAD')
        
    if force_rebuild:
      ctx.Log("Option force rebuild is set.\n")
      ctx.Log("one shot options: "+str(Package.one_shot_options)+"\n")
      if upp + '_REBUILD' in Package.one_shot_options:
        ctx.Log("Rebuild is a one shot option and was performed earlier, disable it this time.\n")
        Package.one_shot_options.remove(upp + '_REBUILD')
        force_rebuild = False
      else:
        Package.one_shot_options.append(upp + '_REBUILD')
        

    #ctx.Log("ctx.env.items: "+str(ctx.env.items()))
    #ctx.Log("ctx.env.items['ENV']: "+str(ctx.env['ENV'].items()))
    
    dependencies_dir = self.determine_dependencies_dir(ctx)
    base_dir = os.path.join(dependencies_dir,self.name.lower())


    ctx.Log("  dependencies:["+dependencies_dir+"]\n")
    ctx.Log("  base_dir:    ["+base_dir+"] (where to download the archive to)\n")
    
    # Create the source directory if it does not already exist.
    if not os.path.exists(base_dir):
      os.makedirs(base_dir)
    ctx.Log("Downloading into " + base_dir + "\n")

    # Setup the filename and build directory name and destination directory.
    if self.download_url == "" or self.download_url is None:
      filename = ""
    else:
      filename = self.download_url[self.download_url.rfind('/') + 1:]
    unpack_dir = "src"
    install_dir = os.path.abspath(os.path.join(base_dir, "install"))
    
    # Set the directory location.
    ctx.env[self.name.upper() + '_DIR'] = install_dir
    
    ctx.Log("Building into " + install_dir + "\n")

    # Change to the source directory.
    old_dir = os.getcwd()
    os.chdir(base_dir)

    ctx.Log("  filename:    ["+filename+"]\n")
    ctx.Log("  old_dir:     ["+old_dir+"]\n")
    ctx.Log("  install_dir: ["+install_dir+"] (where to build)\n")
    ctx.Log("  unpack_dir:  ["+unpack_dir+"] (where to unpack)\n")

    # Download if the file is not already available.
    if filename != "":
      if not os.path.exists(filename) or force_redownload:
        if not self.download_package(ctx, filename):
          os.chdir(old_dir)
          return (0, '')

    # Unpack if there is not already a build directory by the same name.
    if filename != "":
      if not os.path.exists(unpack_dir) or force_redownload:
        if not self.unpack_package(ctx, filename, unpack_dir):
          os.chdir(old_dir)
          return (0, '')

    # Move into the build directory. Most archives will place themselves
    # in a single directory which we should then move into.
    try:
      os.chdir(unpack_dir)
    except:
      pass
    entries = os.listdir('.')
    if len(entries) == 1:
      if (os.path.isdir(entries[0])):
        os.chdir(entries[0])
    source_dir = os.getcwd()
    ctx.Log("  source_dir:  ["+source_dir+"] (where the unpacked sources are)\n\n")
    ctx.env['current_source_dir'] = source_dir

    ctx.Log(" Only build package if at least on of the following three conditions is True:\n")
    ctx.Log(" force_redownload: "+str(force_redownload)+", force_rebuild: "+str(force_rebuild)+", not success: "+str(not os.path.exists('scons_build_success'))+"\n")
    ctx.Log(" (\"success\" means that the file \"scons_build_success\" exists in \""+source_dir+"\")\n")
    
    # Build the package.
    if (not os.path.exists('scons_build_success')) or force_redownload or force_rebuild:
      ctx.Log(" => Build package\n\n")
      if not self.build_package(ctx, install_dir, source_dir, dependencies_dir):
        os.chdir(old_dir)
        return (0, '')

    ctx.Log(' Configuring with downloaded package ... \n')
    os.chdir(old_dir)

    return (1, '')

  def download_package(self, ctx, filename):
    sys.stdout.write('\n  Downloading ... ')
    sys.stdout.flush()

    if os.path.exists(filename):
      os.remove(filename)

    ctx.Log("Downloading file from " + self.download_url + "\n")
    try:
      # normal download with python3
      import urllib.request, urllib.parse, urllib.error
      
    except Exception as e:
      # try with python2
      ctx.Log("Failed to download file: Could not import urllib.request\n")
      ctx.Log('Warning! You are running scons with python2.7, but python3 is required.\nPlease use the scons installation under\ndependencies/scons/scons.py\n\nNow trying to use libraries for python2.7 ... ')
      try:
        import urllib
        urllib.urlretrieve(self.download_url, filename)
        sys.stdout.write('done.\n')
        return True
      
      except Exception as e:
        ctx.Log("Failed to download file with python2.7.\nPlease use the scons installation under\ndependencies/scons/scons.py\n")
        return False
        
    try:
      urllib.request.urlretrieve(self.download_url, filename)
      sys.stdout.write('done.\n')
      return True
      
    except Exception as e:
      sys.stdout.write('failed.\n')
      print(e)
      ctx.Log("Failed to download file - retry in 5s\n")
      time.sleep(5)
      
      try:
        urllib.request.urlretrieve(self.download_url, filename)
        sys.stdout.write('done.\n')
        return True
      except Exception as e:
        sys.stdout.write('failed.\n')
        print(e)
        ctx.Log("Failed to download file again\n")
        return False
      
      return False

  def unpack_package(self, ctx, filename, unpack_dir):
    sys.stdout.write('  Extracting ... ')
    sys.stdout.flush()

    if os.path.exists(unpack_dir):
      shutil.rmtree(unpack_dir)

    # TODO: DRY
    ctx.Log("Extracting contents of " + filename + "\n")
    if os.path.splitext(filename)[1] == '.zip':
      ctx.Log("Using zip\n")
      try:
        import zipfile
        zf = zipfile.ZipFile(filename)
        zf.extractall(unpack_dir)
        zf.close()
      except:
        shutil.rmtree(unpack_dir, True)
        sys.stdout.write('failed.\n')
        ctx.Log("Failed to extract file\n")
        return False
    elif os.path.splitext(filename)[1] == '.whl':
      ctx.Log("Python wheel detected - do not extract\n")
      # create unpack_dir
      os.makedirs(unpack_dir)
      return True
    else:
      ctx.Log("Using tar\n")
      try:
        import tarfile
        tf = tarfile.open(filename)
        tf.extractall(unpack_dir)
        
        # get name of extracted directory
        entries = os.listdir(unpack_dir)
        ctx.Log("top-level files: {}".format(entries))
      except:
        shutil.rmtree(unpack_dir, True)
        try:
          
          filename_base = os.path.splitext(filename)[0]
          
          if ".tar" in filename_base:
            filename_base = os.path.splitext(filename_base)[0]
          
          cmd = "pwd && ls -al && tar -xf "+filename+" && ls -al"
          ctx.Log("tar package failed, trying to use system tar\n")
          ctx.Log(cmd+"\n")
          
          # Make a file to log stdout from the commands.
          stdout_log = open('stdout.log', 'w')
          subprocess.check_call(cmd, stdout=stdout_log, stderr=subprocess.STDOUT, shell=True)
          
          ctx.Log(stdout_log)
          try:
            os.rename(filename_base, unpack_dir)
          except Exception as e:
            ctx.Log("tar succeeded but failed to rename {} to {}: {}".format(filename_base, unpack_dir, str(e)))
          
        except Exception as e:
          shutil.rmtree(unpack_dir, True)          
          sys.stdout.write('failed. '+str(e)+'\n')
          ctx.Log("Failed to extract file\n")
          return False

    # If there is a patch, try to patch code.
    patch = os.path.join(utils.get_data_prefix(), 'patches', self.name.lower() + '.patch')
    if os.path.exists(patch):
      ctx.Log('Trying to apply patch.')
      try:
        utils.apply_patch(unpack_dir, patch)
      except:
        shutil.rmtree(unpack_dir, True)
        ctx.Log('failed to apply patch\n')
        return False

    sys.stdout.write('done.\n')
    return True

  def monitor_progress(self, number_of_lines):
    if not number_of_lines:
      return
    
    while self.command_running:
      time.sleep(1)
  
      if os.path.isfile('stdout.log'):
        try:
          with open('stdout.log') as f:
            output = f.read()

        # if output contains non-asci characters, try to open in binary mode
        except UnicodeDecodeError:
          try:
            with open('stdout.log', 'rb') as f:
              output = f.read()
          except:
            output = ''
            pass
        n = output.count('\n')
      else:
        continue
        
      p = int(float(n) / number_of_lines * 100.0)
      if p > 100.0:
        sys.stdout.write(str(p)+"% (miscounted, sorry)"+"\b"*(len(str(p))+21))
      sys.stdout.write(str(p)+"%"+"\b"*(len(str(p))+1))
      sys.stdout.flush()

  def build_package(self, ctx, install_dir, source_dir, dependencies_dir):
    sys.stdout.write('  Building package {}, this could take a while ... \n'.format(self.name))
      
    sys.stdout.flush()
    ctx.Log("Building package in " + install_dir + "\n")

    # Remove any existing file used to indicate successful builds.
    import os
    if os.path.exists('scons_build_success'):
      os.remove('scons_build_success')

    # Remove the installation directory.
    try:
      if os.path.exists(install_dir):
        if os.path.islink(install_dir):
          os.unlink(install_dir)
        else:
          shutil.rmtree(install_dir)
    except:
      ctx.Log("Failed to remove installation directory {}".format(install_dir))

    # Hunt down the correct build handler.
    handler = self.get_build_handler()
    if handler is None:
      sys.stdout.write('failed.\n Inspect the log file "config.log" for more information.\n')
      ctx.Log("Failed to locate build handler\n")
      return False


    # Process each command in turn.
    for cmd in handler:

      # Make a file to log stdout from the commands.
      stdout_log = open('stdout.log', 'w')
      
      # It's possible to have a tuple, indicating a function and arguments.
      if isinstance(cmd, tuple):
        ctx.Log("Command is a Python function\n")
        func = cmd[0]
        args = cmd[1:]

        # Perform substitutions.
        args = [ctx.env.subst(a.replace('${PREFIX}', install_dir)) for a in args]
        args = [ctx.env.subst(a.replace('${SOURCE_DIR}', source_dir)) for a in args]
        args = [ctx.env.subst(a.replace('${DEPENDENCIES_DIR}', dependencies_dir)) for a in args]

        # Call the function.
        func(*args)

      else:

        # If the first character in a command is an "!", then it means we allow
        # errors from this command.
        allow_errors = False
        if cmd[0] == '!':
          allow_errors = True
          cmd = cmd[1:]
          
        # If the first character in a command is a "$", then it means environment variables are not substituted
        substitute_environment_variables = True
        if cmd[0] == '$':
          substitute_environment_variables = False
          cmd = cmd[1:]

        # Perform substitutions.
        cmd = cmd.replace('${PREFIX}', install_dir)
        cmd = cmd.replace('${SOURCE_DIR}', source_dir)
        cmd = cmd.replace('${DEPENDENCIES_DIR}', dependencies_dir)
        path_environment_variable = ""
        if os.environ.get("PATH") is not None:
          path_environment_variable = os.environ.get("PATH")
        cmd = cmd.replace('${PATH}', path_environment_variable)
        
        
        #ctx.Log("MPI_DIR in env: {}\n".format("MPI_DIR" in ctx.env))
        #ctx.Log("env[MPI_DIR]={}, replace in cmd\n".format(ctx.env["MPI_DIR"]))
        if "MPI_DIR" in ctx.env:
          cmd = cmd.replace('${MPI_DIR}', ctx.env["MPI_DIR"])
        
        if substitute_environment_variables:
          cmd = ctx.env.subst(cmd)
          
        sys.stdout.write("    $ "+cmd+"  ")
        sys.stdout.flush()
    
        ctx.Log("  $ "+cmd+"\n")

        output = ""
        try:
          self.command_running = True
          t = Thread(target=self.monitor_progress, args=(self.number_output_lines,))
          t.start()
          
          subprocess.check_call(cmd, stdout=stdout_log, stderr=subprocess.STDOUT, shell=True)
          
          self.command_running = False
          t.join()
          sys.stdout.write("     \n")
          sys.stdout.flush()
    
          # get output
          if os.path.exists('stdout.log'):
            with open('stdout.log') as f:
              output = f.read()    
            stdout_log.close()
            os.remove('stdout.log')
          self.last_build_log = output
          ctx.Log(output+"\n")
          
        except:
          self.command_running = False
          stdout_log.close()
          if os.path.exists('stdout.log'):
            with open('stdout.log') as f:
              output = f.read()
            stdout_log.close()
            os.remove('stdout.log')
          
          if not allow_errors:
            sys.stdout.write('failed.\n')
            ctx.Log("Command failed: \n"+output)
            self.last_build_log = output
            return False
          else:
            ctx.Log("Command failed (but allowed): \n"+output)
            self.last_build_log = output

    sys.stdout.write('  done.\n')
    return True

  def env_setup_libs(self, ctx, libs):
    defaults = {
      'prepend': True,
      'libraries': libs,
    }
    env = ctx.env

    # If we were given a dictionary update our defaults.
    if len(libs) and isinstance(libs[0], dict):
      defaults.update(libs[0])
      defaults['libraries'] = conv.to_iter(defaults['libraries'])

    # Remove any pre-existing libraries.
    defaults['libraries'] = [d for d in defaults['libraries'] if d not in env.get('LIBS', [])]
    
    # Prepend or append?
    if defaults['prepend']:
      libs = defaults['libraries'] + conv.to_iter(env.get('LIBS', []))
    else:
      libs = conv.to_iter(env.get('LIBS', [])) + defaults['libraries']

    return env_setup(env, LIBS=libs)

  ## try to compile (self.run=0) or compile and run (self.run=1) the given code snippet in self.check_text
  # Returns (1,'program output') on success and (0,'') on failure
  def try_link(self, ctx, **kwargs):
    
    text = self.check_text+'\n'   # ensure that file ends with newline for extra picky cray compiler
    #text = self.check_text+'//'+"{:%Y/%m/%d %H:%m:%S}".format(datetime.datetime.now())+'\n'   # ensure that file ends with newline for extra picky cray compiler
    
    
    bkp = env_setup(ctx.env, **kwargs)
    ctx.env.PrependUnique(CCFLAGS = self.build_flags)
    
    # add sources directly to test program
    for source in self.sources:
      
      if ctx.env["BUILD_TYPE"] == "debug":
        ctx.Log("debug mode detected, using debug version of {}".format(source))
        object_filename = os.path.join(self.base_dir, os.path.join("src", os.path.splitext(source)[0]+".o"))
      else:
        ctx.Log("not debug mode detected, using release version of {}".format(source))
        object_filename = os.path.join(self.base_dir, os.path.join("src", os.path.splitext(source)[0]+".release.o"))
      
      if os.path.isfile(object_filename):
        ctx.env.AppendUnique(LINKFLAGS = object_filename)
      else:
        sys.stdout.write("Error: Could not find file \"{}\".".format(object_filename))
        ctx.Log("Error: Could not find file \"{}\".".format(object_filename))
        
      #with open(filename) as file:
      #  text += "\n"+file.read()
    
    if self.static:
      ctx.env.PrependUnique(CCFLAGS = '-static')
      ctx.env.PrependUnique(LINKFLAGS = '-static')
     
    if self.link_flags is not None:
      ctx.Log("  link_flags is set, use additional link flags: {}\n".format(self.link_flags))
      ctx.env.PrependUnique(LINKFLAGS = self.link_flags)
      
    # compile with C++14 for cpp test files
    if 'cpp' in self.ext:
      if os.environ.get("PE_ENV") != "CRAY":
        ctx.env.PrependUnique(CCFLAGS = "-std=c++14")
      
    #ctx.Log(ctx.env.Dump())
    if "LIBS" in ctx.env: ctx.Log("  LIBS:     "+str(ctx.env["LIBS"])+"\n")
    if "LINKFLAGS" in ctx.env: ctx.Log("  LINKFLAGS:"+str(ctx.env["LINKFLAGS"])+"\n")
    if "LIBPATH" in ctx.env: ctx.Log("  LIBPATH:  "+str(ctx.env["LIBPATH"])+"\n")
    if "CC" in ctx.env: ctx.Log("  CC:       "+str(ctx.env["CC"])+"\n")
    if "CXX" in ctx.env: ctx.Log("  CXX:      "+str(ctx.env["CXX"])+"\n")
    if "CCFLAGS" in ctx.env: 
      ccflags = ctx.env["CCFLAGS"]
      for i in ccflags:
        ctx.Log("  CCFLAGS: "+str(i)+"\n")
    ctx.Log("=============\n")
    #ctx.Log("  CCFLAGS:  "+str(ctx.env["CCFLAGS"])+"\n")    # cannot do str(..CCFLAGS..) when it is a tuple
       
    if ctx.env.get('DISABLE_RUN', []):
      ctx.Log('Do not run executable, only link because DISABLE_RUN is set.\n')
      self.run = False

    name = self.name
    upp = name.upper()
    if ctx.env.get(upp + '_DISABLE_RUN', []):
      ctx.Log('Do not run executable, only link because '+upp + '_DISABLE_RUN is set.\n')
      self.run = False
  
 
    # compile / run test program
    if self.run:
      if os.environ.get("PE_ENV") is not None or "pg" in ctx.env["CC"]:
         ctx.Log("Do not run test, only link\n")
         ctx.Log("Reason: either on compute cluster ($PE_ENV is set) or PGI compiler (env[\"CC\"] contains \"pg\")\n")

         # on cluster or with PGI compiler, do not run program, only link
         res = (ctx.TryLink(text, self.ext), '')    
      else:
        
          # compile and run test program normally
          res = ctx.TryRun(text, self.ext)
          
          # if this failed, prepend 'mpirun -n 1' and if that fails, 'srun -n 1'
          if not res[0]:
            ctx.Log("Running the test program failed, now try with prefix 'mpirun -n 1'\n")
            
            ok = ctx.TryLink(text, self.ext) 
            if(ok): 
              prog = ctx.lastTarget 
              pname = prog.get_internal_path() 
              output = ctx.sconf.confdir.File(os.path.basename(pname)+'.out') 
              node = ctx.env.Command(output, prog, [ [ 'mpirun', '-n', '1', pname, ">", "${TARGET}"] ]) 
              ok = ctx.sconf.BuildNodes(node) 
              if ok: 
                outputStr = output.get_contents()
                res = (1, outputStr)
                print("(If you are seeing MPI errors, this is because programs can only be run with `mpirun`. This is not a problem.)")
              else:
                ctx.Log("Running the test program with mpirun failed, now try with prefix 'srun -n 1'\n")
                node = ctx.env.Command(output, prog, [ [ 'srun', '-n', '1', pname, ">", "${TARGET}"] ]) 
                ok = ctx.sconf.BuildNodes(node) 
                if ok: 
                  outputStr = output.get_contents()
                  res = (1, outputStr)
                
    else:
      res = (ctx.TryLink(text, self.ext), '')
        
    # remove C++11 and C++14 flags, only for the test program
    if 'cpp' in self.ext:
      ccflags = ctx.env["CCFLAGS"]
      ccflags_new = []
      for entry in ccflags:
        if entry != "-std=c++14" and entry != "-std=c++11":
          ccflags_new.append(entry)
      #ctx.Log("recovered ccflags:")
      #ctx.Log(str(ccflags_new)+"\n")
      ctx.env.Replace(CCFLAGS = ccflags_new)
      
    if not res[0]:
      ctx.Log("Compile/Run failed.\n");
      
      if os.path.isfile(str(ctx.lastTarget)+".out"):
        with open(str(ctx.lastTarget)+".out", "rb") as f:
          try:
            output = f.read()
            if isinstance(output, str):
              ctx.Log("Program output: \""+output+"\"\n")
            else:
              ctx.Log("Program output: \""+output.decode('utf-8')+"\"\n")
          except:
              ctx.Log("(Could not load output)\n")
 
      disable_checks = False
      if ctx.env.get('DISABLE_CHECKS', []):
        ctx.Log('Disable checks because DISABLE_CHECKS is set.\n')
        disable_checks = True
      name = self.name
      upp = name.upper()
      if ctx.env.get(upp + '_DISABLE_CHECKS', []):
        ctx.Log('Disable checks because '+upp + '_DISABLE_CHECKS is set.\n')
        disable_checks = True
      if disable_checks:
        res = (True,'')
      else:
        env_restore(ctx.env, bkp)
    else:
      ctx.Log("Compile/Run succeeded.\n")
      if isinstance(res[1], str):
        ctx.Log("Program output: \""+res[1]+"\"\n")
      else:
        ctx.Log("Program output: \""+res[1].decode('utf-8')+"\"\n")
        
    return res

  def try_libs(self, ctx, libs, extra_libs=[], **kwargs):
    if not libs:
      libs = [[]]
    if not extra_libs:
      extra_libs = [[]]
      
    ctx.Log("try to compile with one of the following libs: "+str(libs)+"\n")
    if len(extra_libs) > 0:
      ctx.Log("also always try to include one of the following extra libs: "+str(extra_libs)+"\n")
    
    for l in libs:
      l = [l]
      ctx.Log("try library "+str(l)+"\n")
      
      l_bkp = self.env_setup_libs(ctx, l)
      for e in extra_libs:
        e = [e]
        # add extra lib
        linkflags = None
        if 'LINKFLAGS' in ctx.env:
          linkflags = ctx.env['LINKFLAGS']
        e_bkp = env_setup(ctx.env, LIBS=ctx.env.get('LIBS', []) + e, LINKFLAGS=linkflags)
        
        # try to link or run program
        res = self.try_link(ctx, **kwargs)
        
        # if succeeded
        if res[0]:
          if not self.auto_add_libs:
            env_restore(ctx.env, e_bkp)
            env_restore(ctx.env, l_bkp)
          self._used_libs = l
          break
        env_restore(ctx.env, e_bkp)
      if res[0]:
        ctx.Log("this was successful (using "+str(l)+"), done with list of libraries\n")
        break
      env_restore(ctx.env, l_bkp)
    return res

  ## Look if the header files self.headers are present in one of inc_dirs
  # return True or False
  def try_headers(self, ctx, inc_dirs, **kwargs):
    ctx.Log('Trying to find headers in %s\n'%repr(inc_dirs))
    new_inc_dirs = []
    # if there are no headers to be found it is successful
    if not self.headers:
      return True
    
    # if the entries in headers is only a single list, make it a list of listst
    if self.headers:
      if not isinstance(self.headers[0], list):
        self.headers = [self.headers]
    
    # here, self.headers is e.g. [["lapacke.h"], ["mkl_lapacke.h"]]
    found_headers = False
    for headers_option in self.headers:
      found_headers = True
    
      ctx.Log('Try headers ' + str(headers_option) + '\n')
    
      for (i,hdr) in enumerate(headers_option):
        found = False
        for path in inc_dirs:
          hdr_path = os.path.join(path, hdr)
          ctx.Log(' ' + hdr_path + ' ... ')
          if os.path.exists(hdr_path) and not os.path.isfile(hdr_path):
            ctx.Log('(is directory) ')
            
          if os.path.exists(hdr_path) and os.path.isfile(hdr_path):
            ctx.Log('yes.\n')
            found = True
            break
            
          # remove leading "../" and see if file is there
          if hdr_path.find("../") == 0:
            new_hdr_path = hdr_path[3:]
            new_path = path[3:]
            ctx.Log('no.\n')
            ctx.Log(' ' + new_hdr_path + ' ... ')
            if os.path.exists(new_hdr_path):
              #new_inc_dirs.append(new_path)
              ctx.Log('(yes, here it is, but this directory is not considered)\n')
              #found = True
              break
            
          ctx.Log('no.\n')
          
          # look in subdirectories
          for (subpath, subdirectories, files) in os.walk(path):
              
            new_path = os.path.join(path, subpath)
            hdr_path = os.path.join(new_path, hdr)
            ctx.Log(' ' + hdr_path + ' ... ')
            
            if os.path.exists(hdr_path) and not os.path.isfile(hdr_path):
              ctx.Log('(is directory) ')
              
            if os.path.exists(hdr_path) and os.path.isfile(hdr_path):
              ctx.Log('yes.\n')
              new_inc_dirs.append(new_path)
              found = True
              break
            ctx.Log('no.\n')
          
          if found:
            break
        
        if not found:
          ctx.Log('Failed to find ' + hdr + '\n')
          found_headers = False
          break
        
      if found_headers:
        break
        
      
      
    if new_inc_dirs:
      ctx.Log('add more inc_dirs: '+str(inc_dirs))
    
    inc_dirs += new_inc_dirs    # append path with new directories
    
    if new_inc_dirs:
      ctx.Log(' -> '+str(inc_dirs)+'\n')
    return found_headers

  def try_location(self, ctx, base_dirs, **kwargs):
    if type(base_dirs) is not list:
      base_dirs = [base_dirs]
    for base in base_dirs:
      ctx.Log('Checking for %s in %s.\n'%(self.name, base))
      ctx.Log('Searching for headers: '+str(self.headers)+", libs: "+str(self.libs)+'\n')
      
      loc_callback = kwargs.get('loc_callback', None)
      libs = copy.deepcopy(conv.to_iter(self.libs))
      extra_libs = copy.deepcopy(conv.to_iter(self.extra_libs))

      sub_dirs = conv.to_iter(self.sub_dirs)
      if not sub_dirs:
        sub_dirs = [[]]

      ctx.Log("Try the following combinations of (include, lib) directories: "+str(sub_dirs)+"\n")
      res = (False, None)
      for inc_sub_dirs, lib_sub_dirs in sub_dirs:
        inc_sub_dirs = [inc_sub_dirs]
        lib_sub_dirs = [lib_sub_dirs]
        #ctx.Log(inc_sub_dirs)
        #ctx.Log("\n")
        #ctx.Log(lib_sub_dirs)
        #ctx.Log("\n")
        #inc_sub_dirs = list(inc_sub_dirs)
        #lib_sub_dirs = list(lib_sub_dirs)

        for i in range(len(inc_sub_dirs)):
          if not os.path.isabs(inc_sub_dirs[i]):
            inc_sub_dirs[i] = os.path.join(base, inc_sub_dirs[i])
        for i in range(len(lib_sub_dirs)):
          if not os.path.isabs(lib_sub_dirs[i]):
            lib_sub_dirs[i] = os.path.join(base, lib_sub_dirs[i])

        # Remove any directories that can already be found
        # in their respective lists.
        #inc_sub_dirs = [d for d in inc_sub_dirs if d not in ctx.env.get('CPPPATH', [])]
        #lib_sub_dirs = [d for d in lib_sub_dirs if d not in ctx.env.get('LIBPATH', [])]

        if loc_callback:
          loc_callback(ctx, base, inc_sub_dirs, lib_sub_dirs, libs, extra_libs)

        ctx.Log('Trying include directories: "' + str(inc_sub_dirs) + '" and library directories: "' + str(lib_sub_dirs) + '"\n')

        # Before continuing, try and find all of the sample headers.
        if not self.try_headers(ctx, inc_sub_dirs, **kwargs):
          continue

        system_inc_dirs = []
        for inc_dir in inc_sub_dirs:
          if "pgcc" in ctx.env["cc"] or "pgcc" in ctx.env["CC"] or "/usr/include" in inc_dir:     
            system_inc_dirs.append(('-I', inc_dir))
          else:
            system_inc_dirs.append(('-isystem', inc_dir))     # -isystem is the same is -I for gcc, except it suppresses warning (useful for dependencies)
            # also, -isystem /usr/include does not work properly

        if self.set_rpath:
          bkp = env_setup(ctx.env,
                  #CPPPATH=ctx.env.get('CPPPATH', []) + inc_sub_dirs,
                  LIBPATH=ctx.env.get('LIBPATH', []) + lib_sub_dirs,
                  RPATH=ctx.env.get('RPATH', []) + lib_sub_dirs,
                  CCFLAGS=ctx.env.get('CCFLAGS', []) + system_inc_dirs)
        else:
          bkp = env_setup(ctx.env,
                  #CPPPATH=ctx.env.get('CPPPATH', []) + inc_sub_dirs,
                  LIBPATH=ctx.env.get('LIBPATH', []) + lib_sub_dirs,
                  #RPATH=ctx.env.get('RPATH', []) + lib_sub_dirs,
                  CCFLAGS=ctx.env.get('CCFLAGS', []) + system_inc_dirs)
        
        self.base_dir = base  # set base directory (is needed by try_libs)
        res = self.try_libs(ctx, libs, extra_libs, **kwargs)
        if res[0]:
          self.base_dir = base # set base directory
          
          ctx.Log("Combination of (include, lib) directories "+str(inc_sub_dirs) + "," + str(lib_sub_dirs)+" was successful, done with combinations.\n")
          return res
          #break
        env_restore(ctx.env, bkp)
    return res

  def have_option(self, env, name):
    if name in env and env[name] is not False:
      return True

    # Only check the OS environment if no other options for this package were given.
    for opt in self.options:
      if opt in env:
        return False
    return name in env['ENV']

  def have_all_opts(self, env, *names):
    for n in names:
      if not self.have_option(env, n):
        return False
    return True

  def have_any_options(self, env, *names):
    for n in names:
      if self.have_option(env, n):
        return True
    return False

  def get_option(self, env, name):
    if name in env:
      return env[name]
    for opt in self.options:
      if opt in env:
        return ''
    return env['ENV'].get(name, '')

  def check_options(self, env):
    name = self.name.upper()

    # Either base or include/library paths.
    if self.have_option(env, name + '_DIR') and self.have_any_options(env, name + '_INC_DIR', name + '_LIB_DIR'):
      print('\n')
      print(('Please specify either %s_DIR or either of %s_INC_DIR or'%(name, name)))
      print(('%s_LIB_DIR.\n'%name))
      env.Exit(1)

    # Either download or location. (both is allowed, download overrides given directory)
    #elif self.have_option(env, name + '_DOWNLOAD') and self.have_any_options(env, name + '_DIR', name + '_INC_DIR', name + '_LIB_DIR'):
    #  print('\n')
    #  print('Cannot specify to download %s and also give a system location.'%self.name)
    #  env.Exit(1)

  def need_cmake(self, env):
    if not self.have_cmake():
      print('\n')
      print(('%s requires CMake to be installed to autobuild.'%self.name))
      print("")
      env.Exit(1)

  def have_cmake(self):
    if getattr(self, '_cmake', False):
      return self._cmake
    try:
      subprocess.check_call('cmake', stdout=subprocess.PIPE, stderr=subprocess.PIPE)
      self._cmake = True
    except:
      self._cmake = False
    return self._cmake
