import subprocess, os, re, tempfile
from SCons.Script import *

def multiget(dicts, key, default=None):
    for d in dicts:
        if d.has_key(key):
            return d[key]
    else:
        return default

def exists(env):
    return env.Detect(env['MEXCC']) and env.Detect('matlab')

def mex_action(env, target, source, **kw):
    cppdefines = env['MEXCPPDEFINES'] + multiget([kw, env], 'CPPDEFINES', [])
    cpppath = env['MEXCPPPATH'] + multiget([kw, env], 'CPPPATH', [])
    linkflags = env['MEXLINKFLAGS'] + multiget([kw, env], 'LINKFLAGS', [])
    libpath = env['MEXLIBPATH'] + multiget([kw, env], 'LIBPATH', [])
    libs = env['MEXLIBS'] + multiget([kw, env], 'LIBS', [])
    env.SharedLibrary(target, source, CPPPATH=cpppath, CPPDEFINES=cppdefines,
                      LINKFLAGS=linkflags, LIBPATH=libpath, LIBS=libs,
                      SHLIBPREFIX='', SHLIBSUFFIX=env['MEXSUFFIX'])

def generate(env):
    env.SetDefault(MEXCC='mex')
    if not exists(env):
        print 'Error: Could not find Matlab\'s %s.'%repr(env['MEXCC'])
        env.Exit(1)
        return

    # Get the plugin extension.
    cmd = 'matlab -nodisplay -nojvm'
    subp = subprocess.Popen(cmd, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    os.write(subp.stdin.fileno(), "mexext\n")
    os.write(subp.stdin.fileno(), "quit\n")
    stdout, stderr = subp.communicate()
    mex_ext = re.search('>>.*\n(.*)\n', stdout).group(1)
    # Set the default.
    env.SetDefault(MEXSUFFIX='.' + mex_ext)

    # Get the mex options.
    tmp_cc = tempfile.NamedTemporaryFile(mode='w', suffix='.cc')
    name = tmp_cc.name
    cmd = 'mex -n ' + name
    subp = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = subp.communicate()
    tmp_cc.close()
    try:
        os.remove(name)
    except:
        pass
    # Compiling.
    match = re.search(r'->.*?(?=-)(.*)\n', stdout)
    mex_opts = env.ParseFlags(match.group(1))
    stdout = stdout[match.end(1):]
    env.SetDefault(MEXCPPDEFINES=[d for d in mex_opts['CPPDEFINES'] if d != 'NDEBUG'])
    env.SetDefault(MEXCPPPATH=mex_opts['CPPPATH'])
    # Linking.
    match = re.search(r'->.*?(?=-)(.*)\n', stdout)
    mex_opts = env.ParseFlags(match.group(1))
    env.SetDefault(MEXLINKFLAGS=mex_opts['LINKFLAGS'])
    env.SetDefault(MEXLIBPATH=mex_opts['LIBPATH'])
    env.SetDefault(MEXLIBS=[l for l in mex_opts['LIBS'] if isinstance(l, str)])

    env['BUILDERS']['MatlabPlugin'] = mex_action
