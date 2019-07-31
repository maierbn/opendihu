import os, subprocess
from SCons.Script import *

utbot_tests = []

def multiget(dicts, key, default=None):
    for d in dicts:
        if d.has_key(key):
            return d[key]
    else:
        return default

def UTBotHarness(env, target, source, **kw):
    source = File(source[0])
    name = os.path.splitext(os.path.basename(str(source)))[0]
    path = os.path.dirname(source.srcnode().path)
    # Build the binary.
    bin = env.Program(name, [source] + kw.get('OBJS', []), LIBS=multiget([kw, env], 'LIBS', []) + ['utbot'])
    # Build the LLVM bit-code.
    if multiget([kw, env], 'BUILD_LLVM', None):
        o = env.LlvmObject(source,
                           CFLAGS=multiget([kw, env], 'CFLAGS', []) + ['--emit-llvm'],
                           CPPDEFINES=multiget([kw, env], 'CPPDEFINES', []) + ['HAVE_KLEE'])
        bc = env.LlvmProgram(name, o + kw.get('LLVMOBJS', []),
                             LLVMLINKFLAGS=multiget([kw, env], 'LLVMLINKFLAGS', []) + ['--disable-opt'],
                             LIBS=multiget([kw, env], 'LLVMLIBS', []) + ['utbot_llvm'])
    else:
        bc = 'none'
    # Convert the binary and bit-code locations.
    bin = bin[0].path
    if bc is not 'none':
        bc = bc[0].path + '.bc'
    # Add the test and information to the list.
    utbot_tests.append((name, path, bin, bc))

def sync(target, source, env):
    for test in utbot_tests:
        name, path, bin, bc = test
        cmd = 'utbot query -t %s'%name
        subp = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = subp.communicate()
        if subp.returncode != 0:
            # Initialise the test.
            cmd = 'utbot init -t %s %s %s %s'%(name, path, bin, bc)
        else:
            # Test exists, update it.
            cmd = 'utbot update -t %s'%name
        subp = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = subp.communicate()

def purge(target, source, env):
    cmd = 'utbot query -n'
    subp = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = subp.communicate()
    existing = stderr.splitlines()
    names = [test[0] for test in utbot_tests]
    for name in existing:
        if name not in names:
            cmd = 'utbot remove -t %s'%name
            subp = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            stdout, stderr = subp.communicate()

def exists(env):
    return env.Detect(env['UTBOT'])

def generate(env):
    env.SetDefault(UTBOT='utbot')
    env.SetDefault(UTBOT_CHECK_TARGET='utbot-check')
    env.SetDefault(UTBOT_SYNC_TARGET='utbot-sync')
    env.SetDefault(UTBOT_PURGE_TARGET='utbot-purge')
    if not exists(env):
        print 'Error: Unable to locate UTBot command %s.'%repr(env['UTBOT'])
        env.Exit(1)
    env.Append(BUILDERS={'UTBotHarness': UTBotHarness})
    env.Alias(env['UTBOT_CHECK_TARGET'], None, env.Action('@utbot run'))
    env.AlwaysBuild(env['UTBOT_CHECK_TARGET'])
    env.Alias(env['UTBOT_SYNC_TARGET'], None, env.Action(sync))
    env.AlwaysBuild(env['UTBOT_SYNC_TARGET'])
    env.Alias(env['UTBOT_PURGE_TARGET'], None, env.Action(purge))
    env.AlwaysBuild(env['UTBOT_PURGE_TARGET'])
