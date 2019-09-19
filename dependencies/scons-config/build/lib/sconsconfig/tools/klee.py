from SCons.Script import *

def generate(env):
    pass

def exists(env):
    return env.Detect('klee')
