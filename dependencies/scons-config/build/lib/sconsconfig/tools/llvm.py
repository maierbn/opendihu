from SCons.Script import *

def exists(env):
    return env.Detect('llvm-gcc') and env.Detect('llvm-ld')

def generate(env):
    env.SetDefault(LLVMCC='llvm-gcc')
    env.SetDefault(LLVMLINK='llvm-ld')
    if not exists(env):
        print 'Error: Could not find either or both of %s and %s.'%(repr(env['LLVMCC']), repr(env['LLVMLINK']))
        env.Exit(1)
        return

    env['BUILDERS']['LlvmObject'] = SCons.Builder.Builder(
        action=SCons.Action.Action("$LLVMCCCOM", "$LLVMCCCOMSTR"),
        emitter=SCons.Defaults.StaticObjectEmitter,
        prefix='$LLVMOBJPREFIX',
        suffix='$LLVMOBJSUFFIX',
        src_builder=['CFile', 'CXXFile'],
        source_scanner=SourceFileScanner,
        single_source=1)

    env.SetDefault(LLVMOBJSUFFIX='.bc')
    env['LLVMCCCOM'] = '$LLVMCC -o $TARGET -c $CFLAGS $CCFLAGS $_CCCOMCOM $SOURCES'

    env['BUILDERS']['LlvmProgram'] = SCons.Builder.Builder(
        action=SCons.Action.Action("$LLVMLINKCOM", "$LLVMLINKCOMSTR"),
        emitter='$PROGEMITTER',
        prefix='$PROGPREFIX',
        suffix='$LLVMPROGSUFFIX',
        src_suffix='$LLVMOBJSUFFIX',
        src_builder='LlvmObject',
        target_scanner=ProgramScanner)

    env['LLVMLINKCOM'] = '$LLVMLINK -o $TARGET $LLVMLINKFLAGS $SOURCES $_LIBDIRFLAGS $_LIBFLAGS'
    env['LLVMLINKFLAGS'] = []
    env['LLVMPROGSUFFIX'] = '.llvm'
