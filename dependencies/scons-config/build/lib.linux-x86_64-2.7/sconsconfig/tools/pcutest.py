from SCons.Script import *

class ToolPCUTestWarning(SCons.Warnings.Warning):
    pass

SCons.Warnings.enableWarningClass(ToolPCUTestWarning)

def to_list(var):
    if isinstance(var, str):
        return [var]
    elif var is None:
        return []
    elif isinstance(var, list):
        return var
    elif isinstance(var, tuple):
        return list(var)
    else:
        return [var]

def multiget(dicts, key, default=None):
    for d in dicts:
        if d.has_key(key):
            return d[key]
    else:
        return default

def build_suite_runner(env, target, suites, **kw):
    proto_txt = ''
    suite_txt = ''
    init = multiget([kw, env], 'PCU_INIT', '_init')
    name_suffix = multiget([kw, env], 'PCU_NAME_SUFFIX', '_t')
    mpi_init = multiget([kw, env], 'PCU_MPI_INIT', True)
    mpi_setup = multiget([kw, env], 'PCU_MPI_SETUP', '')
    if mpi_init:
        mpi_init = '1'
    else:
        mpi_init = '0'

    setup = multiget([kw, env], 'PCU_SETUP', '')
    if setup:
        setup = '\n   ' + setup

    teardown = multiget([kw, env], 'PCU_TEARDOWN', '')
    if teardown:
        teardown = '\n   ' + teardown

    for s in suites:
        name = os.path.splitext(os.path.basename(s.path))[0]
        suite_txt += '   pcu_runner_add_suite(%s, %s);\n'%(name + name_suffix, name + init)
        proto_txt += 'void %s(pcu_suite_t *);\n'%(name + init)

    src = '''#include <stdlib.h>
#include <mpi.h>
#include <pcu/pcu.h>

%s
int main(int argc, char* argv[]) {
   int mpi_init = %s;
   int flag;
   pcu_listener_t* lsnr;

   if(mpi_init) {
      MPI_Initialized(&flag);
      if(!flag)
         MPI_Init(&argc, &argv);
   }
   %s
   pcu_runner_init(argc, argv);%s

   %s
   lsnr = pcu_text_output_create();
   pcu_runner_run(lsnr);
   pcu_text_output_destroy(lsnr);
%s
   pcu_runner_finalise();
   if(mpi_init) {
      MPI_Initialized(&flag);
      if(flag)
         MPI_Finalize();
   }
   return EXIT_SUCCESS;
}
'''%(proto_txt, mpi_init, mpi_setup, setup, suite_txt, teardown)

    dir_path = os.path.dirname(target.abspath)
    if not os.path.exists(dir_path):
        os.makedirs(dir_path)
    f = open(target.abspath, 'w')
    f.write(src)
    f.close()
    return File(target.abspath)

def generate(env, **kw):
    env.SetDefault(PCUTEST_TARGET='check')

    def PCUSuite(env, target, source, **kw):
        '''Create an object/header pair out of a
        *Suite.c/*Suite.h pair. The target should just
        be the name of the suite. So, if target were
        'Happy', the sources would be 'HappySuite.c' and
        'HappySuite.h' '''
        obj = env.SharedObject(target[0], source[0])
        return [obj + [File(os.path.splitext(source[0].abspath)[0] + '.h')]]

    def PCUTest(env, target, suites, np=1, source=[], ext='.c', **kw):
        # Generate a list of headers, one for each suite source.
        objs = []
        for s in suites:
            objs.append(s)
        objs.extend(to_list(source))

        # Generate the program source.
        prog_src = build_suite_runner(env, File(str(target[0]) + ext), suites, **kw)

        # Build everything.
        objs = env.SharedObject(os.path.splitext(prog_src.abspath)[0], prog_src) + objs
        libs = multiget([kw, env], 'LIBS', [])
#        if 'pcu' not in libs:
#            kw['LIBS'] = list(libs) + ['pcu']
        test = env.Program(target[0], objs, **kw)

        mpi_cmd = '%s %s %d '%(env.get('MPIEXEC', 'mpiexec'), env.get('MPINP', '-np'), np)
        runner = env.Action('@' + mpi_cmd + test[0].path)
        env.Alias(env['PCUTEST_TARGET'], test, runner)
        env.AlwaysBuild(env['PCUTEST_TARGET'])
        return test

    env.Append(BUILDERS={'PCUSuite': PCUSuite, 'PCUTest': PCUTest})

def exists(env):
    # Should probably have this search for the pcu
    # libraries/source or something.
    return True
