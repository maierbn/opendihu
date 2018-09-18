
import sys, os
from Package import Package
import subprocess

class MPI(Package):

  def __init__(self, **kwargs):
    defaults = {
      'download_url': 'http://www.mcs.anl.gov/research/projects/mpich2/downloads/tarballs/1.4.1p1/mpich2-1.4.1p1.tar.gz',
    }
    defaults.update(kwargs)
    super(MPI, self).__init__(**defaults)
    
    if False:
      defaults = {
        'download_url': 'http://www.mcs.anl.gov/research/projects/mpich2/downloads/tarballs/1.4.1p1/mpich2-1.4.1p1.tar.gz',
      }
      defaults.update(kwargs)
      super(MPI, self).__init__(**defaults)
    
      self.mpi_compilers = ['mpicc', 'mpic++', 'mpicxx',
                            'mpif77', 'mpif90', 'mpif']
      self.headers = ['mpi.h']
      self.libs=[
        ['mpich'],                   # mpich/mpich2
        ['pmpich', 'mpich'],         # mpich2
        ['mpich', 'mpl'],            # mpich2
        ['pmpich', 'mpich', 'mpl'],  # mpich2
        ['mpi', 'mpi_cxx'],          # openmpi
        ['mpi'],                     # openmpi
      ]
      self.extra_libs=[
        [],
        ['rt'],
        ['pthread', 'rt'],
        ['dl'],
        ['dl', 'rt'],
        ['dl', 'pthread'],
        ['dl', 'pthread', 'rt']
      ]
    self.check_text = r'''
#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
int main(int argc, char* argv[])
{
  int rank;
  MPI_Init(&argc, &argv);
  printf("%d\n", MPI_VERSION);
  printf("%d\n", MPI_SUBVERSION);

  /* This will catch OpenMPI libraries. */
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  MPI_Finalize();
  return EXIT_SUCCESS;
}
'''

    # MPI on a cluster usually won't run properly, so don't
    # try to.
    # self.run = False

    # Setup the build handler. I'm going to assume this will work for all architectures.
    self.set_build_handler([
        './configure --prefix=${PREFIX} --enable-shared --disable-fc --disable-f77',
        'make',
        'make install'
    ])

  def check(self, ctx):
    env = ctx.env
    ctx.Message('Checking for MPI ... ')
    self.check_options(env)

    # remove specified flags
    try:
      #cflags = subprocess.check_output("python-config --cflags", shell=True)+' '
      cflags = subprocess.check_output("mpicc --showme:compile", shell=True)
      ldflags = subprocess.check_output("mpicc --showme:link", shell=True)

      #for flag_to_remove in flags_to_remove:
      #  while flag_to_remove in cflags:
      #    startpos = cflags.index(flag_to_remove)
      #    length = len(flag_to_remove)
      #    while cflags[startpos+length] == ' ':
      #      length += 1
      #    cflags = cflags[0:startpos] + cflags[startpos+length:]

      # remove trailing newline
      if cflags[-1] == '\n':
        cflags = cflags[:-1]
      if ldflags[-1] == '\n':
        ldflags = ldflags[:-1]

      ctx.Log("extracted cflags from mpicc: \n{}\n\n".format(cflags))
      ctx.Log("extracted ldflags from mpicc: \n{}\n\n".format(ldflags))

      env.MergeFlags(cflags)
      env.MergeFlags(ldflags)
      
      res = self.try_link(ctx)
      
    except Exception as e: 
      ctx.Message("MPI mpicc --showme failed: "+str(e)+", now using MPI_DIR")

      res = super(MPI, self).check(ctx)
    
    
    self.check_required(res[0], ctx)
    ctx.Result(res[0])
    return res[0]
