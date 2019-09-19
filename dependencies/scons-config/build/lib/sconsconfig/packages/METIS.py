import sys, os
import Package
from Package import have_any_opts, try_link

metis_text = r'''
#include <stdlib.h>
#include <stdio.h>
#include <metis.h>
int main(int argc, char* argv[]) {
   return EXIT_SUCCESS;
}
'''

metis_libs=[['metis']]

def CheckMETIS(ctx, required=True):
    env = ctx.env
    ctx.Message('Checking for METIS ... ')
    Package.check_options(env, 'METIS')

    res = Package.CheckPkg(ctx, 'METIS', metis_text, metis_libs)

    Package.Required('METIS', res[0], required)
    ctx.Result(res[0])
    return res[0]

def AddOptions(vars):
    Package.AddOptions(vars, 'METIS')
