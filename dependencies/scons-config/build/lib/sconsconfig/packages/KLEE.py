import sys, os
import Package
from Package import have_any_opts, try_link

klee_text = r'''
#include <stdlib.h>
#include <stdio.h>
#include <klee/klee.h>
int main(int argc, char* argv[]) {
   int var;
   klee_make_symbolic(&var, sizeof(int), "var");
   return EXIT_SUCCESS;
}
'''

klee_libs=[['kleeRuntest']]

def CheckKLEE(ctx, required=True):
    env = ctx.env
    ctx.Message('Checking for KLEE ... ')
    Package.check_options(env, 'KLEE')

    res = Package.CheckPkg(ctx, 'KLEE', klee_text, klee_libs, auto_add_libs=False, run=False)

    Package.Required('KLEE', res[0], required)
    ctx.Result(res[0])
    return res[0]

def AddOptions(vars):
    Package.AddOptions(vars, 'KLEE')
