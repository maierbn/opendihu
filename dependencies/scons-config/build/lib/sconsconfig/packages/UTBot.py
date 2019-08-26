import sys, os
import Package
from Package import have_any_opts, try_link

utbot_text = r'''
#include <stdlib.h>
#include <stdio.h>
#include <utbot/utbot.hh>
int main(int argc, char* argv[]) {
   return EXIT_SUCCESS;
}
'''

utbot_libs=[['utbot']]

def CheckUTBot(ctx, required=True):
    env = ctx.env
    ctx.Message('Checking for UTBot ... ')
    Package.check_options(env, 'UTBot')

    res = Package.CheckPkg(ctx, 'UTBot', utbot_text, utbot_libs, ext='.cc')

    Package.Required('UTBot', res[0], required)
    ctx.Result(res[0])
    return res[0]

def AddOptions(vars):
    Package.AddOptions(vars, 'UTBot')
