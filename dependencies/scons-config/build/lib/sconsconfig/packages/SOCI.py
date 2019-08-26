import sys, os
from Package import Package
import sconsconfig as config

class SOCI(Package):

    def __init__(self, **kwargs):
        self.backends = set(kwargs.pop('backends', []))
        defaults = {
            'download_url': 'http://downloads.sourceforge.net/project/soci/soci/soci-3.1.0/soci-3.1.0.zip',
        }
        defaults.update(kwargs)
        super(SOCI, self).__init__(**defaults)
        self.ext = '.cc'
        self.sub_dirs = [
            (('include', 'include/soci'), 'lib'),
            (('include', 'include/soci'), 'lib64'),
        ]
        self.libs=[
            ['soci_core']
        ]
        self.extra_libs=[
        ]
        self.check_text = r'''
#include <stdlib.h>
#include <stdio.h>
#include <soci/soci.h>
int main(int argc, char* argv[]) {
   return EXIT_SUCCESS;
}
'''

    def check(self, ctx):
        env = ctx.env
        ctx.Message('Checking for SOCI ... ')
        self.check_options(env)

        # SOCI can use Boost, so check to see if Boost is in our
        # set of configuration options and set accordingly.
        cmake = 'cmake -DCMAKE_INSTALL_PREFIX:PATH=${PREFIX}'
        boost = config.package(config.packages.boost)
        # if boost and boost.found and boost.base_dir:
        #     cmake += ' -DBOOST_DIR:PATH=' + boost.base_dir
        cmake += ' -DWITH_BOOST=off'

        # Turn on release mode to prevent a very odd bug
        # in SOCI that produces an assertion whenever you
        # try to read data from PostgreSQL into a vector
        # larger than numeric_limits<unsigned short>::max().
        cmake += ' -DCMAKE_BUILD_TYPE:STRING=release'

        # Check for sqlite3, like boost.
        sqlite = config.package(config.packages.sqlite3)
        if sqlite and sqlite.found and sqlite.base_dir:
            cmake += ' -DSQLITE3_INCLUDE_DIR:PATH=' + sqlite.include_directories()
            cmake += ' -DSQLITE3_LIBRARY:FILEPATH=' + sqlite.libraries()

        # Check for MySQL.
        pkg = config.package(config.packages.MySQL)
        if pkg and pkg.found and pkg.base_dir:
            cmake += ' -DMYSQL_INCLUDE_DIR:PATH=' + pkg.include_directories()
            cmake += ' -DMYSQL_LIBRARY:FILEPATH=' + pkg.libraries()

        # Check for PostgreSQL.
        pkg = config.package(config.packages.PostgreSQL)
        if pkg and pkg.found and pkg.base_dir:
            cmake += ' -DPOSTGRESQL_INCLUDE_DIR:PATH=' + pkg.include_directories()
            cmake += ' -DPOSTGRESQL_LIBRARY:FILEPATH=' + pkg.libraries()

        # For some reason SOCI is incompatible with gcc 4.7.1.
        # Need to switch off testing and the empty thingy.
        cmake += ' -DSOCI_TEST:BOOL=off -DSOCI_EMPTY:BOOL=off'

        cmake += ' .'
        self.set_build_handler([
            cmake,
            'make',
            'make install'
        ])

        # Check each backend for existence.
        backends = {
            'sqlite3': {
                'libs': ['soci_core', 'soci_sqlite3'],
            },
            'mysql': {
                'libs': ['soci_core', 'soci_mysql'],
                'extra_libs': ['dl', 'mysqlclient'],
            },
            'postgresql': {
                'libs': ['soci_core', 'soci_postgresql'],
                'extra_libs': ['dl', 'pq'],
            },
        }
        found = False
        for be, opts in backends.iteritems():
            self.libs = [opts['libs']]
            self.extra_libs = [opts.get('extra_libs', [])]
            res = super(SOCI, self).check(ctx)
            if res[0]:
                found = True
                env.MergeFlags('-DHAVESOCI' + be.upper())

            # If this is a required backend, terminate with failure.
            elif be in self.backends:
                found = False
                break

        self.check_required(found, ctx)
        ctx.Result(found)
        return found
