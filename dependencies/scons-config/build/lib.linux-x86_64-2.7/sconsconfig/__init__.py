from utils import conv
import packages

enabled_packages = []
package_map = {}
custom_tests = {}

def select(*args):
    packages = args
    for pkg in packages:
        if pkg in enabled_packages:
            continue
        enabled_packages.append(pkg)
        package_map[pkg.__module__] = pkg
        custom_tests.update(pkg.custom_tests)

def package(klass):
    return package_map.get(klass.__module__, None)

def add_options(vars):
    for pkg in enabled_packages:
        pkg.add_options(vars)

def check(sconf):
    for pkg in enabled_packages:
        for name in conv.to_iter(pkg.test_names):
            getattr(sconf, name)()

def configure(env, vars):

    # Before I save the issued options, I need to strip out one-shot options
    # so they don't get repeated later.
    bkp = {}
    for pkg in enabled_packages:
        for opt_name in pkg.one_shot_options:
            bkp[opt_name] = env[opt_name]
            del env[opt_name]

    # Save the configuration options here, but not again. This is because certain
    # options will set other options in the background, and we don't want to mistakenly
    # think they have been set by the user.
    vars.Save('config.py', env)

    # Restore any one-shot options we stripped out earlier.
    for k,v in bkp.iteritems():
        env[k] = v

    # Perform the configuration.
    sconf = env.Configure(custom_tests=custom_tests)
    check(sconf)
    sconf.Finish()
