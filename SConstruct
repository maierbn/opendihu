# This is the top-level SConstruct
# It calls the SConscript files for building the core and the examples

# build core
SConscript('core/SConstruct')

# build tests
if 'no_tests' not in ARGUMENTS or 'travis_ci' in ARGUMENTS:
  SConscript('testing/unit_testing/SConstruct')

# build examples
if 'no_examples' not in ARGUMENTS:
  SConscript('examples/SConstruct')
