cd /store/software/opendihu && scons BUILD_TYPE=debug no_tests=True && (cd - && echo opendihu build without tests succeeded) || (echo opendihu build without tests failed && cd - && exit -1)
