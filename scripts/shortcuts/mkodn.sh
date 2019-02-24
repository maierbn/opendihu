cd $OPENDIHU_HOME && scons BUILD_TYPE=debug no_tests=True && (cd - && echo opendihu debug build without tests succeeded) || (echo opendihu debug build without tests failed && cd - && exit -1)
