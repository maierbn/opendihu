cd $OPENDIHU_HOME && time scons BUILD_TYPE=debug no_tests=True -j 4 && (cd - && echo opendihu debug build without tests succeeded) || (echo opendihu debug build without tests failed && cd - && exit -1)
