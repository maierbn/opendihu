cd .. && scons BUILD_TYPE=debug && (cd - && echo debug build succeeded) || (echo debug build failed && cd - && exit -1)
