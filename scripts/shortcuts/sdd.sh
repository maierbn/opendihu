cd .. && scons BUILD_TYPE=debug -j 4 && (cd - && echo debug build succeeded) || (echo debug build failed && cd - && exit -1)
