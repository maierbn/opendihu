cd .. && time scons BUILD_TYPE=release -j 4 && (cd - && echo release build succeeded) || (echo release build failed && cd - && exit -1)
