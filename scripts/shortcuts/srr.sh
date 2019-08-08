cd .. && scons BUILD_TYPE=release && (cd - && echo release build succeeded) || (echo release build failed && cd - && exit -1)
