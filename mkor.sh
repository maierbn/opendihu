cd /store/software/opendihu && scons BUILD_TYPE=release && (cd - && echo opendihu build succeeded) || (echo opendihu build failed && cd - && exit -1)
