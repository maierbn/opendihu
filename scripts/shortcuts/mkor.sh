cd $OPENDIHU_HOME && scons BUILD_TYPE=release -j 4 && (cd - && echo opendihu release build succeeded) || (echo opendihu release build failed && cd - && exit -1)
