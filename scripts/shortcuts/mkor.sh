cd $OPENDIHU_HOME && scons BUILD_TYPE=release && (cd - && echo opendihu release build succeeded) || (echo opendihu release build failed && cd - && exit -1)
