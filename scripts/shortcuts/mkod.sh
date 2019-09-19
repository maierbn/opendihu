cd $OPENDIHU_HOME && time scons BUILD_TYPE=debug -j 4 && (cd - && echo opendihu debug succeeded) || (echo opendihu debug failed && cd - && exit -1)
