cd $OPENDIHU_HOME && scons BUILD_TYPE=debug && (cd - && echo opendihu debug succeeded) || (echo opendihu debug failed && cd - && exit -1)
