cd $OPENDIHU_HOME && time $OPENDIHU_HOME/dependencies/scons/scons.py BUILD_TYPE=debug -j 4 $* && (cd - && echo opendihu debug succeeded) || (echo opendihu debug failed && cd - && exit -1)
