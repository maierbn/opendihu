cd .. && time python2.7 $OPENDIHU_HOME/dependencies/scons/scons.py BUILD_TYPE=debug -j 4 $* && (cd - && echo debug build succeeded) || (echo debug build failed && cd - && exit -1)
