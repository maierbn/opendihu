cd .. && time python2.7 $OPENDIHU_HOME/dependencies/scons/scons.py BUILD_TYPE=release -j 4 $* && (cd - && echo release build succeeded) || (echo release build failed && cd - && exit -1)
