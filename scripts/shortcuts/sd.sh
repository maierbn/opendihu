time python2.7 $OPENDIHU_HOME/dependencies/scons/scons.py BUILD_TYPE=debug -j 4 $* && (echo debug build succeeded) || (echo debug build failed && exit -1)
