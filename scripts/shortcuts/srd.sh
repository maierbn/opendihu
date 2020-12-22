time $OPENDIHU_HOME/dependencies/scons/scons.py BUILD_TYPE=releasewithdebuginfo -j 4 $* && (echo release with debug info build succeeded) || (echo release with debug info build failed && exit -1)
