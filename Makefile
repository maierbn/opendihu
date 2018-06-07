.PHONY: clean all

all: debug release

#ifneq ("$(wildcard ./dependencies/python/install/bin/python)","")
#python := ./dependencies/python/install/bin/python
#else
python := python2.7
#endif

debug:
	$(python) dependencies/scons/scons.py BUILD_TYPE=DEBUG

release:
	$(python) dependencies/scons/scons.py BUILD_TYPE=RELEASE

clean:
	rm -rf .sconf_temp
	rm .sconsign.dblite	

purge: clean
	rm -rf core/build_debug
	rm -rf core/build_release

purge_dependencies:
	cd dependencies; rm -rf base64/ bzip2/ cython/ easyloggingpp/ googletest/ lapack/ matplotlib/ numpyc/ petsc/ python/ scipy/ semt/; cd -

rebuild: purge_dependencies purge clean debug release

system_testing:
	cd testing/system_testing && ./run.sh

solid_mechanics:
	cd testing/system_testing/tests/solid_mechanics && python ../../../../dependencies/scons/scons.py BUILD_TYPE=DEBUG

multiple_fibers:
	cd examples/multiple_fibers && python ../../dependencies/scons/scons.py BUILD_TYPE=DEBUG

streamline_tracer:
	cd testing/system_testing/tests/fibers && python ../../../../dependencies/scons/scons.py BUILD_TYPE=DEBUG

diffusion:
	cd testing/system_testing/tests/diffusion &&  python ../../../../dependencies/scons/scons.py BUILD_TYPE=DEBUG

laplace:
	cd testing/system_testing/tests/laplace &&  python ../../../../dependencies/scons/scons.py BUILD_TYPE=DEBUG

quadrature:
	cd examples/quadrature/own && python ../../../dependencies/scons/scons.py BUILD_TYPE=DEBUG

fibers:
	cd testing/system_testing/tests/fibres &&  python ../../../../dependencies/scons/scons.py BUILD_TYPE=DEBUG
