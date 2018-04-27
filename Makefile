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

system_testing:
	cd testing/system_testing && $(MAKE) default

solid_mechanics:
	cd testing/system_testing/tests/solid_mechanics && python ../../../../dependencies/scons/scons.py BUILD_TYPE=DEBUG

multiple_fibers:
	cd examples/multiple_fibers && python ../../dependencies/scons/scons.py BUILD_TYPE=DEBUG
