.PHONY: clean all

all: debug release

ifneq ("$(wildcard ./dependencies/python/install/bin/python)","")
python := ./dependencies/python/install/bin/python
else
python := python
endif

debug:
	$(python) dependencies/scons/scons.py -Q BUILD_TYPE=DEBUG

release:
	$(python) dependencies/scons/scons.py -Q BUILD_TYPE=RELEASE

clean:
	rm -rf .sconf_temp
	rm .sconsign.dblite	
	$(python) dependencies/scons/scons.py -Q -c BUILD_TYPE=DEBUG
	$(python) dependencies/scons/scons.py -Q -c BUILD_TYPE=RELEASE

functional_testing:
	cd testing/functional_testing && $(MAKE) default

solid_mechanics:
	cd testing/functional_testing/tests/solid_mechanics && python ../../../../dependencies/scons/scons.py -Q BUILD_TYPE=DEBUG
