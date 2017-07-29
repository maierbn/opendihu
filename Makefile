.PHONY: clean all

all: debug release

debug:
	python dependencies/scons/scons.py -Q BUILD_TYPE=DEBUG

release:
	python dependencies/scons/scons.py -Q BUILD_TYPE=RELEASE

clean:
	python dependencies/scons/scons.py -Q -c BUILD_TYPE=DEBUG
	python dependencies/scons/scons.py -Q -c BUILD_TYPE=RELEASE
