.PHONY: clean all doc

all: debug_without_tests release

#ifneq ("$(wildcard ./dependencies/$(python)/install/bin/python3)","")
#	python=./dependencies/$(python)/install/bin/python3
#else
python=python3
#endif

debug:
	$(python) dependencies/scons/scons.py BUILD_TYPE=DEBUG -j $(shell expr `nproc --all` / 2)

release:
	$(python) dependencies/scons/scons.py BUILD_TYPE=RELEASE -j $(shell expr `nproc --all` / 2)

clean:
	rm -rf .sconf_temp
	rm .sconsign.dblite	

purge: clean
	rm -rf core/build_debug
	rm -rf core/build_release

purge_dependencies:
	cd dependencies; rm -rf adios/ base64/ easyloggingpp/ googletest/ libxml2/ opencor/ petsc/ precice/ python/ pythonpackages/ semt/ std_simd/ vc/ xbraid/; cd -

rebuild: purge_dependencies purge clean debug release

# on hazel hen rebuild everying including dependencies
rebuild_hazelhen:
	unset OPENDIHU_HOME; rm -rf dependencies/easyloggingpp/install dependencies/easyloggingpp/src dependencies/python/install dependencies/python/src  dependencies/base64/install dependencies/base64/src dependencies/numpyc/install dependencies/numpyc/src dependencies/semt/src dependencies/semt/install && rm -rf core/build_release && $(python) dependencies/scons/scons.py BUILD_TYPE=RELEASE; cd dependencies/matplotlib && ../python/install/bin/pip3 install *.whl

doc:
	cd doc/doxygen; doxygen
	
# the following targets are just for convenience and could also be deleted
release_without_tests:
	$(python) dependencies/scons/scons.py BUILD_TYPE=RELEASE no_tests=True -j $(shell expr `nproc --all` / 2)

debug_without_tests:
	$(python) dependencies/scons/scons.py BUILD_TYPE=DEBUG no_tests=True -j $(shell expr `nproc --all` / 2)

travis_ci:
	$(python) dependencies/scons/scons.py BUILD_TYPE=RELEASE no_tests=True travis_ci=True

system_testing:
	cd testing/system_testing && ./run.sh

solid_mechanics:
	cd testing/system_testing/tests/solid_mechanics && $(python) ../../../../dependencies/scons/scons.py BUILD_TYPE=DEBUG

multiple_fibers_system_testing:
	cd testing/system_testing/tests/multiple_fibers && $(python) ../../../../dependencies/scons/scons.py BUILD_TYPE=DEBUG

streamline_tracer:
	cd testing/system_testing/tests/fibers && $(python) ../../../../dependencies/scons/scons.py BUILD_TYPE=DEBUG

diffusion:
	cd testing/system_testing/tests/diffusion &&  $(python) ../../../../dependencies/scons/scons.py BUILD_TYPE=DEBUG

fibers:
	cd testing/system_testing/tests/fibers &&  $(python) ../../../../dependencies/scons/scons.py BUILD_TYPE=DEBUG

hodgkin_huxley:
	cd examples/electrophysiology/monodomain/hodgkin_huxley && $(python) ../../../../dependencies/scons/scons.py BUILD_TYPE=DEBUG

shorten:
	cd examples/electrophysiology/monodomain/shorten && $(python) ../../../../dependencies/scons/scons.py BUILD_TYPE=DEBUG

cellml:
	cd examples/electrophysiology/cellml && $(python) ../../../dependencies/scons/scons.py BUILD_TYPE=DEBUG

multidomain:
	cd examples/electrophysiology/multidomain3d && $(python) ../../../dependencies/scons/scons.py BUILD_TYPE=DEBUG

parallel_fiber_estimation:
	cd examples/parallel_fiber_estimation && $(python) ../../dependencies/scons/scons.py BUILD_TYPE=DEBUG

load_balancing:
	cd examples/load_balancing && $(python) ../../dependencies/scons/scons.py BUILD_TYPE=DEBUG

multiple_fibers:
	cd examples/electrophysiology/multiple_fibers && $(python) ../../../dependencies/scons/scons.py BUILD_TYPE=DEBUG
	
multiple_fibers_cubes_partitioning:
	cd examples/electrophysiology/multiple_fibers_cubes_partitioning && $(python) ../../../dependencies/scons/scons.py BUILD_TYPE=DEBUG
	
fibers_emg:
	cd examples/electrophysiology/fibers/fibers_emg && $(python) ../../../../dependencies/scons/scons.py BUILD_TYPE=DEBUG

laplace2d:
	cd examples/laplace/laplace2d && $(python) ../../../dependencies/scons/scons.py BUILD_TYPE=DEBUG

laplace3d:
	cd examples/laplace/laplace3d_triangular_prism && $(python) ../../../dependencies/scons/scons.py BUILD_TYPE=DEBUG

laplace_surface:
	cd examples/laplace/laplace3d_surface && $(python) ../../../dependencies/scons/scons.py BUILD_TYPE=DEBUG

linear_elasticity:
	cd examples/solid_mechanics/linear_elasticity && $(python) ../../../dependencies/scons/scons.py BUILD_TYPE=DEBUG

fibers_linear_elasticity:
	cd examples/electrophysiology/fibers/fibers_emg && $(python) ../../../../dependencies/scons/scons.py BUILD_TYPE=DEBUG

mooney_rivlin_transiso:
	cd examples/solid_mechanics/mooney_rivlin_transiso  && $(python) ../../../dependencies/scons/scons.py BUILD_TYPE=DEBUG

mooney_rivlin:
	cd examples/solid_mechanics/mooney_rivlin_isotropic  && $(python) ../../../dependencies/scons/scons.py BUILD_TYPE=DEBUG

dynamic_mooney_rivlin:
	cd examples/solid_mechanics/dynamic_mooney_rivlin  && $(python) ../../../dependencies/scons/scons.py BUILD_TYPE=DEBUG

static_biceps_emg:
	cd examples/electrophysiology/static_biceps_emg  && $(python) ../../../dependencies/scons/scons.py BUILD_TYPE=DEBUG

precice0:	
	cd examples/electrophysiology/biceps_contraction/opendihu_precice_opendihu && $(python) ../../../../dependencies/scons/scons.py BUILD_TYPE=DEBUG

biceps_contraction:	
	cd examples/electrophysiology/biceps_contraction/opendihu_opendihu && $(python) ../../../../dependencies/scons/scons.py BUILD_TYPE=DEBUG

multidomain_contraction:	
	cd examples/electrophysiology/multidomain_contraction && $(python) ../../../dependencies/scons/scons.py BUILD_TYPE=DEBUG

laplace2d_split:
	cd examples/laplace/laplace2d_split && $(python) ../../../dependencies/scons/scons.py BUILD_TYPE=DEBUG

linear_elasticity_with_activation:
	cd examples/solid_mechanics/linear_elasticity/with_fiber_activation && $(python) ../../../../dependencies/scons/scons.py BUILD_TYPE=DEBUG

tensile_test:
	cd examples/solid_mechanics/tensile_test && $(python) ../../../dependencies/scons/scons.py BUILD_TYPE=DEBUG

precice1:	
	cd examples/electrophysiology/fibers/fibers_contraction/with_tendons_precice && ../../../../../dependencies/scons/scons.py BUILD_TYPE=DEBUG

multidomain_neuromuscular:
	cd examples/electrophysiology/multidomain_neuromuscular && $(python) ../../../dependencies/scons/scons.py BUILD_TYPE=DEBUG

only_neurons:
	cd examples/electrophysiology/neuromuscular/only_neurons_flat && $(python) ../../../../dependencies/scons/scons.py BUILD_TYPE=DEBUG

contraction_no_precice:
	cd examples/electrophysiology/fibers/fibers_contraction/no_precice && $(python) ../../../../../dependencies/scons/scons.py BUILD_TYPE=DEBUG

with_precice_volume_coupling:
	cd examples/electrophysiology/fibers/fibers_contraction/with_precice_volume_coupling && $(python) ../../../../../dependencies/scons/scons.py BUILD_TYPE=DEBUG
