import sys, os, multiprocessing
from Package import Package
import subprocess

##
##  Handles Chaste library
##
class Chaste(Package):

  def __init__(self, **kwargs):
    defaults = {
        'download_url': 'https://github.com/Chaste/Chaste/archive/2018.1.tar.gz',
    }
    defaults.update(kwargs)
    super(Chaste, self).__init__(**defaults)
    self.ext = '.cpp'
    self.sub_dirs = [('include', 'lib/chaste')]
    
    self.headers = [
      'ChasteBuildRoot.hpp',      # Chaste/global/src
      'ArchiveLocationInfo.hpp',  # Chaste/global/src/checkpointing
      'AbstractFileComparison.hpp',           # Chaste/global/src/fortests
      'PetscTools.hpp',                       # Chaste/global/src/parallel
      'Boost165ExponentialDistribution.hpp',  # Chaste/global/src/random
      'GenericEventHandler.hpp',         # Chaste/global/src/timing
      'AbstractHdf5Access.hpp',  # Chaste/io/src/common;
      'AbstractFileComparison.hpp',  # Chaste/io/src/fortests;
      'AbstractDataReader.hpp',  # Chaste/io/src/reader;
      'AbstractDataWriter.hpp',  # Chaste/io/src/writer;
      'AbstractNonlinearSolver.hpp',  # Chaste/linalg/src;
      'predicates.hpp',  # Chaste/mesh/src/3rdparty/tetgen1.4.2;
      'triangle.h',  # Chaste/mesh/src/3rdparty/triangle;
      'AbstractChasteRegion.hpp',  # Chaste/mesh/src/common;
      'Cylindrical2dMesh.hpp',  # Chaste/mesh/src/mutable;
      'AbstractMeshReader.hpp',  # Chaste/mesh/src/reader;
      'AbstractPerElementWriter.hpp',  # Chaste/mesh/src/utilities;
      'VertexMesh.hpp',  # Chaste/mesh/src/vertex;
      'AbstractMeshWriter.hpp',  # Chaste/mesh/src/writer;
      'AbstractCvodeSystem.hpp',  # Chaste/ode/src/common;
      'OdeSystemForCoupledHeatEquation.hpp',  # Chaste/ode/src/fortests;
      'CvodeAdaptor.hpp',  # Chaste/ode/src/solver;
      'AbstractBoundaryCondition.hpp',  # Chaste/pde/src/common;
      'AbstractHdf5Converter.hpp',  # Chaste/pde/src/postprocessing;
      'AbstractLinearEllipticPde.hpp',  # Chaste/pde/src/problem;
      'LinearBasisFunction.hpp',  # Chaste/pde/src/solver;
      'AbstractBoundaryCondition.hpp',  # Chaste/continuum_mechanics/src/common;
      'ContinuumMechanicsProblemDefinition.hpp',  # Chaste/continuum_mechanics/src/problem;
      'AbstractCompressibleMaterialLaw.hpp',  # Chaste/continuum_mechanics/src/problem/material_laws;
      'AbstractContinuumMechanicsAssembler.hpp',  # Chaste/continuum_mechanics/src/solver
      'NonlinearElasticityTools.hpp',  # Chaste/continuum_mechanics/src/common
    ]
    self.libs = [['chaste_ode', 'chaste_continuum_mechanics', 'chaste_pde', 'chaste_io', 'chaste_mesh', 'chaste_linalg', 'chaste_global']]

    #self.extra_libs = ['lapack', 'blas']
    #self.set_rpath = False    # do not use dynamic linkage
    self.check_text = r'''
int main(int argc, char* argv[]) {
}
'''

  def check(self, ctx):
      
    # get number of available processors
    p = multiprocessing.cpu_count()

    # Setup the build handler.
    self.set_build_handler([
      #'pip install --user "python-dateutil==1.5"',
      #'easy_install --user "Amara==1.2.0.2"',    # this fails because Amara is soo old (actually the dependency 4suite-xml fails)
      #'!pip install --user "rdflib==2.4.2"',      # this requires Cython
      #'pip install --user lxml',
      "mkdir -p build; cd build; \
        METIS_ROOT=${DEPENDENCIES_DIR}/petsc/install/ \
        XSD_ROOT=${DEPENDENCIES_DIR}/xsd/install \
        SUNDIALS_ROOT=${DEPENDENCIES_DIR}/petsc/install/ \
        XERCESC_INCLUDE_DIR=${DEPENDENCIES_DIR}/xercesc/install/include \
        XERCESC_LIBRARY_DIR=${DEPENDENCIES_DIR}/xercesc/install/lib \
        HDF5_ROOT=${DEPENDENCIES_DIR}/hdf5/install \
        BOOST_ROOT=${DEPENDENCIES_DIR}/boost/install \
      cmake -DCMAKE_POLICY_DEFAULT_CMP0074=NEW \
        -DVTK_DIR=${DEPENDENCIES_DIR}/vtk/install/lib/cmake/vtk-8.2/ \
        -DCMAKE_INSTALL_PREFIX=${DEPENDENCIES_DIR}/chaste/install \
        -DChaste_USE_PETSC_HDF5=OFF \
        -DPETSC_DIR=${DEPENDENCIES_DIR}/petsc/install \
        -DPARMETIS_INCLUDE_DIR=${DEPENDENCIES_DIR}/petsc/install/include \
        -DPARMETIS_LIBRARY=${DEPENDENCIES_DIR}/petsc/install/lib/parmetis.so \
        -DMETIS_LIBRARY=${DEPENDENCIES_DIR}/petsc/install/lib/libmetis.so \
        -DLAPACK_lapack_LIBRARY=${DEPENDENCIES_DIR}/lapack/install/lib/libopenblas.so \
        -DLAPACK_LIBRARIES=${DEPENDENCIES_DIR}/lapack/install/lib/libopenblas.so \
        -DBLAS_blas_LIBRARY=${DEPENDENCIES_DIR}/lapack/install/lib/libopenblas.so \
        -DBLAS_LIBRARIES=${DEPENDENCIES_DIR}/lapack/install/lib/libopenblas.so \
        -DBUILD_SHARED_LIBS=OFF \
        -DChaste_ENABLE_TESTING=OFF \
        -DBoost_USE_STATIC_LIBS=OFF \
        -DBoost_USE_STATIC_RUNTIME=FALSE \
        -DPETSC_ARCH=arch-linux2-c-opt \
        ..",
      "!cd build; make core -j{}".format(p),   # '!' means this command is allowed to fail and fails because Amara could not be installed
      #"cd build; make install",  # this does not really work
      "mkdir -p ${PREFIX}/lib/chaste && mkdir -p ${PREFIX}/include",
      "$cp -ft ${PREFIX}/lib/chaste $(find ${SOURCE_DIR}/build/ -name *.a)",   # install manually
      "$cp -ft ${PREFIX}/include $(find ${SOURCE_DIR} -name *.hpp)",
      "$cp -ft ${PREFIX}/include $(find ${SOURCE_DIR}/mesh/src -name *.h)",
    ])
    self.number_output_lines = 1539
    self.number_output_lines = 696
    
    env = ctx.env
    ctx.Message('Checking for Chaste ...        ')
    self.check_options(env)

    res = super(Chaste, self).check(ctx)

    self.check_required(res[0], ctx)
    ctx.Result(res[0])
    return res[0]
