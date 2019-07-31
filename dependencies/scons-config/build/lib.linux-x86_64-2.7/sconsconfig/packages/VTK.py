import sys, os, multiprocessing
from Package import Package
import subprocess

##
##  Handles VTK library
##
class VTK(Package):

  def __init__(self, **kwargs):
    defaults = {
        'download_url': 'https://www.vtk.org/files/release/8.2/VTK-8.2.0.zip',
    }
    defaults.update(kwargs)
    super(VTK, self).__init__(**defaults)
    self.ext = '.cpp'
    #self.sub_dirs = [
    #    ('include/mysql', 'lib'),
    #    ('include/mysql', 'lib64'),
    #]
    self.headers = ['vtkLine.h']
    self.libs = [['vtkViewsInfovis-8.2', 'vtkChartsCore-8.2', 'vtkIOExportOpenGL2-8.2', 'vtkImagingStencil-8.2', 'vtkGeovisCore-8.2', 'vtkDomainsChemistryOpenGL2-8.2', 'vtkCommonColor-8.2', 'vtkIOExportPDF-8.2', 'vtkIOExport-8.2', 'vtkRenderingContextOpenGL2-8.2', 'vtkViewsContext2D-8.2', 'vtkRenderingContext2D-8.2', 'vtkCommonComputationalGeometry-8.2', 'vtkIOInfovis-8.2', 'vtkInfovisCore-8.2', 'vtkDomainsChemistry-8.2', 'vtkIOMINC-8.2', 'vtkIOAMR-8.2', 'vtkFiltersParallelImaging-8.2', 'vtkFiltersGeneric-8.2', 'vtkFiltersHybrid-8.2', 'vtkFiltersAMR-8.2', 'vtkIOParallel-8.2', 'vtkFiltersParallel-8.2', 'vtkFiltersExtraction-8.2', 'vtkIOAsynchronous-8.2', 'vtkInteractionImage-8.2', 'vtkIOCityGML-8.2', 'vtkIOExodus-8.2', 'vtkIOGeometry-8.2', 'vtkIOImport-8.2', 'vtkIOLSDyna-8.2', 'vtkIONetCDF-8.2', 'vtkIOParallelXML-8.2', 'vtkIOPLY-8.2', 'vtkIOSQL-8.2', 'vtkIOVideo-8.2', 'vtksys-8.2', 'vtkRenderingGL2PSOpenGL2-8.2', 'vtkRenderingImage-8.2', 'vtkRenderingLabel-8.2', 'vtkRenderingLOD-8.2', 'vtkRenderingVolumeOpenGL2-8.2', 'vtkRenderingVolume-8.2', 'vtkRenderingCore-8.2', 'vtkFiltersPoints-8.2', 'vtkFiltersModeling-8.2', 'vtkFiltersSMP-8.2', 'vtkFiltersTexture-8.2', 'vtkCommonTransforms-8.2', 'vtkFiltersGeneral-8.2', 'vtkFiltersFlowPaths-8.2', 'vtkFiltersHyperTree-8.2', 'vtkFiltersImaging-8.2', 'vtkFiltersProgrammable-8.2', 'vtkFiltersSelection-8.2', 'vtkFiltersTopology-8.2', 'vtkFiltersVerdict-8.2', 'vtkImagingColor-8.2', 'vtkImagingHybrid-8.2', 'vtkImagingMath-8.2', 'vtkImagingMorphological-8.2', 'vtkImagingStatistics-8.2', 'vtkIOEnSight-8.2', 'vtkIOMovie-8.2', 'vtkIOSegY-8.2', 'vtkIOTecplotTable-8.2', 'vtkIOVeraOut-8.2', 'vtkCommonExecutionModel-8.2', 'vtkCommonDataModel-8.2', 'vtkCommonCore-8.2', 'vtksys-8.2', 'vtkCommonSystem-8.2', 'vtkCommonMisc-8.2', 'vtkCommonTransforms-8.2', 'vtkCommonMath-8.2', 'vtkCommonCore-8.2', 'vtkDICOMParser-8.2', 'vtksys-8.2', 'vtkIOXMLParser-8.2', 'vtkFiltersSources-8.2', 'vtkFiltersCore-8.2', 'vtkexodusII-8.2', 'vtkdoubleconversion-8.2', 'vtkRenderingOpenGL2-8.2', 'vtkexpat-8.2', 'vtkNetCDF-8.2', 'vtkParallelCore-8.2', 'vtkCommonDataModel-8.2', 'vtkCommonExecutionModel-8.2', 'vtkCommonCore-8.2', 'vtkCommonMisc-8.2', 'vtkCommonTransforms-8.2', 'vtkCommonSystem-8.2', 'vtkCommonMath-8.2', 'vtkFiltersStatistics-8.2', 'vtkFiltersGeometry-8.2', 'vtklz4-8.2', 'vtkIOCore-8.2', 'vtkCommonComputationalGeometry-8.2', 'vtkCommonCore-8.2', 'vtkCommonDataModel-8.2', 'vtkCommonExecutionModel-8.2', 'vtkFiltersCore-8.2', 'vtkImagingSources-8.2', 'vtkImagingGeneral-8.2', 'vtkIOLegacy-8.2', 'vtkImagingCore-8.2', 'vtkFiltersGeneral-8.2', 'vtkImagingFourier-8.2', 'vtkCommonMisc-8.2', 'vtkfreetype-8.2', 'vtkverdict-8.2', 'vtkgl2ps-8.2', 'vtkhdf5_hl-8.2', 'vtkhdf5-8.2', 'vtkzlib-8.2', 'vtkIOXML-8.2', 'vtkViewsCore-8.2', 'vtkInteractionStyle-8.2', 'vtkInfovisLayout-8.2', 'vtkInteractionWidgets-8.2', 'vtkproj-8.2', 'vtkglew-8.2', 'vtkpng-8.2', 'vtkCommonCore-8.2', 'vtkCommonExecutionModel-8.2', 'vtkCommonDataModel-8.2', 'vtkCommonMath-8.2', 'vtkCommonTransforms-8.2', 'vtkImagingCore-8.2', 'vtksys-8.2', 'vtkImagingSources-8.2', 'vtkCommonCore-8.2', 'vtkCommonExecutionModel-8.2', 'vtkImagingCore-8.2', 'vtkCommonDataModel-8.2', 'vtkInfovisCore-8.2', 'vtkImagingHybrid-8.2', 'vtkFiltersCore-8.2', 'vtkFiltersSources-8.2', 'vtkCommonSystem-8.2', 'vtkFiltersGeneral-8.2', 'vtkCommonComputationalGeometry-8.2', 'vtkRenderingCore-8.2', 'vtkFiltersExtraction-8.2', 'vtkRenderingAnnotation-8.2', 'vtkFiltersModeling-8.2', 'vtkFiltersHybrid-8.2', 'vtkImagingGeneral-8.2', 'vtkRenderingFreeType-8.2', 'vtkIOImage-8.2', 'vtkpugixml-8.2', 'vtklz4-8.2', 'vtklzma-8.2', 'vtkdoubleconversion-8.2', 'vtklibharu-8.2', 'vtkCommonCore-8.2', 'vtkCommonExecutionModel-8.2', 'vtkCommonDataModel-8.2', 'vtksys-8.2', 'vtkDICOMParser-8.2', 'vtkCommonSystem-8.2', 'vtkCommonTransforms-8.2', 'vtkmetaio-8.2', 'vtkCommonMath-8.2', 'vtkpng-8.2', 'vtktiff-8.2', 'vtkjpeg-8.2', 'vtkzlib-8.2', 'vtklibxml2-8.2', 'vtkIOCore-8.2', 'vtktheora-8.2', 'vtkogg-8.2', 'vtkjsoncpp-8.2', 'vtksqlite-8.2', 'vtkIOXMLParser-8.2', 'vtkCommonMisc-8.2', 'vtkCommonCore-8.2', 'vtkIOCore-8.2', 'vtkCommonDataModel-8.2', 'vtkexpat-8.2', 'vtksys-8.2', 'vtkzlib-8.2', 'vtkRenderingCore-8.2', 'vtkFiltersSources-8.2', 'vtkFiltersCore-8.2', 'vtkFiltersGeneral-8.2', 'vtkCommonColor-8.2', 'vtkFiltersGeometry-8.2', 'vtkCommonComputationalGeometry-8.2', 'vtkCommonCore-8.2', 'vtkCommonExecutionModel-8.2', 'vtkCommonSystem-8.2', 'vtkCommonDataModel-8.2', 'vtkCommonMath-8.2', 'vtkCommonTransforms-8.2', 'vtksys-8.2', 'vtkCommonMisc-8.2', 'vtkfreetype-8.2']]
    #self.extra_libs = ['lapack', 'blas']
    #self.set_rpath = False    # do not use dynamic linkage
    self.check_text = r'''
// source: https://lorensen.github.io/VTKExamples/site/Cxx/SimpleOperations/DistancePointToLine/
#include <vtkSmartPointer.h>
#include <vtkLine.h>
#include <vtkPoints.h>

int main(int argc, char* argv[]) {

  double lineP0[3] = {0.0, 0.0, 0.0};
  double lineP1[3] = {2.0, 0.0, 0.0};

  double p0[3] = {1.0, 0, 0};
  double p1[3] = {1.0, 2.0, 0};

  /*
  vtkSmartPointer<vtkLine> line = 
    vtkSmartPointer<vtkLine>::New();
  line->GetPoints()->SetPoint(0, lineP0);
  line->GetPoints()->SetPoint(0, lineP1);
  */

  /*{
  double dist0 = vtkLine::DistanceToLine(p0, lineP0, lineP1);
  std::cout << "Dist0: " << dist0 << std::endl;

  double dist1 = vtkLine::DistanceToLine(p1, lineP0, lineP1);
  std::cout << "Dist1: " << dist1 << std::endl;
  }

  {
  double t;
  double closest[3];
  double dist0 = vtkLine::DistanceToLine(p0, lineP0, lineP1, t, closest);
  std::cout << "Dist0: " << dist0 << " closest point: " << closest[0] << " " << closest[1] << " " << closest[2] << std::endl;

  double dist1 = vtkLine::DistanceToLine(p1, lineP0, lineP1, t, closest);
  std::cout << "Dist1: " << dist1 << " closest point: " << closest[0] << " " << closest[1] << " " << closest[2] << std::endl;
  } */
}
'''

  def check(self, ctx):
      
    # get number of available processors
    p = multiprocessing.cpu_count()

    # Setup the build handler.
    self.set_build_handler([
      "mkdir -p build; cd build; cmake -DCMAKE_INSTALL_PREFIX=${PREFIX} -DBUILD_SHARED_LIBS=ON -DCMAKE_BUILD_TYPE=Release ..",
      "cd build; make -j{}".format(p),
      "cd build; make install",
    ])
    self.number_output_lines = 4768
    
    env = ctx.env
    ctx.Message('Checking for VTK ...           ')
    self.check_options(env)

    res = super(VTK, self).check(ctx)

    self.check_required(res[0], ctx)
    ctx.Result(res[0])
    return res[0]
