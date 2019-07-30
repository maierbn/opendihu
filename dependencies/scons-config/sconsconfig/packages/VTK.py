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
    self.libs = [
      'vtkChartsCore-8.2.so',
      'vtkCommonColor-8.2.so',
      'vtkCommonComputationalGeometry-8.2.so',
      'vtkCommonCore-8.2.so',
      'vtkCommonDataModel-8.2.so',
      'vtkCommonExecutionModel-8.2.so',
      'vtkCommonMath-8.2.so',
      'vtkCommonMisc-8.2.so',
      'vtkCommonSystem-8.2.so',
      'vtkCommonTransforms-8.2.so',
      'vtkDICOMParser-8.2.so',
      'vtkDomainsChemistry-8.2.so',
      'vtkDomainsChemistryOpenGL2-8.2.so',
      'vtkdoubleconversion-8.2.so',
      'vtkexodusII-8.2.so',
      'vtkexpat-8.2.so',
      'vtkFiltersAMR-8.2.so',
      'vtkFiltersCore-8.2.so',
      'vtkFiltersExtraction-8.2.so',
      'vtkFiltersFlowPaths-8.2.so',
      'vtkFiltersGeneral-8.2.so',
      'vtkFiltersGeneric-8.2.so',
      'vtkFiltersGeometry-8.2.so',
      'vtkFiltersHybrid-8.2.so',
      'vtkFiltersHyperTree-8.2.so',
      'vtkFiltersImaging-8.2.so',
      'vtkFiltersModeling-8.2.so',
      'vtkFiltersParallel-8.2.so',
      'vtkFiltersParallelImaging-8.2.so',
      'vtkFiltersPoints-8.2.so',
      'vtkFiltersProgrammable-8.2.so',
      'vtkFiltersSelection-8.2.so',
      'vtkFiltersSMP-8.2.so',
      'vtkFiltersSources-8.2.so',
      'vtkFiltersStatistics-8.2.so',
      'vtkFiltersTexture-8.2.so',
      'vtkFiltersTopology-8.2.so',
      'vtkFiltersVerdict-8.2.so',
      'vtkfreetype-8.2.so',
      'vtkGeovisCore-8.2.so',
      'vtkgl2ps-8.2.so',
      'vtkglew-8.2.so',
      'vtkhdf5-8.2.so',
      'vtkhdf5_hl-8.2.so',
      'vtkImagingColor-8.2.so',
      'vtkImagingCore-8.2.so',
      'vtkImagingFourier-8.2.so',
      'vtkImagingGeneral-8.2.so',
      'vtkImagingHybrid-8.2.so',
      'vtkImagingMath-8.2.so',
      'vtkImagingMorphological-8.2.so',
      'vtkImagingSources-8.2.so',
      'vtkImagingStatistics-8.2.so',
      'vtkImagingStencil-8.2.so',
      'vtkInfovisCore-8.2.so',
      'vtkInfovisLayout-8.2.so',
      'vtkInteractionImage-8.2.so',
      'vtkInteractionStyle-8.2.so',
      'vtkInteractionWidgets-8.2.so',
      'vtkIOAMR-8.2.so',
      'vtkIOAsynchronous-8.2.so',
      'vtkIOCityGML-8.2.so',
      'vtkIOCore-8.2.so',
      'vtkIOEnSight-8.2.so',
      'vtkIOExodus-8.2.so',
      'vtkIOExport-8.2.so',
      'vtkIOExportOpenGL2-8.2.so',
      'vtkIOExportPDF-8.2.so',
      'vtkIOGeometry-8.2.so',
      'vtkIOImage-8.2.so',
      'vtkIOImport-8.2.so',
      'vtkIOInfovis-8.2.so',
      'vtkIOLegacy-8.2.so',
      'vtkIOLSDyna-8.2.so',
      'vtkIOMINC-8.2.so',
      'vtkIOMovie-8.2.so',
      'vtkIONetCDF-8.2.so',
      'vtkIOParallel-8.2.so',
      'vtkIOParallelXML-8.2.so',
      'vtkIOPLY-8.2.so',
      'vtkIOSegY-8.2.so',
      'vtkIOSQL-8.2.so',
      'vtkIOTecplotTable-8.2.so',
      'vtkIOVeraOut-8.2.so',
      'vtkIOVideo-8.2.so',
      'vtkIOXML-8.2.so',
      'vtkIOXMLParser-8.2.so',
      'vtkjpeg-8.2.so',
      'vtkjsoncpp-8.2.so',
      'vtklibharu-8.2.so',
      'vtklibxml2-8.2.so',
      'vtklz4-8.2.so',
      'vtklzma-8.2.so',
      'vtkmetaio-8.2.so',
      'vtkNetCDF-8.2.so',
      'vtkogg-8.2.so',
      'vtkParallelCore-8.2.so',
      'vtkpng-8.2.so',
      'vtkproj-8.2.so',
      'vtkpugixml-8.2.so',
      'vtkRenderingAnnotation-8.2.so',
      'vtkRenderingContext2D-8.2.so',
      'vtkRenderingContextOpenGL2-8.2.so',
      'vtkRenderingCore-8.2.so',
      'vtkRenderingFreeType-8.2.so',
      'vtkRenderingGL2PSOpenGL2-8.2.so',
      'vtkRenderingImage-8.2.so',
      'vtkRenderingLabel-8.2.so',
      'vtkRenderingLOD-8.2.so',
      'vtkRenderingOpenGL2-8.2.so',
      'vtkRenderingVolume-8.2.so',
      'vtkRenderingVolumeOpenGL2-8.2.so',
      'vtksqlite-8.2.so',
      'vtksys-8.2.so',
      'vtktheora-8.2.so',
      'vtktiff-8.2.so',
      'vtkverdict-8.2.so',
      'vtkViewsContext2D-8.2.so',
      'vtkViewsCore-8.2.so',
      'vtkViewsInfovis-8.2.so',
      'vtkzlib-8.2.so'
    ]
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
      "mkdir -p build; cd build; cmake -DCMAKE_INSTALL_PREFIX=${PREFIX} ..",
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
