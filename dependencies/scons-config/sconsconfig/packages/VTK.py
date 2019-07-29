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
    #self.libs = ['mysqlclient']
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
    ctx.Message('Checking for VTK ...        ')
    self.check_options(env)

    res = super(VTK, self).check(ctx)

    self.check_required(res[0], ctx)
    ctx.Result(res[0])
    return res[0]
