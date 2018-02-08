#include "output_writer/exfile/exfile.h"

#include <fstream>
#include <iostream>
#include <thread>
#include <chrono>

#include "easylogging++.h"

#include <utility/python_utility.h>
#include <utility/petsc_utility.h>
#include <mesh/regular_fixed.h>
#include <mesh/structured_deformable.h>
#include <mesh/unstructured_deformable.h>
#include <mesh/mesh.h>

namespace OutputWriter
{

Exfile::Exfile(PyObject *settings) : Generic(settings)
{
}


};