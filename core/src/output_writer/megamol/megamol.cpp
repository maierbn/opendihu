#include "output_writer/megamol/megamol.h"

#include <fstream>
#include <iostream>
#include <thread>
#include <chrono>

#include "easylogging++.h"
#include <utility/python_utility.h>

namespace OutputWriter
{

BoundingBox::BoundingBox():
  min(Vec3({0.0,0.0,0.0})),
  max(Vec3({0.0,0.0,0.0}))
{

}

MegaMol::MegaMol(DihuContext context, PythonConfig settings) :
  Generic(context, settings)
{
}

} // namespace
