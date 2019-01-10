#include "output_writer/megamol/megamol.h"

#include <fstream>
#include <iostream>
#include <thread>
#include <chrono>

#include "easylogging++.h"
#include <utility/python_utility.h>

namespace OutputWriter
{

MegaMOL::MegaMOL(DihuContext context, PythonConfig settings) :
  Generic(context, settings)
{
}

};  // namespace
