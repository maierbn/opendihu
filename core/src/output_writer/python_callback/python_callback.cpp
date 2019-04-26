#include "output_writer/python_callback/python_callback.h"

#include <Python.h>  // has to be the first included header
#include <iostream>

#include "easylogging++.h"
#include "utility/python_utility.h"

namespace OutputWriter
{

PythonCallback::PythonCallback(DihuContext context, PythonConfig settings, std::shared_ptr<Partition::RankSubset> rankSubset) :
  Generic(context, settings, rankSubset)
{
  callback_ = settings.getOptionPyObject("callback");
  onlyNodalValues_ = settings.getOptionBool("onlyNodalValues", true);
}

}  // namespace
