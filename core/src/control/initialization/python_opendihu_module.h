#pragma once

#include <iostream>
#include <string>
#include <Python.h>

/** This is the definition of a python module named "opendihu",
 *  which can be accessed from within the python settings as follows:
 *
 *  import opendihu
 *  print("program_name: {}".format(opendihu.program_name))
 *
 */
namespace PythonOpendihu
{

extern PyModuleDef opendihuModule;

extern std::string programName;      //< the name of the executable, as given in argv[0]
extern std::string versionText;      //< the version text of the DihuContext class
extern std::string metaText;         //< the meta text of the DihuContext class

PyMODINIT_FUNC PyInit_opendihu(void);


} // namespace PythonOpendihu

