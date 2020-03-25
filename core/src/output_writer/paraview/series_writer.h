#pragma once

#include <Python.h>  // has to be the first included header

#include <iostream>
#include <memory>
#include <string>
#include <map>

namespace OutputWriter
{

/** Writes a "*.vtk.series" JSON file that contains all written files and their timestamps.
 * This is the way to tell paraview about different files at specific time steps.
 * Example for such a file:
 * {
 *  "file-series-version" : "1.0",
 *  "files" : [
 *   { "name" : "foo1.vtk", "time" : 0 },
 *   { "name" : "foo2.vtk", "time" : 5.5 },
 *   { "name" : "foo3.vtk", "time" : 11.2 }
 *  ]
 * }
 **/
class SeriesWriter
{
public:

  //! add a new output file with timestamp to the JSON series file
  void registerNewFile(std::string filename, double currentTime);

protected:

  struct SeriesFile
  {
    std::string filename;       //< name of the series file
    std::string fileContents;   //< file contents except the last "]}"
  };

  std::map<std::string,SeriesFile> seriesFiles_;   // all series files  
};

} // namespace
