#include "output_writer/paraview/series_writer.h"

#include <fstream>
#include <sstream>

#include "output_writer/generic.h"

namespace OutputWriter
{

void SeriesWriter::registerNewFile(std::string filename, double currentTime)
{
  // remove 0 character at the end
  if (filename[filename.size()-1] == 0)
    filename = filename.substr(0,filename.size()-1);

  LOG(DEBUG) << "registerNewFile(" << filename << ", " << currentTime << ")";
  
  // extract basename of filename to store, e.g. "out/output_0000000.vtu" -> "output_0000000"
  std::string basename = StringUtility::extractBasename(filename);
  if (basename.find("_") != std::string::npos)
    basename = basename.substr(0, basename.rfind("_"));

  bool firstEntry = false;

  // if filename was not yet initialized
  if (seriesFiles_[basename].filename.empty())
  {
    // get path
    std::string path = "";
    if (filename.find("/") != std::string::npos)
    {
      path = filename.substr(0,filename.rfind("/")+1);
    }

    // get file extension
    std::string extension = ".vtk";
    if (filename.find(".") != std::string::npos)
    {
      extension = filename.substr(filename.rfind("."));
    }

    // assemble filename 
    seriesFiles_[basename].filename = 
      path + basename + extension + std::string(".series");

    // initialize file contents with the header of the file
    seriesFiles_[basename].fileContents = "{\n\t\"file-series-version\" : \"1.0\",\n\t\"files\" : [\n";
    firstEntry = true;

    // open file once to create path if necessary
    std::ofstream file;
    Generic::openFile(file, seriesFiles_[basename].filename, false);
  }

  std::stringstream s;

  // add comma to last entry if there is already an entry
  if (!firstEntry)
  {
    s << ",\n";
  }

  // remove path from file name and add filename with time to JSON file
  if (filename.find("/") != std::string::npos)
  {
    filename = filename.substr(filename.rfind("/")+1);
  }
  s << "\t\t{ \"name\" : \"" << filename << "\", \"time\" : " << currentTime << " }";
  
  seriesFiles_[basename].fileContents += s.str();
  
  // open series vtk file
  std::ofstream file;
  file.open(seriesFiles_[basename].filename.c_str(), std::ios::out);
  if (!file.is_open())
  {
    LOG(WARNING) << "Could not write to vtk series file \"" << seriesFiles_[basename].filename.c_str() << "\".";
  }
  else 
  {
    // write all entries
    file << seriesFiles_[basename].fileContents << "\n\t]\n}\n";

    file.close();
  }
}

} // namespace
