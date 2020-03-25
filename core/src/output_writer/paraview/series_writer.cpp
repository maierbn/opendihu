#include "output_writer/paraview/series_writer.h"

#include <fstream>
#include <sstream>

#include "output_writer/generic.h"

namespace OutputWriter
{

void SeriesWriter::registerNewFile(std::string filename, double currentTime)
{
  // store new entry
  SeriesEntry entry;
  entry.filename = filename;
  entry.time = currentTime;
  seriesEntries_.push_back(entry);

  // if filename was not yet initialize
  if (filename_ == "")
  {
    filename_ = "output.vtk.series";
    fileContents_ = "{\n\t\"file-series-version\" : \"1.0\",\n\t\"files\" : [\"\n";
  }

  // add current entry to file contents
  std::stringstream s;

  // if there are already multiple entries in the file contents, add a comma
  if (seriesEntries_.size() > 1)
    s << ",\n";

  // add current entry to file contents
  s << "\t\t{ \"name\" : \"" << filename << "\", \"time\" : " << currentTime << " }";
  fileContents_ += s.str();

  // open series vtk file
  std::ofstream file;
  Generic::openFile(file, filename_, false);

  // write all entries
  file << fileContents_ << "\n\t]\n}\n";

  file.close();
}

} // namespace
