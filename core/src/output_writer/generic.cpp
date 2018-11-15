#include "output_writer/generic.h"

#include <chrono>
#include <thread>

namespace OutputWriter
{

Generic::Generic(DihuContext context, PythonConfig specificSettings) :
  context_(context), specificSettings_(specificSettings)
{
  // get the rank subset of all processes that collectively call the write methods
  rankSubset_ = this->context_.partitionManager()->rankSubsetForCollectiveOperations();
  VLOG(1) << "OutputWriter::Generic constructor, rankSubset: " << *rankSubset_;
}

Generic::~Generic()
{
}

std::ofstream Generic::openFile(std::string filename, bool append)
{
  // open file
  std::ofstream file;
  if (append)
  {
    file.open(filename.c_str(), std::ios::out | std::ios::binary | std::ios::app);
  }
  else
  {
    file.open(filename.c_str(), std::ios::out | std::ios::binary);
  }

  if (!file.is_open())
  {
    // try to create directories
    if (filename.rfind("/") != std::string::npos)
    {
      // extract directory from filename
      std::string path = filename.substr(0, filename.rfind("/"));

      // create directory and wait until system has created it
      int ret = system((std::string("mkdir -p ")+path).c_str());
      if (ret != 0)
        LOG(WARNING) << "Creation of directory \"" <<path<< "\" failed.";
      std::this_thread::sleep_for (std::chrono::milliseconds(500));

      file.clear();
      file.open(filename.c_str(), std::ios::out | std::ios::binary);
    }
  }

  if (!file.is_open())
  {
    LOG(WARNING) << "Could not open file \"" <<filename << "\" for writing!";
  }

  return file;
}

void Generic::appendRankNo(std::stringstream &str, int nRanks, int ownRankNo)
{
  int nCharacters = 1 + int(std::log10(nRanks));

  str << "." << std::setw(nCharacters) << std::setfill('0') << ownRankNo;
}

};
