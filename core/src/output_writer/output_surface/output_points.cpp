#include "output_writer/output_surface/output_points.h"

#include "output_writer/generic.h"
#include "output_writer/paraview/paraview.h"

#include <algorithm>

namespace OutputWriter
{

void OutputPoints::
writeCsvFile(std::string filename, double currentTime, const std::vector<double> &geometry, const std::vector<double> &values, bool writeGeometry)
{
  std::ofstream file;
  Generic::openFile(file, filename, true);  // append to file

  int nPointsGlobal = geometry.size()/3;

  // in first timestep, write header
  if (currentTime <= 1e-5)
  {
    if (!writeGeometry)
    {
      file << "#electrode positions (x0,y0,z0,x1,y1,z1,...);\n#; ";
      for (int i = 0; i < nPointsGlobal; i++)
      {
        file << ";" << geometry[3*i+0]
          << ";" << geometry[3*i+1]
          << ";" << geometry[3*i+2];
      }
      file << std::endl;
    }

    file << "#timestamp;t;n_points";

    if (writeGeometry)
    {
      for (int pointNo = 0; pointNo < nPointsGlobal; pointNo++)
      {
        file << ";p" << pointNo << "_x;p" << pointNo << "_y;p" << pointNo << "_z";
      }
    }

    if (!values.empty())
    {
      for (int pointNo = 0; pointNo < nPointsGlobal; pointNo++)
      {
        file << ";p" << pointNo << "_value";
      }
    }
    file << std::endl;
  }

  // write timestamp and time
  // time stamp
  auto t = std::time(nullptr);
  auto tm = *std::localtime(&t);
  file << StringUtility::timeToString(&tm) << ";"
    << currentTime << ";" << nPointsGlobal;

  // write geometry
  if (writeGeometry)
  {
    for (int i = 0; i < nPointsGlobal; i++)
    {
      file << ";" << geometry[3*i+0]
        << ";" << geometry[3*i+1]
        << ";" << geometry[3*i+2];
    }
  }

  // write values
  if (!values.empty())
  {
    for (int i = 0; i < nPointsGlobal; i++)
    {
      file << ";" << values[i];
    }
  }
  file << std::endl;
  file.close();
}

void OutputPoints::
writeVtpFile(std::string filename, double currentTime, const std::vector<double> &geometry,
             const std::vector<double> &values, int nComponents,
             const std::vector<int> &partitioning, std::string fieldVariableName)
{
  int nPointsGlobal = geometry.size()/3;

  std::ofstream file;
  Generic::openFile(file, filename, false);  // do not append to existing file

  // transform current time to string
  std::vector<double> time(1, currentTime);
  std::string stringTime = Paraview::encodeBase64Float(time.begin(), time.end());

  file << "<?xml version=\"1.0\"?>" << std::endl
    << "<!-- " << DihuContext::versionText() << " " << DihuContext::metaText()
    << ", currentTime: " << currentTime << " -->" << std::endl
    << "<VTKFile type=\"PolyData\" version=\"1.0\" byte_order=\"LittleEndian\">" << std::endl    // intel cpus are LittleEndian
    << std::string(1, '\t') << "<PolyData>" << std::endl
    << std::string(2, '\t') << "<FieldData>" << std::endl
    << std::string(3, '\t') << "<DataArray type=\"Float32\" Name=\"Time\" NumberOfTuples=\"1\" format=\"binary\" >" << std::endl
    << std::string(4, '\t') << stringTime << std::endl
    << std::string(3, '\t') << "</DataArray>" << std::endl
    << std::string(2, '\t') << "</FieldData>" << std::endl;

  file << std::string(2, '\t') << "<Piece NumberOfPoints=\"" << nPointsGlobal << "\" NumberOfVerts=\"0\" "
    << "NumberOfLines=\"0\" NumberOfStrips=\"0\" NumberOfPolys=\"0\">" << std::endl
    << std::string(3, '\t') << "<PointData>" << std::endl;


  // write emg values
  if (!values.empty())
  {
    file << std::string(4, '\t') << "<DataArray "
        << "Name=\"" << fieldVariableName << "\" "
        << "type=\"Float32\" "
        << "NumberOfComponents=\"" << nComponents << "\" "
        << "format=\"binary\" >" << std::endl << std::string(5, '\t')
        << Paraview::encodeBase64Float(values.begin(), values.end())
        << std::endl << std::string(4, '\t') << "</DataArray>" << std::endl;
  }

  // write partitioning values
  file << std::string(4, '\t') << "<DataArray "
      << "Name=\"partitioning\" "
      << "type=\"Int32\" "
      << "NumberOfComponents=\"1\" "
      << "format=\"binary\" >" << std::endl << std::string(5, '\t')
      << Paraview::encodeBase64Int32(partitioning.begin(), partitioning.end())
      << std::endl << std::string(4, '\t') << "</DataArray>" << std::endl;

  file << std::string(3, '\t') << "</PointData>" << std::endl
    << std::string(3, '\t') << "<CellData>" << std::endl
    << std::string(3, '\t') << "</CellData>" << std::endl
    << std::string(3, '\t') << "<Points>" << std::endl
    << std::string(4, '\t') << "<DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"binary"
    << "\" >" << std::endl << std::string(5, '\t');

  // write geometry values
  file << Paraview::encodeBase64Float(geometry.begin(), geometry.end());

  file << std::endl << std::string(4, '\t') << "</DataArray>" << std::endl
    << std::string(3, '\t') << "</Points>" << std::endl
    << std::string(3, '\t') << "<Verts></Verts>" << std::endl
    << std::string(3, '\t') << "<Lines></Lines>" << std::endl
    << std::string(3, '\t') << "<Strips></Strips>" << std::endl
    << std::string(3, '\t') << "<Polys></Polys>" << std::endl
    << std::string(2, '\t') << "</Piece>" << std::endl
    << std::string(1, '\t') << "</PolyData>" << std::endl
    << "</VTKFile>" << std::endl;

  file.close();
}

}  // namespace OutputWriter
