#pragma once

#include <Python.h>  // has to be the first included header
#include <iostream>
#include <vector>

#include "control/types.h"
#include "output_writer/generic.h"
#include "output_writer/paraview/poly_data_properties_for_mesh.h"

namespace OutputWriter
{

class Paraview : public Generic
{
public:

  //! constructor
  Paraview(DihuContext context, PythonConfig specificSettings, std::shared_ptr<Partition::RankSubset> rankSubset = nullptr);

  //! write out solution to given filename, if timeStepNo is not -1, this value will be part of the filename
  template<typename DataType>
  void write(DataType &data, int timeStepNo = -1, double currentTime = -1, int callCountIncrement = 1);

  //! write the given field variable as VTK <DataArray> element to file, if onlyParallelDatasetElement write the <PDataArray> element
  template<typename FieldVariableType>
  static void writeParaviewFieldVariable(FieldVariableType &fieldVariable, std::ofstream &file,
                                         bool binaryOutput, bool fixedFormat, bool onlyParallelDatasetElement=false);


  //! write the a field variable indicating which ranks own which portion of the domain as VTK <DataArray> element to file, if onlyParallelDatasetElement write the <PDataArray> element
  template<typename FieldVariableType>
  static void writeParaviewPartitionFieldVariable(FieldVariableType &geometryField, std::ofstream &file,
                                                  bool binaryOutput, bool fixedFormat, bool onlyParallelDatasetElement=false);
  
  //! write a single *.vtp file that contains all data of all 1D field variables. This is uses MPI IO. It can be enabled with the "combineFiles" option.
  //! on return, combinedMeshesOut contains the 1D mesh names that were written to the vtp file.
  template<typename FieldVariablesForOutputWriterType>
  void writePolyDataFile(const FieldVariablesForOutputWriterType &fieldVariables, std::set<std::string> &combinedMeshesOut);


  //! write a single *.vtu file that contains all data of all 3D or 2D field variables (depending on output3DMeshes).
  //! This is uses MPI IO. It can be enabled with the "combineFiles" option.
  //! on return, combinedMeshesOut contains the 3D or 2D mesh names that were written to the vtu file.
  template<typename FieldVariablesForOutputWriterType>
  void writeCombinedUnstructuredGridFile(const FieldVariablesForOutputWriterType &fieldVariables, std::set<std::string> &combinedMeshesOut,
                                         bool output3DMeshes);

  //! encode a Petsc vector in Base64,
  //! @param withEncodedSizePrefix if the length of the vector should be added as encoded prefix
  static std::string encodeBase64Vec(const Vec &vector, bool withEncodedSizePrefix=true);

  //! encode a std::vector<double> or std::list<double> as base64
  template <typename Iter>
  static std::string encodeBase64Float(Iter iterBegin, Iter iterEnd, bool withEncodedSizePrefix=true);

  //! encode a std::vector<int> as 32bit base64 values
  template <typename Iter>
  static std::string encodeBase64Int32(Iter iterBegin, Iter iterEnd, bool withEncodedSizePrefix=true);

  //! encode a vector as UInt64 values
  template <typename Iter>
  std::string encodeBase64UInt8(Iter iterBegin, Iter iterEnd, bool withEncodedSizePrefix=true);

  //! convert to a string with space separated values
  static std::string convertToAscii(const Vec &vector, bool humanReadable);
  
  //! convert to a string with space separated values
  static std::string convertToAscii(const std::vector<double> &vector, bool humanReadable);
  
  //! convert to a string with space separated values
  static std::string convertToAscii(const std::vector<element_no_t> &vector, bool humanReadable);
  
  // only if Petsc uses long long int for PetscInt, element_no_t is not the same as int and we need to define convertToAscii for int vectors
#ifdef PETSC_USE_64BIT_INDICES
  //! convert to a string with space separated values
  static std::string convertToAscii(const std::vector<int> &vector, bool humanReadable);
#endif

protected:

  /** one VTKPiece is the XML element that will be output as <Piece></Piece>. It is created from one or multiple opendihu meshes
   */
  struct VTKPiece
  {
    std::set<std::string> meshNamesCombinedMeshes;   ///< the meshNames of the combined meshes, or only one meshName if it is not a merged mesh
    PolyDataPropertiesForMesh properties;   ///< the properties of the merged mesh

    std::string firstScalarName;   ///< name of the first scalar field variable of the mesh
    std::string firstVectorName;   ///< name of the first non-scalar field variable of the mesh

    //! constructor, initialize nPoints and nCells to 0
    VTKPiece();

    //! assign the correct values to firstScalarName and firstVectorName, only if properties has been set
    void setVTKValues();
  };

  //! write some ascii data to the file as a collective shared operation. Only rank 0 writes, but the other ranks wait and the shared file pointer is incremented.
  void writeAsciiDataShared(MPI_File fileHandle, int ownRankNo, std::string writeBuffer);

  //! write the values vector combined to the file, correctly encoded, identifier is an id to access cached values
  template<typename T>
  void writeCombinedValuesVector(MPI_File fileHandle, int ownRankNo, const std::vector<T> &values, int identifier, bool writeFloatsAsInt=false);

  //! write a vector containing nValues "12" (if output3DMeshes) or "9" (if !output3DMeshes) values for the types for an unstructured grid
  void writeCombinedTypesVector(MPI_File fileHandle, int ownRankNo, int nValues, bool output3DMeshes, int identifier);

  bool binaryOutput_;  ///< if the data output should be binary encoded using base64
  bool fixedFormat_;   ///< if non-binary output is selected, if the ascii values should be written with a fixed precision, like 1.000000e5

  bool combineFiles_;   ///< if the output data should be combined for 1D meshes into a single PolyData output file (*.vtp) and for 2D and 3D meshes to normal *.vtu,*.vts or *.vtr files. This is needed when the number of output files should be reduced.

  std::vector<int> globalValuesSize_;   ///< cached values used in writeCombinedValuesVector
  std::vector<int> nPreviousValues_;    ///< cached values used in writeCombinedValuesVector

  std::map<std::string, PolyDataPropertiesForMesh> meshPropertiesUnstructuredGridFile2D_;    ///< mesh information for a combined unstructured grid file (*.vtu), for 2D data
  std::map<std::string, PolyDataPropertiesForMesh> meshPropertiesUnstructuredGridFile3D_;    ///< mesh information for a combined unstructured grid file (*.vtu), for 3D data
  std::map<std::string, PolyDataPropertiesForMesh> meshPropertiesPolyDataFile_;    ///< mesh information for a poly data file (*.vtp), for 1D data
  VTKPiece vtkPiece_;   ///< the VTKPiece data structure used for PolyDataFile

  int nCellsPreviousRanks1D_ = 0;   ///< sum of number of cells on other processes with lower rank no., for vtp file
  int nPointsPreviousRanks1D_ = 0;  ///< sum of number of points on other processes with lower rank no., for vtp file
  int nPointsGlobal1D_ = 0;       ///< total number of points on all ranks, for vtp file
  int nLinesGlobal1D_ = 0;       ///< total number of lines on all ranks, for vtp file
  int nCellsPreviousRanks3D_ = 0;   ///< sum of number of cells on other processes with lower rank no., for vtu file
  int nPointsPreviousRanks3D_ = 0;  ///< sum of number of points on other processes with lower rank no., for vtu file
  int nPointsGlobal3D_ = 0;       ///< total number of points on all ranks, for vtu file
};

} // namespace

#include "output_writer/paraview/paraview.tpp"
#include "output_writer/paraview/paraview_write_combined_values.tpp"
#include "output_writer/paraview/paraview_write_combined_file_1D.tpp"
#include "output_writer/paraview/paraview_write_combined_file_2D3D.tpp"
