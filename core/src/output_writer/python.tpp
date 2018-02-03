#include "output_writer/python.h"

#include <fstream>
#include <iostream>
#include <thread>
#include <chrono>

#include "easylogging++.h"
#include <Python.h>

#include <utility/python_utility.h>
#include <utility/petsc_utility.h>
#include <mesh/regular_fixed.h>
#include <mesh/mesh.h>

namespace OutputWriter
{

template<typename DataType>
void Python::write(DataType& data, int timeStepNo, double currentTime)
{
  // check if output should be written in this timestep and prepare filename
  if (!Generic::prepareWrite(data, timeStepNo, currentTime))
  {
    return;
  }
  
  const int dimension = data.mesh()->dimension();
  
  // solution and rhs vectors in mesh shape
  switch(dimension)
  {
  case 1:
    writeSolutionDim<1>(data);
    break;
  case 2:
    writeSolutionDim<2>(data);
    break;
  case 3:
    writeSolutionDim<3>(data);
    break;
  };
}

template <int dimension, typename DataType>
void Python::writeSolutionDim(DataType &data)
{
  LOG(TRACE) << "writeSolution<"<<dimension<<">()";
  
  // solution and rhs vectors in mesh shape
  // if mesh is regular fixed
  if (std::dynamic_pointer_cast<Mesh::RegularFixed<dimension>>(data.mesh()) != NULL)
  {
    // cast mesh to its special type
    std::shared_ptr<Mesh::RegularFixed<dimension>> mesh = std::dynamic_pointer_cast<Mesh::RegularFixed<dimension>>(data.mesh());
      
    // determine file names
    std::stringstream s[2];
    s[0] << filename_ << "_solution.npy";
    s[1] << filename_ << "_solution_shaped.npy";
    std::string filenameSolution = s[0].str();
    std::string filenameSolutionShaped = s[1].str();
    
    // get data of solution vector
    std::vector<double> vectorValues;
    PetscUtility::getVectorEntries(data.solution().values(), vectorValues);
    
    // determine number of entries in each dimension
    std::vector<long int> nEntries(dimension);
    for (int i=0; i<dimension; i++)
    {
      nEntries[i] = (mesh->nElements(i) + 1) * data.nComponentsPerNode();
    }
    std::vector<long int> singleEntry({(long)vectorValues.size()});
    
    // write as numpy file
    writeToNumpyFile(vectorValues, filenameSolution, 1, singleEntry);
    writeToNumpyFile(vectorValues, filenameSolutionShaped, dimension, nEntries);
    
    // if data is of class Data::FiniteElements
    if (dynamic_cast<Data::FiniteElements<typename DataType::BasisOnMesh> *>(&data) != NULL)
    {
      writeRhsMatrix<dimension>(*dynamic_cast<Data::FiniteElements<typename DataType::BasisOnMesh> *>(&data));
    }
  }
}

template <int dimension, typename DataType>
void Python::writeRhsMatrix(DataType &data)
{
  // solution and rhs vectors in mesh shape
  if (std::dynamic_pointer_cast<Mesh::RegularFixed<dimension>>(data.mesh()) != NULL)
  {
    // cast mesh to its special type
    std::shared_ptr<Mesh::RegularFixed<dimension>> mesh = std::dynamic_pointer_cast<Mesh::RegularFixed<dimension>>(data.mesh());
      
    // determine file names
    std::stringstream s[3];
    s[0] << filename_ << "_rhs.npy";
    s[1] << filename_ << "_rhs_shaped.npy";
    s[2] << filename_ << "_stiffness.npy";
    std::string filenameRhs = s[0].str();
    std::string filenameRhsShaped = s[1].str();
    std::string filenameStiffness = s[2].str();
    
    // get data of rhs vector
    int vectorSize = 0;
    VecGetSize(data.rightHandSide().values(), &vectorSize);

    std::vector<int> indices(vectorSize);
    std::iota(indices.begin(), indices.end(), 0);
    std::vector<double> vectorValues(vectorSize);

    // determine number of entries in each dimension
    std::vector<long int> nEntries(dimension);
    for (int i=0; i<dimension; i++)
    {
      nEntries[i] = (mesh->nElements(i) + 1);
    }
    std::vector<long int> singleEntry({(long)vectorValues.size()});
    
    VecGetValues(data.rightHandSide().values(), vectorSize, indices.data(), vectorValues.data());
    
    // write as numpy file
    writeToNumpyFile(vectorValues, filenameRhs, 1, singleEntry);
    writeToNumpyFile(vectorValues, filenameRhsShaped, dimension, nEntries);
    
    // get stiffness matrix
    int nRows, nColumns;
    MatGetSize(data.stiffnessMatrix(), &nRows, &nColumns);
    std::vector<int> rowIndices(nRows);
    std::iota(rowIndices.begin(), rowIndices.end(), 0);
    std::vector<int> columnIndices(nColumns);
    std::iota(columnIndices.begin(), columnIndices.end(), 0);
    std::vector<double> matrixValues(nRows*nColumns);
    
    nEntries = {nRows, nColumns};
    
    MatGetValues(data.stiffnessMatrix(), nRows, rowIndices.data(), nColumns, columnIndices.data(), matrixValues.data());
    
    // write as numpy file
    writeToNumpyFile(matrixValues, filenameStiffness, 2, nEntries);
    
  }
}

}; // namespace