#include "output_writer/callback.h"

#include <iostream>

#include "easylogging++.h"
#include <Python.h>

#include <utility/python_utility.h>
#include <utility/petsc_utility.h>
#include <mesh/regular_fixed.h>
#include <mesh/structured_deformable.h>
#include <mesh/mesh.h>

namespace OutputWriter
{

template<typename DataType>
void Callback::writeSolution(DataType& data, int timeStepNo, double currentTime)
{
  if (!data.mesh())
  {
    LOG(FATAL) << "Mesh is not set!";
  }
  const int dimension = data.mesh()->dimension();
  
  // solution and rhs vectors in mesh shape
  switch(dimension)
  {
  case 1:
    writeSolutionDim<1>(data, timeStepNo, currentTime);
    break;
  case 2:
    writeSolutionDim<2>(data, timeStepNo, currentTime);
    break;
  case 3:
    writeSolutionDim<3>(data, timeStepNo, currentTime);
    break;
  };
}

template <int dimension, typename DataType>
void Callback::writeSolutionDim(DataType &data, int timeStepNo, double currentTime)
{
  // solution and rhs vectors in mesh shape
  // if mesh is regular fixed
  if (std::dynamic_pointer_cast<Mesh::RegularFixed<dimension>>(data.mesh()) != NULL)
  {
    // cast mesh to its special type
    std::shared_ptr<Mesh::RegularFixed<dimension>> mesh = std::dynamic_pointer_cast<Mesh::RegularFixed<dimension>>(data.mesh());
      
    // get data of solution vector
    std::vector<double> vectorValues;
    PetscUtility::getVectorEntries(data.solution(), vectorValues);
    
    // determine number of entries in each dimension
    std::vector<long int> nEntries(dimension);
    for (int i=0; i<dimension; i++)
    {
      nEntries[i] = (mesh->nElements(i) + 1) * data.nDegreesOfFreedomPerNode();
    }

    callCallback(vectorValues, nEntries, timeStepNo, currentTime);
  }
}

};