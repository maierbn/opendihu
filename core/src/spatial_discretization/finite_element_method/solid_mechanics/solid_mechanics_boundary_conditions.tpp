#include "spatial_discretization/finite_element_method/solid_mechanics/02_stiffness_matrix_incompressible.h"

#include <Python.h>  // has to be the first included header
#include "utility/python_utility.h"

namespace SpatialDiscretization
{
  
template<typename BasisOnMeshType,typename Term>
SolidMechanicsBoundaryConditions<BasisOnMeshType,Term>::TractionBoundaryCondition::
TractionBoundaryCondition(PyObject *specificSettings, std::shared_ptr<typename BasisOnMeshType::HighOrderBasisOnMesh> mesh)
{
  typedef typename BasisOnMeshType::HighOrderBasisOnMesh BasisOnMesh;  // the high order BasisOnMesh for mixed formulation or the BasisOnMesh itself for non-mixed formulation
  const int D = BasisOnMesh::dim();
 
  // store global element no
  elementGlobalNo = PythonUtility::getOptionInt(specificSettings, "element", 0, PythonUtility::NonNegative);
  
  // store face
  std::string faceStr = PythonUtility::getOptionString(specificSettings, "face", "0+");
  face = Mesh::parseFace(faceStr);
  
  if (PythonUtility::hasKey(specificSettings, "constantValue") && PythonUtility::hasKey(specificSettings, "constantVector"))
  {
    LOG(ERROR) << "Specified both \"constantValue\" and \"constantVector\".";
  }
  
  // parse dof vectors
  if (PythonUtility::hasKey(specificSettings, "constantValue") || PythonUtility::hasKey(specificSettings, "constantVector"))
  {
    VecD<D> constantVector;
    
    if (PythonUtility::hasKey(specificSettings, "constantValue"))
    {
      double constantValue = PythonUtility::getOptionDouble(specificSettings, "constantValue", 0.0);
      
      VecD<D-1> xiSurface;
      for (int i = 0; i < D-1; i++)
        xiSurface[i] = 0.5;
       
      VecD<D> xi = Mesh::getXiOnFace(face, xiSurface);
      constantVector = MathUtility::transformToD<D,3>(mesh->getNormal(face, elementGlobalNo, xi) * constantValue);
    }
    else if (PythonUtility::hasKey(specificSettings, "constantVector"))
    {
      constantVector = PythonUtility::getOptionArray<double,D>(specificSettings, "constantVector", 0.0);
    }
    
    // get dofs indices within element that correspond to the selected face
    const int D = BasisOnMesh::dim();
    const int nDofs = ::BasisOnMesh::BasisOnMeshBaseDim<D-1,typename BasisOnMesh::BasisFunction>::nDofsPerElement();
    std::array<dof_no_t,nDofs> dofIndices;
    BasisOnMesh::getFaceDofs(face, dofIndices);

    for (int i = 0; i < nDofs; i++)
    {
      dofVectors.push_back(std::pair<dof_no_t, VecD<D>>(dofIndices[i], constantVector));
    }
  }
  else if (PythonUtility::hasKey(specificSettings, "dofVectors"))
  {
    std::pair<dof_no_t, PyObject *> dofVectorItem;
    
    // loop over dofVectors
    for (dofVectorItem = PythonUtility::getOptionDictBegin<dof_no_t, PyObject *>(specificSettings, "dofVectors");
      !PythonUtility::getOptionDictEnd(specificSettings, "dofVectors"); 
      PythonUtility::getOptionDictNext<dof_no_t, PyObject *>(specificSettings, "dofVectors", dofVectorItem))
    {
      dof_no_t dofIndex = dofVectorItem.first;
      VecD<D> dofVector = PythonUtility::convertFromPython<double,D>(dofVectorItem.second);
      
      dofVectors.push_back(std::pair<dof_no_t, VecD<D>>(dofIndex, dofVector));
    }
  }
  else
  {
    LOG(ERROR) << "Traction on element has not specified \"constantValue\", \"constantVector\" nor \"dofVectors\".";
  }
}

template<typename BasisOnMeshType,typename Term>
void SolidMechanicsBoundaryConditions<BasisOnMeshType,Term>::
initializeBoundaryConditions(bool &externalVirtualWorkIsConstant, const int nUnknowns, PyObject *specificSettings, Data::FiniteElements<BasisOnMeshType,Term> &data)
{
  LOG(TRACE)<<"initializeBoundaryConditions";
 
  const int D = BasisOnMeshType::dim();
  std::shared_ptr<typename BasisOnMeshType::HighOrderBasisOnMesh> mesh = data.mesh();
  
  // ----- Neumann BC ----------
  // parse values for traction and body force
  externalVirtualWorkIsConstant = true;
  if (PythonUtility::hasKey(specificSettings, "tractionReferenceConfiguration"))
  {
    PyObject *listItem = PythonUtility::getOptionListBegin<PyObject*>(specificSettings, "tractionReferenceConfiguration"); 
    
    // loop over items in tractionReferenceConfiguration list
    for (;
         !PythonUtility::getOptionListEnd(specificSettings, "tractionReferenceConfiguration"); 
         PythonUtility::getOptionListNext<PyObject*>(specificSettings, "tractionReferenceConfiguration", listItem))
    {
      tractionReferenceConfiguration_.emplace_back(listItem, mesh);
    }
  }
  
  if (PythonUtility::hasKey(specificSettings, "tractionCurrentConfiguration"))
  {
    externalVirtualWorkIsConstant = false;
    
    PyObject *listItem = PythonUtility::getOptionListBegin<PyObject*>(specificSettings, "tractionReferenceConfiguration"); 
    
    // loop over items in tractionCurrentConfiguration list
    for (; 
         !PythonUtility::getOptionListEnd(specificSettings, "tractionCurrentConfiguration"); 
         PythonUtility::getOptionListNext<PyObject*>(specificSettings, "tractionCurrentConfiguration", listItem))
    {
      tractionCurrentConfiguration_.emplace_back(listItem, mesh);
    }
  }
  
  if (PythonUtility::hasKey(specificSettings, "bodyForceReferenceConfiguration"))
  {
   // example entry: {0: [tmax,0,0], 5: [tmax,tmax,tmax]},   # {<element global no.>: <vector>, ...}
   // loop over dict items
   std::pair<element_no_t, PyObject *> dofVectorItem;
   for (dofVectorItem = PythonUtility::getOptionDictBegin<element_no_t, PyObject *>(specificSettings, "bodyForceReferenceConfiguration");
      !PythonUtility::getOptionDictEnd(specificSettings, "bodyForceReferenceConfiguration"); 
      PythonUtility::getOptionDictNext<element_no_t, PyObject *>(specificSettings, "bodyForceReferenceConfiguration", dofVectorItem))
    {
      element_no_t elementGlobalNo = dofVectorItem.first;
      VecD<D> vector = PythonUtility::convertFromPython<VecD<D>>(dofVectorItem.second);
     
      bodyForceReferenceConfiguration_.push_back(std::pair<element_no_t,VecD<D>>(elementGlobalNo, vector));
    }
  }
  
  if (PythonUtility::hasKey(specificSettings, "bodyForceCurrentConfiguration"))
  {
    externalVirtualWorkIsConstant = false;
     
    // example entry: {0: [tmax,0,0], 5: [tmax,tmax,tmax]},   # {<global dof no.>: <vector>, ...}
    // loop over dict items
    std::pair<element_no_t, PyObject *> dofVectorItem;
    for (dofVectorItem = PythonUtility::getOptionDictBegin<element_no_t, PyObject *>(specificSettings, "bodyForceCurrentConfiguration");
       !PythonUtility::getOptionDictEnd(specificSettings, "bodyForceCurrentConfiguration");
       PythonUtility::getOptionDictNext<element_no_t, PyObject *>(specificSettings, "bodyForceCurrentConfiguration", dofVectorItem))
    {
      element_no_t elementGlobalNo = dofVectorItem.first;
      VecD<D> vector = PythonUtility::convertFromPython<VecD<D>>(dofVectorItem.second);
     
      bodyForceReferenceConfiguration_.push_back(std::pair<element_no_t,VecD<D>>(elementGlobalNo, vector));
    }
  }
  
  // --------- Dirichlet BC --------------
  // iterate over Dirichlet boundary conditions
  
  // get the first dirichlet boundary condition from the list
  std::pair<node_no_t, double> boundaryCondition 
    = PythonUtility::getOptionDictBegin<dof_no_t, double>(specificSettings, "dirichletBoundaryCondition");
  
  // loop over Dirichlet boundary conditions
  for (; !PythonUtility::getOptionDictEnd(specificSettings, "dirichletBoundaryCondition"); 
       PythonUtility::getOptionDictNext<dof_no_t, double>(specificSettings, "dirichletBoundaryCondition", boundaryCondition))
  {
    dof_no_t boundaryConditionUnknownsIndex = boundaryCondition.first;
    double boundaryConditionValue = boundaryCondition.second;
    
    if (boundaryConditionUnknownsIndex < 0)
      continue;
    
    if (boundaryConditionUnknownsIndex > nUnknowns) 
    {
      LOG(WARNING) << "Boundary condition specified for degree of freedom no. "<<boundaryConditionUnknownsIndex
       <<", but scenario has only "<<nUnknowns<<" degrees of freedom.";
       continue;
    }
    
    // store values
    dirichletIndices_.push_back(boundaryConditionUnknownsIndex);
    dirichletValues_.push_back(boundaryConditionValue);
  }
  zeros_.resize(dirichletValues_.size(), 0.0);
  
}

template<typename BasisOnMeshType,typename Term>
void SolidMechanicsBoundaryConditions<BasisOnMeshType,Term>::
printBoundaryConditions()
{
  const int D = BasisOnMeshType::dim();
  /*
   * 
  std::vector<dof_no_t> dirichletIndices_;  ///< the indices of unknowns (not dofs) for which the displacement is fixed
  std::vector<double> dirichletValues_;     ///< the to dirichletIndices corresponding fixed values for the displacement
  std::vector<double> rhsValues_;           ///< the entries in rhs for which for u Dirichlet values are specified
  
  //TODO split into boundary conditions class
  struct TractionBoundaryCondition
  {
    element_no_t elementGlobalNo;
    
    Mesh::face_t face;
    std::vector<std::pair<dof_no_t, Vec3>> dofVectors;  //<element-local dof no, value>
    
    // parse values from python config, e.g. {"element": 1, "face": "0+", "dofVectors:", {0: [tmax,0,0], 1: [tmax,0,0], 2: [tmax,0,0], 3: [tmax,0,0]}}
    TractionBoundaryCondition(PyObject *specificSettings);
  };
  
  std::vector<TractionBoundaryCondition> tractionReferenceConfiguration_;   //< tractions for elements
  std::vector<TractionBoundaryCondition> tractionCurrentConfiguration_;
  
  std::vector<std::pair<dof_no_t, Vec3>> bodyForceReferenceConfiguration_;  //< <global dof no, vector>
  std::vector<std::pair<dof_no_t, Vec3>> bodyForceCurrentConfiguration_;    //< <global dof no, vector>
  */
  LOG(DEBUG) << "============ Dirichlet BC ============";
  LOG(DEBUG) << "dirichletIndices_: " << dirichletIndices_;
  LOG(DEBUG) << "dirichletValues_: " << dirichletValues_ << std::endl;
  
  LOG(DEBUG) << "============ Neumann BC ============";
  LOG(DEBUG) << "tractionReferenceConfiguration_: size: " << tractionReferenceConfiguration_.size();
  for (int i = 0; i < tractionReferenceConfiguration_.size(); i++)
  {
    LOG(DEBUG) << "no. " << i << ", element " << tractionReferenceConfiguration_[i].elementGlobalNo 
      << ", face: " << tractionReferenceConfiguration_[i].face << ", " << tractionReferenceConfiguration_[i].dofVectors.size() << " dofVectors: ";
    for (int j = 0; j < tractionReferenceConfiguration_[i].dofVectors.size(); j++)
    {
      std::pair<dof_no_t, VecD<D>> &dofVector = tractionReferenceConfiguration_[i].dofVectors[j];
      LOG(DEBUG) << "   " << j << ". dof index: " << dofVector.first << ", vector: " << dofVector.second;
    }
  }
  
  LOG(DEBUG) << "tractionCurrentConfiguration_: size: " << tractionCurrentConfiguration_.size();
  for (int i = 0; i < tractionCurrentConfiguration_.size(); i++)
  {
    LOG(DEBUG) << "no. " << i << ", element " << tractionCurrentConfiguration_[i].elementGlobalNo 
      << ", face: " << tractionCurrentConfiguration_[i].face << ", " << tractionCurrentConfiguration_[i].dofVectors.size() << " dofVectors: ";
    for (int j = 0; j < tractionCurrentConfiguration_[i].dofVectors.size(); j++)
    {
      std::pair<dof_no_t, VecD<D>> &dofVector = tractionCurrentConfiguration_[i].dofVectors[j];
      LOG(DEBUG) << "   " << j << ". dof index: " << dofVector.first << ", vector: " << dofVector.second;
    }
  }
  LOG(DEBUG) << "bodyForceReferenceConfiguration_: size: " << bodyForceReferenceConfiguration_.size();
  
  for (int i = 0; i < bodyForceReferenceConfiguration_.size(); i++)
  {
    LOG(DEBUG) << "   " << i << ". dof global no: " << bodyForceReferenceConfiguration_[i].first
      << ", vector: " << bodyForceReferenceConfiguration_[i].second;
  }
      
  LOG(DEBUG) << "bodyForceCurrentConfiguration_: size: " << bodyForceCurrentConfiguration_.size();
  
  for (int i = 0; i < bodyForceCurrentConfiguration_.size(); i++)
  {
    LOG(DEBUG) << "   " << i << ". dof global no: " << bodyForceCurrentConfiguration_[i].first
      << ", vector: " << bodyForceCurrentConfiguration_[i].second;
  }
  LOG(DEBUG) << "============ ============";
}

template<typename BasisOnMeshType,typename Term>
void SolidMechanicsBoundaryConditions<BasisOnMeshType,Term>::
applyDirichletBoundaryConditionsInNonlinearFunction(Vec &f, Data::FiniteElements<BasisOnMeshType,Term> &data)
{
  LOG(TRACE) << "applyDirichletBoundaryConditionsInNonlinearFunction";
  
// set value for rhs for entries of Dirichlet BC
// set f(x) = x0  
  if (!data.computeWithReducedVectors())
  {
    // overwrite the values in f that correspond to entries for which Dirichlet BC are set with zero values,
    // such that the nonlinear equation is satisfied for them
    PetscErrorCode ierr;
    ierr = VecSetValues(f, this->dirichletIndices_.size(), this->dirichletIndices_.data(), this->zeros_.data(), INSERT_VALUES); CHKERRV(ierr);
  }
  
// set f(x) to f(x) = x - x0
#if 0  
  if (!this->data_.computeWithReducedVectors())
  {
    PetscErrorCode ierr;
    
    // get displacements values for indices for with Dirichlet BC are set
    static std::vector<double> displacementValues;
    displacementValues.resize(dirichletIndices_.size());
    ierr = VecGetValues(this->data_.displacements().values(), dirichletIndices_.size(), dirichletIndices_.data(), displacementValues.data()); CHKERRV(ierr);
    
    VLOG(1) << "  got displacement: " << displacementValues;
    
    // subtract dirichletValues
    for (int i = 0; i < displacementValues.size(); i++)
    {
      displacementValues[i] -= dirichletValues_[i];
    }
    
    VLOG(1) << "  subtracted dirichlet: " << displacementValues;
    
    // set f(x) = x - x0 for dofs with dirichlet values x0
    ierr = VecSetValues(f, dirichletIndices_.size(), dirichletIndices_.data(), displacementValues.data(), INSERT_VALUES); CHKERRV(ierr);
  }
#endif  
}

template<typename BasisOnMeshType,typename Term>
void SolidMechanicsBoundaryConditions<BasisOnMeshType,Term>::
applyDirichletBoundaryConditionsInDisplacements(Data::FiniteElements<BasisOnMeshType,Term> &data)
{
  LOG(TRACE) << "applyDirichletBoundaryConditionsInDisplacements";
 
  PetscErrorCode ierr;
  
  // if the nonlinear solver should work with not reduced vectors
  if (false)
  if (!data.computeWithReducedVectors())
  {
    // set entries of Dirichlet BCs to specified values
    ierr = VecSetValues(data.displacements().values(), this->dirichletIndices_.size(), this->dirichletIndices_.data(), this->dirichletValues_.data(), INSERT_VALUES); CHKERRV(ierr);
    
    VLOG(1) << "after applying Dirichlet BC displacements u:" << PetscUtility::getStringVector(data.displacements().values());
  }
}

template<typename BasisOnMeshType,typename Term>
void SolidMechanicsBoundaryConditions<BasisOnMeshType,Term>::
applyDirichletBoundaryConditionsInStiffnessMatrix(Mat &matrix, Data::FiniteElements<BasisOnMeshType,Term> &data)
{
  // if the nonlinear solver should work with not reduced vectors
  if (!data.computeWithReducedVectors())
  {
   
    VLOG(1) << "stiffness matrix before MatZeroRowsColumns: " << PetscUtility::getStringMatrix(matrix);
   
    // zero rows and columns for which Dirichlet BC is set 
    PetscErrorCode ierr;
    ierr = MatZeroRowsColumns(matrix, this->dirichletIndices_.size(), this->dirichletIndices_.data(), 1.0, PETSC_IGNORE, PETSC_IGNORE); CHKERRV(ierr);
    
  }
}

template<typename BasisOnMeshType,typename Term>
void SolidMechanicsBoundaryConditions<BasisOnMeshType,Term>::
reduceVector(Vec &input, Vec &output, const int nUnknownsInputVector)
{ 
  const double *inputData;
  double *outputData;
  VecGetArrayRead(input, &inputData);
  VecGetArray(output, &outputData);
  
  dof_no_t reducedIndex = 0;
  std::vector<dof_no_t>::const_iterator dirichletIndicesIter = dirichletIndices_.begin();
  
  for (dof_no_t currentDofNo = 0; currentDofNo < nUnknownsInputVector; currentDofNo++)
  {
    // exclude variables for which Dirichlet BC are set
    if (dirichletIndicesIter != dirichletIndices_.end())
    {
      if (currentDofNo == *dirichletIndicesIter)
      {
        dirichletIndicesIter++;
        continue;
      }
    }
    
    outputData[reducedIndex++] = inputData[currentDofNo];
  }
  
  VecRestoreArrayRead(input, &inputData);
  VecRestoreArray(output, &outputData);
}

template<typename BasisOnMeshType,typename Term>
void SolidMechanicsBoundaryConditions<BasisOnMeshType,Term>::
expandVector(Vec &input, Vec &output, const int nUnknownsOutputVector)
{
  const double *inputData;
  double *outputData;
  VecGetArrayRead(input, &inputData);
  VecGetArray(output, &outputData);
  
  dof_no_t reducedIndex = 0;
  std::vector<dof_no_t>::const_iterator dirichletIndicesIter = dirichletIndices_.begin();
  std::vector<double>::const_iterator dirichletValuesIter = dirichletValues_.begin();
  
  for (dof_no_t currentDofNo = 0; currentDofNo < nUnknownsOutputVector; currentDofNo++)
  {
    // exclude variables for which Dirichlet BC are set
    if (dirichletIndicesIter != dirichletIndices_.end())
    {
      if (currentDofNo == *dirichletIndicesIter)
      {
        outputData[currentDofNo] = *dirichletValuesIter;
        dirichletIndicesIter++;
        dirichletValuesIter++;
        
        continue;
      }
    }
    
    outputData[currentDofNo] = inputData[reducedIndex++];
  }
  
  VecRestoreArrayRead(input, &inputData);
  VecRestoreArray(output, &outputData);
}

template<typename BasisOnMeshType,typename Term>
void SolidMechanicsBoundaryConditions<BasisOnMeshType,Term>::
expandVectorTo3D(Vec &input, Vec &output, const int nUnknowns3D)
{
  // compute number of reduced dofs
  assert(BasisOnMeshType::dim() == 2);
  
  const double *inputData;
  double *outputData;
  VecGetArrayRead(input, &inputData);
  VecGetArray(output, &outputData);
  
  dof_no_t reducedIndex = 0;
  for (dof_no_t currentDofNo = 0; currentDofNo < nUnknowns3D; currentDofNo++)
  {
    if (currentDofNo % 3 == 2)
    {
      outputData[currentDofNo] = 0.0;
    }
    else 
    {
      outputData[currentDofNo] = inputData[reducedIndex++];
    }
  }
  
  VecRestoreArrayRead(input, &inputData);
  VecRestoreArray(output, &outputData);
}

template<typename BasisOnMeshType,typename Term>
void SolidMechanicsBoundaryConditions<BasisOnMeshType,Term>::
reduceMatrix(Mat &input, Mat &output, const int nUnknownsFull)
{
  // compute number of reduced dofs
  int nUnknownsReduced = nUnknownsFull - this->dirichletValues_.size();
  
  std::vector<int> notConstraintIndices;
  notConstraintIndices.reserve(nUnknownsReduced);
  
  std::vector<dof_no_t>::const_iterator dirichletIndicesIter = dirichletIndices_.begin();
  
  for (dof_no_t currentDofNo = 0; currentDofNo < nUnknownsFull; currentDofNo++)
  {
    // exclude variables for which Dirichlet BC are set
    if (dirichletIndicesIter != dirichletIndices_.end())
    {
      if (currentDofNo == *dirichletIndicesIter)
      {
        dirichletIndicesIter++;
        continue;
      }
    }
    
    notConstraintIndices.push_back(currentDofNo);
  }
  
  // create Petsc IS object of the rows and columns to keep
  IS isrow;
  PetscErrorCode ierr;
  ierr = ISCreateGeneral(PETSC_COMM_WORLD, nUnknownsReduced, notConstraintIndices.data(), PETSC_COPY_VALUES, &isrow); CHKERRV(ierr);
  
  VLOG(1) << "non-reduced tangent stiffness matrix: " << PetscUtility::getStringMatrix(input);
  VLOG(1) << "rows to extract: " << notConstraintIndices << " target matrix object: " << output;
  
  // extract reduced matrix to output
  if (output == PETSC_NULL)
  {
    VLOG(1) << "MatGetSubMatrix MAT_INITIAL_MATRIX";
    ierr = MatGetSubMatrix(input, isrow, isrow, MAT_INITIAL_MATRIX, &output); CHKERRV(ierr);
    VLOG(1) << "output: " << output;
  }
  else
  {
    VLOG(1) << "MatGetSubMatrix MAT_REUSE_MATRIX";
    ierr = MatGetSubMatrix(input, isrow, isrow, MAT_REUSE_MATRIX, &output);// CHKERRV(ierr);
    if (ierr != 0)
    {
      LOG(WARNING) << "Could not reuse previous matrix data structure";
      VLOG(1) << "MatGetSubMatrix MAT_INITIAL_MATRIX";
      
      ierr = MatGetSubMatrix(input, isrow, isrow, MAT_INITIAL_MATRIX, &output); CHKERRV(ierr);
    }
    VLOG(1) << "output: " << output;
  }
  
  VLOG(1) << "reduced tangent stiffness matrix: " << PetscUtility::getStringMatrix(output);
}
  
}; // namespace