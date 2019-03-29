#include "spatial_discretization/boundary_conditions/neumann_boundary_conditions.h"

#include "easylogging++.h"
#include "utility/python_utility.h"
#include "utility/vector_operators.h"
#include "control/types.h"

namespace SpatialDiscretization
{

// 2D,3D initializeRhs
template<typename FunctionSpaceType,typename QuadratureType,int nComponents>
void NeumannBoundaryConditionsInitializeRhs<FunctionSpaceType,QuadratureType,nComponents>::
initializeRhs()
{
  LOG(TRACE) << "NeumannBoundaryConditions::initializeRhs, D=" << FunctionSpaceType::dim() << ", nComponents=" << nComponents;
  // initialize RHS for mesh dimension 2 or 3, this is the same for nComponents == 1 and nComponents > 1

  this->data_.rhs()->setRepresentationGlobal();
  this->data_.rhs()->zeroEntries();
  this->data_.rhs()->startGhostManipulation();
  this->data_.rhs()->zeroGhostBuffer();

  typedef typename NeumannBoundaryConditionsBase<FunctionSpaceType,QuadratureType,nComponents>::ElementWithFaces ElementWithFaces;
  typedef FunctionSpace::FunctionSpace<typename FunctionSpaceType::SurfaceMesh, typename FunctionSpaceType::BasisFunction> FunctionSpaceSurface;

  std::shared_ptr<FunctionSpaceType> functionSpace = this->data_.functionSpace();

  const int D = FunctionSpaceType::dim();  // = 2 or 3
  const int nDofsPerElement = FunctionSpaceType::nDofsPerElement();

  // define shortcuts for quadrature
  typedef Quadrature::TensorProduct<D-1,QuadratureType> QuadratureSurface;

  // define type to hold evaluations of integrand for result vector
  typedef std::array<double, nDofsPerElement*nComponents> EvaluationsType;
  typedef std::array<
            EvaluationsType,
            QuadratureSurface::numberEvaluations()
          > EvaluationsArraySurfaceType;     // evaluations[nGP^D](nDofs*D)

  std::array<std::array<double,D-1>, QuadratureSurface::numberEvaluations()> samplingPointsSurface = QuadratureSurface::samplingPoints();
  EvaluationsArraySurfaceType evaluationsArraySurface;

  // assemble rhs contribution
  // loop over elements that have boundary conditions on their faces
  for (typename std::vector<ElementWithFaces>::const_iterator elementIter = this->boundaryConditionElements_.begin();
       elementIter != this->boundaryConditionElements_.end(); elementIter++)
  {
    // integrate over face

    //  element_no_t elementNoLocal;
    //
    //  Mesh::face_t face;
    //  std::vector<std::pair<dof_no_t, VecD>> dofVectors;  //<element-local dof no, value>

    element_no_t elementNoLocal = elementIter->elementNoLocal;

    LOG(DEBUG) << "element no " << elementNoLocal << ", dofVectors: " << elementIter->dofVectors;

    // check if element no is valid
    if (elementNoLocal < 0 || elementNoLocal > functionSpace->nElementsLocal())
    {
      LOG(ERROR) << "Element local no. " << elementNoLocal << " for which Neumann BC is specified, is invalid (number of local elements: " << functionSpace->nElementsLocal() << ")";
      continue;
    }

    std::array<Vec3,FunctionSpaceType::nDofsPerElement()> geometryVolume;
    std::array<Vec3,FunctionSpaceSurface::nDofsPerElement()> geometrySurface;
    functionSpace->getElementGeometry(elementNoLocal, geometryVolume);
    functionSpace->extractSurfaceGeometry(geometryVolume, elementIter->face, geometrySurface);

    // loop over integration points (e.g. gauss points)
    for (unsigned int samplingPointIndex = 0; samplingPointIndex < samplingPointsSurface.size(); samplingPointIndex++)
    {
      // get parameter values of current sampling point
      std::array<double,D-1> xiSurface = samplingPointsSurface[samplingPointIndex];
      VecD<D> xi = Mesh::getXiOnFace(elementIter->face, xiSurface);

      // compute the 3xD jacobian of the parameter space to world space mapping
      std::array<Vec3,D-1> jacobian = FunctionSpaceSurface::computeJacobian(geometrySurface, xiSurface);
      double integrationFactor = MathUtility::computeIntegrationFactor<D-1>(jacobian);

      // set all entries to 0
      evaluationsArraySurface[samplingPointIndex] = {0.0};

      // compute the value at xi by summing contributions from elemental dofs, weighted with the D-dimensional ansatz functions at xi
      // this corresponds to f(xi) = sum_i psi_i(xi) * f_i
      VecD<nComponents> boundaryConditionValueAtXi({0.0});
      for (typename std::vector<std::pair<dof_no_t, VecD<nComponents>>>::const_iterator dofVectorsIter = elementIter->dofVectors.begin();
           dofVectorsIter != elementIter->dofVectors.end();
           dofVectorsIter++)
      {
        int dofIndex = dofVectorsIter->first;
        VecD<nComponents> fluxValue = dofVectorsIter->second;

        boundaryConditionValueAtXi += fluxValue * functionSpace->phi(dofIndex, xi);
      }

      // loop over dofs of element with given boundary condition value
      for (typename std::vector<std::pair<dof_no_t, VecD<nComponents>>>::const_iterator dofVectorsIter = elementIter->dofVectors.begin();
           dofVectorsIter != elementIter->dofVectors.end();
           dofVectorsIter++)
      {
        int dofIndex = dofVectorsIter->first;

        VecD<nComponents> dofIntegrand = boundaryConditionValueAtXi * functionSpace->phi(dofIndex, xi) * integrationFactor;

        LOG(DEBUG) << "  dofIndex " << dofIndex << ", xi=" << xi << ", BC value: " << boundaryConditionValueAtXi
          << " phi = " << functionSpace->phi(dofIndex, xi);

        // store integrand in evaluations array
        for (int i = 0; i < nComponents; i++)
        {
          evaluationsArraySurface[samplingPointIndex][dofIndex*nComponents + i] = dofIntegrand[i];
        }
      }  // dofIndex

    }  // function evaluations

    // integrate all values for result vector entries at once
    EvaluationsType integratedValues = QuadratureSurface::computeIntegral(evaluationsArraySurface);

    // get indices of element-local dofs
    std::array<dof_no_t,nDofsPerElement> dofNosLocal = functionSpace->getElementDofNosLocal(elementNoLocal);

    // add entries in result vector
    // loop over indices of unknows (dofIndex,dofComponent)
    for (int dofIndex = 0; dofIndex < nDofsPerElement; dofIndex++)
    {
      std::array<double,nComponents> dofValues;
      std::copy(integratedValues.begin() + dofIndex*nComponents, integratedValues.begin() + (dofIndex+1)*nComponents, dofValues.begin());

      //! set a single dof (all components) , after all calls to setValue(s), finishGhostManipulation has to be called to apply the cached changes
      //void setValue(dof_no_t dofLocalNo, const std::array<double,nComponents> &value, InsertMode petscInsertMode=INSERT_VALUES);
      this->data_.rhs()->setValue(dofNosLocal[dofIndex], dofValues, ADD_VALUES);

      //LOG(DEBUG) << "rhs: " << *this->data_.rhs();

#ifndef NDEBUG
      for (int dofComponent = 0; dofComponent < nComponents; dofComponent++)
      {
        // compute index of degree of freedom and component (integrade values vector index)
        const int i = dofIndex*nComponents + dofComponent;

        // get integrated value
        double integratedValue = integratedValues[i];

        LOG(DEBUG) << "  dof " << dofIndex << ", component " << dofComponent << " integrated value: " << integratedValue;
      }  // dofComponent
#endif
    }  // dofIndex
  } // elementGlobalNo

  this->data_.rhs()->finishGhostManipulation();
}

// 1D initializeRhs
template<typename FunctionSpaceType, typename QuadratureType, int nComponents>
void NeumannBoundaryConditions<FunctionSpaceType, QuadratureType, nComponents, Mesh::isDim<1,typename FunctionSpaceType::Mesh>>::
initializeRhs()
{
  LOG(TRACE) << "initializeRhs 1D, nComponents: " << nComponents;
  // initialize RHS for mesh dimension 1

  this->data_.rhs()->setRepresentationGlobal();
  this->data_.rhs()->zeroEntries();
  this->data_.rhs()->startGhostManipulation();
  this->data_.rhs()->zeroGhostBuffer();

  typedef typename NeumannBoundaryConditionsBase<FunctionSpaceType,QuadratureType,1>::ElementWithFaces ElementWithFaces;

  std::shared_ptr<FunctionSpaceType> functionSpace = this->data_.functionSpace();
  const int nDofsPerElement = FunctionSpaceType::nDofsPerElement();

  // assemble rhs contribution
  // loop over elements that have boundary conditions on their faces
  for (typename std::vector<ElementWithFaces>::const_iterator elementIter = this->boundaryConditionElements_.begin();
       elementIter != this->boundaryConditionElements_.end(); elementIter++)
  {
    // integrate over face

    //  element_no_t elementNoLocal;
    //
    //  Mesh::face_t face;
    //  std::vector<std::pair<dof_no_t, VecD>> dofVectors;  //<element-local dof no, value>

    element_no_t elementNoLocal = elementIter->elementNoLocal;

    LOG(DEBUG) << "element no " << elementNoLocal << ", dofVectors: " << elementIter->dofVectors[0];

    // check if element no is valid
    if (elementNoLocal < 0 || elementNoLocal > functionSpace->nElementsLocal())
    {
      LOG(ERROR) << "Element local no. " << elementNoLocal << " for which Neumann BC is specified, is invalid (number of local elements: " << functionSpace->nElementsLocal() << ")";
      continue;
    }

    double sign = 1;
    if (elementIter->face == Mesh::face_t::face0Minus)
      sign = -1;

    assert(!elementIter->dofVectors.empty());
    VecD<1> boundaryConditionValue = elementIter->dofVectors[0].second * sign;

    // get indices of element-local dofs
    std::array<dof_no_t,nDofsPerElement> dofNosLocal = functionSpace->getElementDofNosLocal(elementNoLocal);

    // add entries in result vector
    this->data_.rhs()->setValue(dofNosLocal[elementIter->dofVectors[0].first], boundaryConditionValue, ADD_VALUES);

  } // elementGlobalNo

  this->data_.rhs()->finishGhostManipulation();
}

// 2D,3D, nComponents > 1
template<typename FunctionSpaceType,typename QuadratureType,int nComponents,typename DummyForTraits>
typename NeumannBoundaryConditions<FunctionSpaceType,QuadratureType,nComponents,DummyForTraits>::ElementWithFaces
NeumannBoundaryConditions<FunctionSpaceType,QuadratureType,nComponents,DummyForTraits>::
parseElementWithFaces(PythonConfig specificSettings, std::shared_ptr<FunctionSpaceType> functionSpace)
{
  // extract one item of the Neumann BC specification, e.g. {"element": 1, "face": "0+", "dofVectors:", {0: [tmax,0,0], 1: [tmax,0,0], 2: [tmax,0,0], 3: [tmax,0,0]}}
  // nComponents == FunctionSpaceType::dim() for traction boundary condition or nComponents = 1 for flux BC

  ElementWithFaces result;

  // store global element no
  result.elementNoLocal = specificSettings.getOptionInt("element", 0, PythonUtility::NonNegative);
  bool inputMeshIsGlobal = specificSettings.getOptionBool("inputMeshIsGlobal", true);

  // this is either the local or the global element no, depending on "inputMeshIsGlobal"

  // store face
  std::string faceStr = specificSettings.getOptionString("face", "0+");
  result.face = Mesh::parseFace(faceStr);

  if (specificSettings.hasKey("constantValue") && specificSettings.hasKey("constantVector"))
  {
    LOG(ERROR) << "Specified both \"constantValue\" and \"constantVector\".";
  }

  // parse dof vectors
  // set the constant vector if it is present
  if (specificSettings.hasKey("constantValue") || specificSettings.hasKey("constantVector"))
  {
    VecD<nComponents> constantVector;

    if (specificSettings.hasKey("constantValue"))
    {
      double constantValue = specificSettings.getOptionDouble("constantValue", 0.0);

      if (nComponents == 1)
      {
        // for normal flux BC
        constantVector[0] = constantValue;
      }
      else
      {
        // for traction BC
        const int D = nComponents;

        VecD<D-1> xiSurface;
        for (int i = 0; i < D-1; i++)
          xiSurface[i] = 0.5;

        VecD<D> xi = Mesh::getXiOnFace(result.face, xiSurface);
        if (inputMeshIsGlobal)
        {
          LOG(ERROR) << "For \"inputMeshIsGlobal\" the \"constantValue\" is not implemented for nComponents=" << nComponents << " (!=1).";
        }
        constantVector = MathUtility::transformToD<D,3>(functionSpace->getNormal(result.face, result.elementNoLocal, xi) * constantValue);
      }
    }
    else if (specificSettings.hasKey("constantVector"))
    {
      constantVector = specificSettings.getOptionArray<double,nComponents>("constantVector", 0.0);
    }

    // get dofs indices within element that correspond to the selected face
    const int D = FunctionSpaceType::dim();
    const int nDofs = FunctionSpace::FunctionSpaceBaseDim<D-1,typename FunctionSpaceType::BasisFunction>::nDofsPerElement();
    std::array<dof_no_t,nDofs> dofIndices;
    FunctionSpaceType::getFaceDofs(result.face, dofIndices);

    int stride = 1;
    if (std::is_same<typename FunctionSpaceType::BasisFunction, BasisFunction::Hermite>::value)
      stride = 2;

    for (int i = 0; i < nDofs; i+=stride)
    {
      result.dofVectors.push_back(std::pair<dof_no_t, VecD<nComponents>>(dofIndices[i], constantVector));
    }
  }
  else if (specificSettings.hasKey("dofVectors"))
  {
    std::pair<dof_no_t, PyObject *> dofVectorItem;

    // loop over dofVectors
    for (dofVectorItem = specificSettings.getOptionDictBegin<dof_no_t, PyObject *>("dofVectors");
      !specificSettings.getOptionDictEnd("dofVectors");
      specificSettings.getOptionDictNext<dof_no_t, PyObject *>("dofVectors", dofVectorItem))
    {
      dof_no_t dofIndex = dofVectorItem.first;
      VecD<nComponents> dofVector = PythonUtility::convertFromPython<std::array<double,nComponents>>::get(dofVectorItem.second);

      result.dofVectors.push_back(std::pair<dof_no_t, VecD<nComponents>>(dofIndex, dofVector));
    }
  }
  else
  {
    LOG(ERROR) << "Neumann boundary condition on element " << result.elementNoLocal << " has neither specified \"constantValue\", \"constantVector\" nor \"dofVectors\".";
  }
  return result;
}

// 2D,3D, nComponents == 1
template<typename FunctionSpaceType,typename QuadratureType>
typename NeumannBoundaryConditions<FunctionSpaceType,QuadratureType,1,Mesh::isNotDim<1,typename FunctionSpaceType::Mesh>>::ElementWithFaces
NeumannBoundaryConditions<FunctionSpaceType,QuadratureType,1,Mesh::isNotDim<1,typename FunctionSpaceType::Mesh>>::
parseElementWithFaces(PythonConfig specificSettings, std::shared_ptr<FunctionSpaceType> functionSpace)
{
  // extract one item of the Neumann BC specification, e.g. {"element": 1, "face": "0+", "dofVectors:", {0: [tmax,0,0], 1: [tmax,0,0], 2: [tmax,0,0], 3: [tmax,0,0]}}
  // nComponents = 1 for flux BC

  const int D = FunctionSpaceType::dim();
  //LOG(TRACE) << "NeumannBoundaryConditions::parseElementWithFaces, D=" << D << ", nComponents=1";

  ElementWithFaces result;

  // store global element no
  result.elementNoLocal = specificSettings.getOptionInt("element", 0, PythonUtility::NonNegative);

  // this is either the local or the global element no, depending on "inputMeshIsGlobal"

  // store face
  std::string faceStr = specificSettings.getOptionString("face", "0+");
  result.face = Mesh::parseFace(faceStr);

  if (specificSettings.hasKey("constantValue") && specificSettings.hasKey("constantVector"))
  {
    LOG(ERROR) << "Specified both \"constantValue\" and \"constantVector\".";
  }

  // parse dof vectors
  // set the constant vector if it is present
  if (specificSettings.hasKey("constantValue") || specificSettings.hasKey("constantVector"))
  {
    VecD<1> constantVector;

    if (specificSettings.hasKey("constantValue"))
    {
      double constantValue = specificSettings.getOptionDouble("constantValue", 0.0);

      // determine directionFactor
      // if the surface normal is facing in opposite direction of the normal coordinate direction, this is -1, else 1
      VecD<D-1> xiSurface;
      for (int i = 0; i < D-1; i++)
        xiSurface[i] = 0.5;

      VecD<D> xi = Mesh::getXiOnFace(result.face, xiSurface);

      LOG(DEBUG) << "face: " << Mesh::getString(result.face) << ", elementNoLocal: " << result.elementNoLocal << ", xi: " << xi << ", xiSurface: " << xiSurface;

      Vec3 normal = functionSpace->getNormal(result.face, result.elementNoLocal, xi);
      double directionFactor = normal[(int)(result.face)/2];
      VLOG(1) << "normal: " << normal << ", directionFactor: " << directionFactor;

      // for normal flux BC
      constantVector[0] = constantValue * directionFactor;
    }
    else if (specificSettings.hasKey("constantVector"))
    {
      constantVector = specificSettings.getOptionArray<double,1>("constantVector", 0.0);
    }

    // get dofs indices within element that correspond to the selected face
    const int D = FunctionSpaceType::dim();
    const int nDofs = FunctionSpace::FunctionSpaceBaseDim<D-1,typename FunctionSpaceType::BasisFunction>::nDofsPerElement();
    std::array<dof_no_t,nDofs> dofIndices;
    FunctionSpaceType::getFaceDofs(result.face, dofIndices);

    LOG(DEBUG) << "nDofs on " << D-1 << "D face: " << nDofs << ": " << dofIndices;

    int stride = 1;
    if (std::is_same<typename FunctionSpaceType::BasisFunction, BasisFunction::Hermite>::value)
      stride = FunctionSpaceType::nDofsPerNode();

    LOG(DEBUG) << "stride: " << stride;

    for (int i = 0; i < nDofs; i+=stride)
    {
      result.dofVectors.push_back(std::pair<dof_no_t, VecD<1>>(dofIndices[i], constantVector));
    }
  }
  else if (specificSettings.hasKey("dofVectors"))
  {
    std::pair<dof_no_t, PyObject *> dofVectorItem;

    // loop over dofVectors
    for (dofVectorItem = specificSettings.getOptionDictBegin<dof_no_t, PyObject *>("dofVectors");
      !specificSettings.getOptionDictEnd("dofVectors");
      specificSettings.getOptionDictNext<dof_no_t, PyObject *>("dofVectors", dofVectorItem))
    {
      dof_no_t dofIndex = dofVectorItem.first;
      VecD<1> dofVector = PythonUtility::convertFromPython<std::array<double,1>>::get(dofVectorItem.second);

      result.dofVectors.push_back(std::pair<dof_no_t, VecD<1>>(dofIndex, dofVector));
    }
  }
  else
  {
    LOG(ERROR) << "Neumann boundary condition on element " << result.elementNoLocal << " has neither specified \"constantValue\", \"constantVector\" nor \"dofVectors\".";
  }
  return result;
}

// 1D
template<typename FunctionSpaceType, typename QuadratureType, int nComponents>
typename NeumannBoundaryConditions<FunctionSpaceType,QuadratureType,nComponents,Mesh::isDim<1,typename FunctionSpaceType::Mesh>>::ElementWithFaces
NeumannBoundaryConditions<FunctionSpaceType, QuadratureType, nComponents, Mesh::isDim<1,typename FunctionSpaceType::Mesh>>::
parseElementWithFaces(PythonConfig specificSettings, std::shared_ptr<FunctionSpaceType> functionSpace)
{
  // extract one item of the Neumann BC specification, e.g. {"element": 1, "face": "0+", "dofVectors:", {0: [tmax,0,0], 1: [tmax,0,0], 2: [tmax,0,0], 3: [tmax,0,0]}}
  // nComponents == FunctionSpaceType::dim() for traction boundary condition or nComponents = 1 for flux BC

  ElementWithFaces result;

  // store global element no
  result.elementNoLocal = specificSettings.getOptionInt("element", 0, PythonUtility::NonNegative);

  // this is either the local or the global element no, depending on "inputMeshIsGlobal"

  // store face
  std::string faceStr = specificSettings.getOptionString("face", "0+");
  result.face = Mesh::parseFace(faceStr);

  if (specificSettings.hasKey("constantValue") && specificSettings.hasKey("constantVector"))
  {
    LOG(ERROR) << "Specified both \"constantValue\" and \"constantVector\".";
  }

  // parse dof vectors
  // set the constant vector if it is present
  if (specificSettings.hasKey("constantValue") || specificSettings.hasKey("constantVector"))
  {
    VecD<1> constantVector;

    if (specificSettings.hasKey("constantValue"))
    {
      constantVector[0] = specificSettings.getOptionDouble("constantValue", 0.0);
    }
    else if (specificSettings.hasKey("constantVector"))
    {
      constantVector = specificSettings.getOptionArray<double,1>("constantVector", 0.0);
    }

    // get dofs indices within element that correspond to the selected face
    const int nDofs = FunctionSpace::FunctionSpaceBaseDim<1,typename FunctionSpaceType::BasisFunction>::nDofsPerNode();
    std::array<dof_no_t,nDofs> dofIndices;
    FunctionSpaceType::getFaceDofs(result.face, dofIndices);

    int stride = 1;
    if (std::is_same<typename FunctionSpaceType::BasisFunction, BasisFunction::Hermite>::value)
      stride = 2;

    for (int i = 0; i < nDofs; i+=stride)
    {
      result.dofVectors.push_back(std::pair<dof_no_t, VecD<1>>(dofIndices[i], constantVector));
    }
  }
  else if (specificSettings.hasKey("dofVectors"))
  {
    std::pair<dof_no_t, PyObject *> dofVectorItem;

    // loop over dofVectors
    for (dofVectorItem = specificSettings.getOptionDictBegin<dof_no_t, PyObject *>("dofVectors");
      !specificSettings.getOptionDictEnd("dofVectors");
      specificSettings.getOptionDictNext<dof_no_t, PyObject *>("dofVectors", dofVectorItem))
    {
      dof_no_t dofIndex = dofVectorItem.first;
      VecD<1> dofVector = PythonUtility::convertFromPython<std::array<double,1>>::get(dofVectorItem.second);

      result.dofVectors.push_back(std::pair<dof_no_t, VecD<1>>(dofIndex, dofVector));
    }
  }
  else
  {
    LOG(ERROR) << "Neumann boundary condition on element " << result.elementNoLocal << " has neither specified \"constantValue\", \"constantVector\" nor \"dofVectors\".";
  }
  return result;
}

}  // namespace
