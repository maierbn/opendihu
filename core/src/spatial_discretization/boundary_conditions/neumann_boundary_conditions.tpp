#include "spatial_discretization/boundary_conditions/neumann_boundary_conditions.h"

#include "easylogging++.h"
#include "utility/python_utility.h"
#include "utility/vector_operators.h"
#include "control/types.h"
#include "quadrature/gauss.h"

namespace SpatialDiscretization
{

// 2D,3D initializeRhs
template<typename FunctionSpaceType,typename QuadratureType,int nComponents>
void NeumannBoundaryConditionsInitializeRhs<FunctionSpaceType,QuadratureType,nComponents>::
initializeRhs()
{
  LOG(TRACE) << "NeumannBoundaryConditions::initializeRhs, D=" << FunctionSpaceType::dim() << ", nComponents=" << nComponents;
  // initialize RHS for mesh dimension 2 or 3, this is the same for nComponents == 1 and nComponents > 1

/*
  LOG(DEBUG) << "test";

  const int nComponentsTest = 1;

  std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,nComponentsTest>> test = this->data_.functionSpace()->template createFieldVariable<nComponentsTest>("test");
  test->zeroEntries();
  test->zeroGhostBuffer();
  test->setRepresentationGlobal();
  test->startGhostManipulation();
  LOG(DEBUG) << "inititalized, test: " << *test;


  int ownRankNo = this->data_.functionSpace()->meshPartition()->rankSubset()->ownRankNo();
  if (ownRankNo == 0)
  {
    std::array<double,nComponentsTest> value({10.0});
    test->setValue(0, value, INSERT_VALUES);
  }
  else
  {
    std::array<double,nComponentsTest> value({11.0});
    test->setValue(0, value, INSERT_VALUES);
  }

  LOG(DEBUG) << "before finishGhostManipulation, test: " << *test;
  test->finishGhostManipulation();

  LOG(DEBUG) << "after finishGhostManipulation, test: " << *test;*/

  this->data_.rhs()->setRepresentationGlobal();
  this->data_.rhs()->zeroEntries();
  this->data_.rhs()->startGhostManipulation();
  //this->data_.rhs()->zeroGhostBuffer();

  LOG(DEBUG) << "after startGhostManipulation, rhs: " << *this->data_.rhs();

  typedef typename NeumannBoundaryConditionsBase<FunctionSpaceType,QuadratureType,nComponents>::ElementWithFaces ElementWithFaces;
  typedef FunctionSpace::FunctionSpace<typename FunctionSpaceType::SurfaceMesh,
                                       typename FunctionSpaceType::BasisFunction::BasisFunctionUsingOnlyNodalValues>
                                       FunctionSpaceSurface;

  std::shared_ptr<FunctionSpaceType> functionSpace = this->data_.functionSpace();

  const int D = FunctionSpaceType::dim();  // = 2 or 3
  const int nDofsPerElement = FunctionSpaceType::nDofsPerElement();

  // use gauss quadrature with 3 points for surface BCs, theoretically, QuadratureType could be used, but then there would be a partial specialization necessary for Quadrature::None and regular meshes
  typedef Quadrature::Gauss<3> QuadratureTypeSurface;

  // define shortcuts for quadrature
  typedef Quadrature::TensorProduct<D-1,QuadratureTypeSurface> QuadratureSurface;

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

    // check if element no is valid
    if (elementNoLocal < 0 || elementNoLocal > functionSpace->nElementsLocal())
    {
      LOG(ERROR) << "Element local no. " << elementNoLocal << " for which Neumann BC is specified, is invalid (number of local elements: " << functionSpace->nElementsLocal() << ")";
      continue;
    }

    std::array<Vec3,FunctionSpaceType::nDofsPerElement()> geometryVolume;
    std::array<Vec3,FunctionSpaceSurface::nDofsPerElement()> geometrySurface;  // geometry surface is not Hermite, if volume geometry uses Hermite, surface will be linear Lagrange
    functionSpace->getElementGeometry(elementNoLocal, geometryVolume);
    functionSpace->extractSurfaceGeometry(geometryVolume, elementIter->face, geometrySurface);

    VLOG(1) << "element no " << elementNoLocal << ", dofVectors: " << elementIter->dofVectors << ", geometryVolume: "
      << geometryVolume << ", face " << Mesh::getString(elementIter->face) << ", geometrySurface: " << geometrySurface;

    // show node positions
      /*
    if (VLOG_IS_ON(1))
    {
      std::array<dof_no_t,FunctionSpaceType::nNodesPerElement()> elementNodeNos = functionSpace->getElementNodeNos(elementNoLocal);
      for (int elementalNodeIndex = 0; elementalNodeIndex < FunctionSpaceType::nNodesPerElement(); elementalNodeIndex++)
      {
        dof_no_t nodeNoLocal = elementNodeNos[elementalNodeIndex];
        //LOG(DEBUG) << ", global coordinates of element: " << functionSpace->meshPartition()->getCoordinatesGlobal(nodeNoLocal);
      }
    }*/

    // loop over integration points (e.g. gauss points)
    for (unsigned int samplingPointIndex = 0; samplingPointIndex < samplingPointsSurface.size(); samplingPointIndex++)
    {
      // get parameter values of current sampling point
      std::array<double,D-1> xiSurface = samplingPointsSurface[samplingPointIndex];
      VecD<D> xi = Mesh::getXiOnFace(elementIter->face, xiSurface);

      // compute the 3xD jacobian of the parameter space to world space mapping
      std::array<Vec3,D-1> jacobian = FunctionSpaceSurface::computeJacobian(geometrySurface, xiSurface);
      double integrationFactor = MathUtility::computeIntegrationFactor<D-1>(jacobian);

      VLOG(1) << "   jacobian: " << jacobian;

      // set all entries to 0
      evaluationsArraySurface[samplingPointIndex] = {0.0};

      // integral to solve:
      // int_e phi_i(xi) * f(xi) dxi
      // first, solve f(xi) = sum_i psi_i(xi) * f_i

      // compute the value at xi by summing contributions from elemental dofs, weighted with the (D-1)-dimensional Lagrange ansatz functions on the face at xiSurface
      // this corresponds to f(xi) = sum_i psi_i(xi) * f_i, where psi are the Lagrange ansatz functions on the face/surface
      VecD<nComponents> boundaryConditionValueAtXi({0.0});
      for (typename std::vector<std::pair<dof_no_t, VecD<nComponents>>>::const_iterator dofVectorsIter = elementIter->dofVectors.begin();
           dofVectorsIter != elementIter->dofVectors.end();
           dofVectorsIter++)
      {
        int dofIndex = dofVectorsIter->first;
        VecD<nComponents> fluxValue = dofVectorsIter->second;   // this is the prescribed value, either a scalar of the flux or the traction vector

        boundaryConditionValueAtXi += fluxValue * FunctionSpaceSurface::phi(dofIndex, xiSurface);
      }

      // now add contribution of phi_i(xi) * f(xi)

      // loop over all dofs of element with given boundary condition value
      for (typename std::vector<dof_no_t>::const_iterator surfaceDofIter = elementIter->surfaceDofs.begin();
           surfaceDofIter != elementIter->surfaceDofs.end();
           surfaceDofIter++)
      {
        int surfaceDofIndex = *surfaceDofIter;

        VecD<nComponents> dofIntegrand = boundaryConditionValueAtXi * functionSpace->phi(surfaceDofIndex, xi) * integrationFactor;

        VLOG(2) << "  surfaceDofIndex " << surfaceDofIndex << ", xi=" << xi << ", BC value: " << boundaryConditionValueAtXi
          << " phi = " << functionSpace->phi(surfaceDofIndex, xi) << ", integrationFactor: " << integrationFactor << ", dofIntegrand: " << dofIntegrand;

        // store integrand in evaluations array
        for (int i = 0; i < nComponents; i++)
        {
          evaluationsArraySurface[samplingPointIndex][surfaceDofIndex*nComponents + i] = dofIntegrand[i];
        }
      }  // surfaceDofIndex

    }  // function evaluations

    // integrate all values for result vector entries at once
    EvaluationsType integratedValues = QuadratureSurface::computeIntegral(evaluationsArraySurface);

    VLOG(1) << "evaluationsArraySurface: " << evaluationsArraySurface << ", integratedValues: " << integratedValues;

    // get indices of element-local dofs
    std::array<dof_no_t,nDofsPerElement> dofNosLocal = functionSpace->getElementDofNosLocal(elementNoLocal);

    // add entries in result vector
    // loop over indices of unknows (dofIndex,dofComponent)
    for (int dofIndex = 0; dofIndex < nDofsPerElement; dofIndex++)
    {
      std::array<double,nComponents> dofValues;
      std::copy(integratedValues.begin() + dofIndex*nComponents, integratedValues.begin() + (dofIndex+1)*nComponents, dofValues.begin());

      VLOG(2) << "set values " << dofValues << " at dof " << dofNosLocal[dofIndex];

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

        VLOG(2) << "  dof " << dofIndex << ", component " << dofComponent << " integrated value: " << integratedValue;
      }  // dofComponent
#endif
    }  // dofIndex
  } // elementGlobalNo


  LOG(DEBUG) << "before finishGhostManipulation, rhs: " << *this->data_.rhs();
  this->data_.rhs()->finishGhostManipulation();

  VLOG(1) << "after initializeRhs, rhs: " << *this->data_.rhs();
  LOG(DEBUG) << "after initializeRhs, rhs: " << *this->data_.rhs();
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

    VLOG(1) << "element no " << elementNoLocal << ", dofVectors: " << elementIter->dofVectors[0];

    // check if element no is valid
    if (elementNoLocal < 0 || elementNoLocal > functionSpace->nElementsLocal())
    {
      LOG(ERROR) << "Element local no. " << elementNoLocal << " for which Neumann BC is specified, is invalid (number of local elements: " << functionSpace->nElementsLocal() << ")";
      continue;
    }

    assert(!elementIter->dofVectors.empty());

    // get indices of element-local dofs
    std::array<dof_no_t,nDofsPerElement> dofNosLocal = functionSpace->getElementDofNosLocal(elementNoLocal);

    LOG(DEBUG) << "set value " << elementIter->dofVectors[0].second << " at dof " << dofNosLocal[elementIter->dofVectors[0].first];

    // add entries in result vector
    this->data_.rhs()->setValue(dofNosLocal[elementIter->dofVectors[0].first], elementIter->dofVectors[0].second, ADD_VALUES);

  } // elementGlobalNo

  this->data_.rhs()->finishGhostManipulation();
}

// 2D,3D, nComponents > 1
template<typename FunctionSpaceType,typename QuadratureType,int nComponents,typename DummyForTraits>
typename NeumannBoundaryConditions<FunctionSpaceType,QuadratureType,nComponents,DummyForTraits>::ElementWithFaces
NeumannBoundaryConditions<FunctionSpaceType,QuadratureType,nComponents,DummyForTraits>::
parseElementWithFaces(PythonConfig specificSettings, std::shared_ptr<FunctionSpaceType> functionSpace, element_no_t elementNoLocal)
{
  // extract one item of the Neumann BC specification, e.g. {"element": 1, "face": "0+", "dofVectors:", {0: [tmax,0,0], 1: [tmax,0,0], 2: [tmax,0,0], 3: [tmax,0,0]}}
  // nComponents == FunctionSpaceType::dim() for traction boundary condition or nComponents = 1 for flux BC

  ElementWithFaces result;

  // store local element no
  if (elementNoLocal >= 0)
  {
    // if a value for elementNoLocal was given as parameter, use this value
    result.elementNoLocal = elementNoLocal;
  }
  else
  {
    // parse from config
    int elementNoLocal = specificSettings.getOptionInt("element", 0);

    // negative element nos count from the end
    if (elementNoLocal < 0)
      elementNoLocal += functionSpace->nElementsLocal();

    if (elementNoLocal < 0 || elementNoLocal >= functionSpace->nElementsLocal())
    {
      LOG(FATAL) << "In Neumann boundary conditions, local element no. " << elementNoLocal << " is invalid (local number of elements: " << functionSpace->nElementsLocal() << ")";
    }
    result.elementNoLocal = elementNoLocal;
  }

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
        constantVector = MathUtility::transformToD<D,3>(functionSpace->getNormal(result.face, result.elementNoLocal, xi) * constantValue);
      }
    }
    else if (specificSettings.hasKey("constantVector"))
    {
      constantVector = specificSettings.getOptionArray<double,nComponents>("constantVector", 0.0);
    }

    // get dofs on surface where BC is specified
    typedef FunctionSpace::FunctionSpace<typename FunctionSpaceType::SurfaceMesh,
                                         typename FunctionSpaceType::BasisFunction::BasisFunctionUsingOnlyNodalValues>
                                         FunctionSpaceSurface;

    const int nDofs = FunctionSpaceSurface::nDofsPerElement();
    for (int i = 0; i < nDofs; i++)
    {
      result.dofVectors.push_back(std::pair<dof_no_t, VecD<nComponents>>(i, constantVector));
    }

    // get dofs indices within the numbering of the volume element that correspond to the selected face
    const int D = FunctionSpaceType::dim();
    const int nDofsPerNode = FunctionSpaceType::nDofsPerNode();
    const int nVolumeDofsBorder = FunctionSpace::FunctionSpaceBaseDim<D-1,typename FunctionSpaceType::BasisFunction>::nNodesPerElement() * nDofsPerNode;
    std::array<dof_no_t,nVolumeDofsBorder> dofIndices;
    FunctionSpaceType::getFaceDofs(result.face, dofIndices);

    LOG(DEBUG) << "nVolumeDofsBorder on " << D-1 << "D face: " << nVolumeDofsBorder << ": " << dofIndices;

    for (int i = 0; i < nVolumeDofsBorder; i++)
    {
      result.surfaceDofs.push_back(dofIndices[i]);
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
parseElementWithFaces(PythonConfig specificSettings, std::shared_ptr<FunctionSpaceType> functionSpace, element_no_t elementNoLocal)
{
  // extract one item of the Neumann BC specification, e.g. {"element": 1, "face": "0+", "dofVectors:", {0: [tmax,0,0], 1: [tmax,0,0], 2: [tmax,0,0], 3: [tmax,0,0]}}
  // nComponents = 1 for flux BC

  //const int D = FunctionSpaceType::dim();
  //LOG(TRACE) << "NeumannBoundaryConditions::parseElementWithFaces, D=" << D << ", nComponents=1";

  ElementWithFaces result;

  // store local element no
  if (elementNoLocal >= 0)
  {
    // if a value for elementNoLocal was given as parameter, use this value
    result.elementNoLocal = elementNoLocal;
  }
  else
  {
    // parse from config
    int elementNoLocal = specificSettings.getOptionInt("element", 0);

    // negative element nos count from the end
    if (elementNoLocal < 0)
      elementNoLocal += functionSpace->nElementsLocal();

    if (elementNoLocal < 0 || elementNoLocal >= functionSpace->nElementsLocal())
    {
      LOG(FATAL) << "In Neumann boundary conditions, local element no. " << elementNoLocal << " is invalid (local number of elements: " << functionSpace->nElementsLocal() << ")";
    }
    result.elementNoLocal = elementNoLocal;
  }

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
/*
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
*/
      // for normal flux BC
      constantVector[0] = constantValue;
    }
    else if (specificSettings.hasKey("constantVector"))
    {
      constantVector = specificSettings.getOptionArray<double,1>("constantVector", 0.0);
    }

    // get dofs on surface where BC is specified
    typedef FunctionSpace::FunctionSpace<typename FunctionSpaceType::SurfaceMesh,
                                         typename FunctionSpaceType::BasisFunction::BasisFunctionUsingOnlyNodalValues>
                                         FunctionSpaceSurface;

    // number of dofs for Neumann BC, this is equal to the number of nodes (also for Hermite), because only nodal dofs are considered
    const int nDofsBorder = FunctionSpaceSurface::nDofsPerElement();
    for (int i = 0; i < nDofsBorder; i++)
    {
      result.dofVectors.push_back(std::pair<dof_no_t, VecD<1>>(i, constantVector));
    }

    // get dofs indices within the numbering of the volume element that correspond to the selected face
    const int D = FunctionSpaceType::dim();
    const int nDofsPerNode = FunctionSpaceType::nDofsPerNode();
    const int nVolumeDofsBorder = FunctionSpace::FunctionSpaceBaseDim<D-1,typename FunctionSpaceType::BasisFunction>::nNodesPerElement() * nDofsPerNode;
    std::array<dof_no_t,nVolumeDofsBorder> dofIndices;
    FunctionSpaceType::getFaceDofs(result.face, dofIndices);

    LOG(DEBUG) << "nVolumeDofsBorder on " << D-1 << "D face: " << nVolumeDofsBorder << ": " << dofIndices;

    for (int i = 0; i < nVolumeDofsBorder; i++)
    {
      result.surfaceDofs.push_back(dofIndices[i]);
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
parseElementWithFaces(PythonConfig specificSettings, std::shared_ptr<FunctionSpaceType> functionSpace, element_no_t elementNoLocal)
{
  // extract one item of the Neumann BC specification, e.g. {"element": 1, "face": "0+", "dofVectors:", {0: [tmax,0,0], 1: [tmax,0,0], 2: [tmax,0,0], 3: [tmax,0,0]}}
  // nComponents == FunctionSpaceType::dim() for traction boundary condition or nComponents = 1 for flux BC

  ElementWithFaces result;

  // store local element no
  if (elementNoLocal >= 0)
  {
    // if a value for elementNoLocal was given as parameter, use this value
    result.elementNoLocal = elementNoLocal;
  }
  else
  {
    // parse from config
    int elementNoLocal = specificSettings.getOptionInt("element", 0);

    // negative element nos count from the end
    if (elementNoLocal < 0)
      elementNoLocal += functionSpace->nElementsLocal();

    if (elementNoLocal < 0 || elementNoLocal >= functionSpace->nElementsLocal())
    {
      LOG(FATAL) << "In Neumann boundary conditions, local element no. " << elementNoLocal << " is invalid (local number of elements: " << functionSpace->nElementsLocal() << ")";
    }
    result.elementNoLocal = elementNoLocal;
  }

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
      constantVector[0] = constantValue;
    }
    else if (specificSettings.hasKey("constantVector"))
    {
      constantVector = specificSettings.getOptionArray<double,1>("constantVector", 0.0);
    }

    // get dofs indices within the numbering of the volume element that correspond to the selected face
    const int nDofsPerNode = FunctionSpaceType::nDofsPerNode();
    std::array<dof_no_t,nDofsPerNode> dofIndices;
    FunctionSpaceType::getFaceDofs(result.face, dofIndices);
    // for Hermite we get 2 dofs in dofIndices, only use the first one

    VLOG(1) << "nVolumeDofsBorder on 0D face " << Mesh::getString(result.face) << ": " << dofIndices;

    result.dofVectors.push_back(std::pair<dof_no_t, VecD<1>>(dofIndices[0], constantVector));
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
