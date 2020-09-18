#include "field_variable/08_field_variable_gradient.h"

#include "control/dihu_context.h"

namespace FieldVariable
{

//#define USE_APPROXIMATE_GRADIENT    // if the method of approximating the gradient by difference quotients between neighouring nodes should be used (this is not so good in general)
// the alternative is to compute the gradient from the ansatz functions (this is more accurate but bad for badly conditioned elements)

const int CONDITION_TOLERANCE = 25;    // condition number value, if the condition number is higher, the dof will get the approximated gradient

// structured mesh
//! compute the gradient field
template<typename FunctionSpaceType>
void FieldVariableGradient<FunctionSpaceType,1,::Mesh::isStructured<typename FunctionSpaceType::Mesh>>::
computeGradientField(std::shared_ptr<FieldVariable<FunctionSpaceType, FunctionSpaceType::dim()>> gradientField,
                     std::shared_ptr<FieldVariable<FunctionSpaceType,1>> jacobianConditionNumberField)
{
  LOG(DEBUG) << "computeGradientField (structured), functionSpaceType: " << StringUtility::demangle(typeid(FunctionSpaceType).name());

  this->values_->setRepresentationGlobal();
  this->values_->startGhostManipulation();

  // initialize gradient field variable to 0
  gradientField->zeroEntries();

  gradientField->setRepresentationGlobal();
  gradientField->startGhostManipulation();
  gradientField->zeroGhostBuffer();

  if (jacobianConditionNumberField)
  {
    jacobianConditionNumberField->zeroEntries();

    jacobianConditionNumberField->setRepresentationGlobal();
    jacobianConditionNumberField->startGhostManipulation();
    jacobianConditionNumberField->zeroGhostBuffer();
  }

  std::shared_ptr<FieldVariable<FunctionSpaceType, FunctionSpaceType::dim()>> approximatedGradientField
    = std::make_shared<FieldVariable<FunctionSpaceType, FunctionSpaceType::dim()>>(*gradientField, "approximatedGradient");
  approximatedGradientField->zeroEntries();
  approximatedGradientField->setRepresentationGlobal();
  approximatedGradientField->startGhostManipulation();
  approximatedGradientField->zeroGhostBuffer();

  // define constants
  const int nDofsPerElement = FunctionSpaceType::nDofsPerElement();
  const int D = FunctionSpaceType::dim();

  const dof_no_t nDofsWithGhosts = this->functionSpace_->nDofsLocalWithGhosts();
  std::vector<int> nAdjacentElements(nDofsWithGhosts, 0);   //< the number of elements that are adjacent to the node

  // count number evaluations for every dof
  // loop over elements
  for (element_no_t elementNoLocal = 0; elementNoLocal < this->functionSpace_->nElementsLocal(); elementNoLocal++)
  {
    // get local dof nos of this element
    std::array<dof_no_t,nDofsPerElement> elementDofs = this->functionSpace_->getElementDofNosLocal(elementNoLocal);

    // loop over dofs in element, where to compute the gradient
    for (int dofIndex = 0; dofIndex < nDofsPerElement; dofIndex++)
    {
      dof_no_t dofNo = elementDofs[dofIndex];

      assert(dofNo < nDofsWithGhosts);

      // increase counter of number of summands for that dof
      nAdjacentElements[dofNo]++;
    }
  }

  // compute gradient value divided by number of evaluations
  // loop over elements
  for (element_no_t elementNoLocal = 0; elementNoLocal < this->functionSpace_->nElementsLocal(); elementNoLocal++)
  {
    // get local dof nos of this element
    std::array<dof_no_t,nDofsPerElement> elementDofs = this->functionSpace_->getElementDofNosLocal(elementNoLocal);

    // compute gradient at every dof, as continuous to current element (gradients have discontinuities between elements at dofs)
    std::array<double,nDofsPerElement> solutionValues;
    this->getElementValues(elementNoLocal, solutionValues);

    // get geometry field (which are the node positions for Lagrange basis and node positions and derivatives for Hermite)
    std::array<Vec3,nDofsPerElement> geometryValues;
    this->functionSpace_->getElementGeometry(elementNoLocal, geometryValues);
    double approximateMeshWidth = MathUtility::computeApproximateMeshWidth<double,nDofsPerElement>(geometryValues);

    const int nDofsPerNode = this->functionSpace_->nDofsPerNode();
    std::array<double,D> xi;

    //LOG(DEBUG) << "mesh " << this->functionSpace_->meshName() << " element " << elementNoLocal << ", values: " << solutionValues;

    // loop over dofs in element, where to compute the gradient
    for (int dofIndex = 0; dofIndex < nDofsPerElement; dofIndex++)
    {

      dof_no_t dofNo = elementDofs[dofIndex];

      // compute correct gradient from ansatz functions
      // set xi to dofIndex
      for (int i = 0; i < D; i++)
      {
        if (i == 0)
          xi[i] = double(dofIndex % 2);
        else if (i == 1)
          xi[i] = double((dofIndex % 4) / 2);
        else if (i == 2)
          xi[i] = double(dofIndex / 4);
      }

      //VLOG(2) << "element " << elementNoLocal << " dofIndex " << dofIndex << ", xi " << xi << " g:" << geometryValues;

      // compute the 3xD jacobian of the parameter space to world space mapping
      Tensor2<D> jacobianParameterSpace = MathUtility::transformToDxD<D,D>(FunctionSpaceType::computeJacobian(geometryValues, xi));
      double jacobianDeterminant;
      Tensor2<D> inverseJacobianParameterSpace = MathUtility::template computeInverse<double>(jacobianParameterSpace, approximateMeshWidth, jacobianDeterminant);

      // estimate condition value of jacobian
      if (jacobianConditionNumberField != nullptr)
      {
        double conditionNumber = MathUtility::estimateConditionNumber(jacobianParameterSpace, inverseJacobianParameterSpace);

        conditionNumber /= nAdjacentElements[dofNo];
        jacobianConditionNumberField->setValue(dofNo, conditionNumber, ADD_VALUES);
      }
      //VLOG(1) << "conditionNumber: " << conditionNumber;

      // get gradient at dof
      std::array<double,D> gradPhiWorldSpace = this->functionSpace_->interpolateGradientInElement(solutionValues, inverseJacobianParameterSpace, xi);

      //LOG(DEBUG) << "   dof " << dofIndex << ", dofNo " << dofNo << ", nAdjacentElements: " << nAdjacentElements[dofNo] << ", gradPhiWorldSpace: " << gradPhiWorldSpace << ", inverseJacobianParameterSpace: " << inverseJacobianParameterSpace;

      // scale value
      gradPhiWorldSpace /= nAdjacentElements[dofNo];

      // add value to gradient field variable
      if (VLOG_IS_ON(2))
      {
        int rankNo = DihuContext::ownRankNoCommWorld();

        if ((rankNo == 0 && dofNo == 150) || (rankNo == 1 && dofNo == 0))
        {
          LOG(DEBUG) << "dofNo " << dofNo << " gradPhiWorldSpace: " << gradPhiWorldSpace << ", nAdjacentElements[dofNo]: " << nAdjacentElements[dofNo]
            << ",geometryValues: " << geometryValues << ", for this dof: " << geometryValues[dofIndex];
          LOG(DEBUG) << "solutionValues: " << solutionValues[dofIndex] << ", all: " << solutionValues;
        }
      }

      gradientField->setValue(dofNo, gradPhiWorldSpace, ADD_VALUES);

      // -----------------------------------
      // alternative computation of gradient
      // for every adjacent face to the node, compute the difference quotient times the normalized vector of that direction
      int nodeIndex = dofIndex / nDofsPerNode;

      // for nodes
      if (dofIndex % nDofsPerNode == 0)
      {
        std::array<double,D> gradPhiWorldSpace2({0.0});
        Vec3 nodePosition = geometryValues[dofIndex];
        double solutionValue = solutionValues[dofIndex];

        // loop over all directions, all faces that are adjacent to this node in this element (3 for 3D)
        for (int faceNo = 0; faceNo <= (int)Mesh::face_t::face2Plus; faceNo++)
        {
          int neighbourNodeIndex = this->functionSpace_->getNeighbourNodeIndex(nodeIndex, (Mesh::face_t)faceNo);

          if (neighbourNodeIndex != -1)   // if the neighbouring direction is valid inside the element
          {
            int neighbourDofIndex = neighbourNodeIndex * nDofsPerNode;
            Vec3 neighbourNodePosition = geometryValues[neighbourDofIndex];
            double neighbourSolutionValue = solutionValues[neighbourDofIndex];

            Vec3 neighbourDirection = neighbourNodePosition - nodePosition;
            double faceLengthSquared = MathUtility::normSquared<3>(neighbourDirection);

            gradPhiWorldSpace2 += neighbourDirection * ((neighbourSolutionValue - solutionValue) / faceLengthSquared);
          }
        }

        gradPhiWorldSpace2 /= nAdjacentElements[dofNo] * D;  // D directions per element

        VecD<D> gradPhiWorldSpaceD;
        for (int i = 0; i < D; i++)
        {
          gradPhiWorldSpaceD[i] = gradPhiWorldSpace2[i];
        }

        approximatedGradientField->setValue(dofNo, gradPhiWorldSpaceD, ADD_VALUES);
      }
    }  // dofIndex
  }  // elementNoLocal

  gradientField->finishGhostManipulation();
  approximatedGradientField->finishGhostManipulation();

  LOG(DEBUG) << "gradientField: " << *gradientField;
  LOG(DEBUG) << "approximatedGradientField: " << *approximatedGradientField;

  if (jacobianConditionNumberField)
  {
    jacobianConditionNumberField->finishGhostManipulation();
  }

  //LOG(DEBUG) << "gradientField: " << *gradientField;

  // fix bad values in the gradient field
  if (jacobianConditionNumberField)
  {
    // loop over local nodes
    for (node_no_t dofNoLocal = 0; dofNoLocal < this->functionSpace_->nDofsLocalWithoutGhosts(); dofNoLocal++)
    {
      node_no_t nodeNoLocal = dofNoLocal / this->functionSpace_->nDofsPerNode();

      // get condition number
      double conditionNumber = jacobianConditionNumberField->getValue(dofNoLocal);

      if (conditionNumber > CONDITION_TOLERANCE)
      {
        LOG(DEBUG) << "dof " << dofNoLocal << " has condition number " << conditionNumber << ", use approximated gradient value";

        Vec3 gradientSum;
        int nSummands = 0;

        for (int face = (int)Mesh::face_t::face0Minus; face <= (int)Mesh::face_t::face2Plus; face++)
        {
          int neighbourNodeNo = this->functionSpace_->getNeighbourNodeNoLocal(nodeNoLocal, (Mesh::face_t)face);
          if (neighbourNodeNo != -1)
          {
            int neighbourDofNo = neighbourNodeNo*this->functionSpace_->nDofsPerNode();
            double neighbourJacobianConditionNumber = jacobianConditionNumberField->getValue(neighbourDofNo);

            if (fabs(neighbourJacobianConditionNumber) <= CONDITION_TOLERANCE)
            {
              gradientSum += gradientField->getValue(neighbourDofNo);
              nSummands++;
              LOG(DEBUG) << "  neighbour " << Mesh::getString((Mesh::face_t)face) << ", node " << neighbourNodeNo
                << ", dof " << neighbourDofNo << "  has condition number "
                << neighbourJacobianConditionNumber << ", value " << gradientField->getValue(neighbourDofNo);
            }
            else
            {
              LOG(DEBUG) << "  neighbour " << Mesh::getString((Mesh::face_t)face) << ", node " << neighbourNodeNo
                << ", dof " << neighbourDofNo << "  has condition number "
                << neighbourJacobianConditionNumber << ", do not use gradient";
            }
          }
        }

        gradientSum /= nSummands;
        gradientField->setValue(dofNoLocal, gradientSum, INSERT_VALUES);
        LOG(DEBUG) << "node " << nodeNoLocal << " has " << nSummands << " summands, result: " << gradientSum;
      }
    }
  }
}

// unstructured mesh
//! compute the gradient field
template<typename FunctionSpaceType>
void FieldVariableGradient<FunctionSpaceType,1,::Mesh::UnstructuredDeformableOfDimension<FunctionSpaceType::dim()>>::
computeGradientField(std::shared_ptr<FieldVariable<FunctionSpaceType, FunctionSpaceType::dim()>> gradientField,
                     std::shared_ptr<FieldVariable<FunctionSpaceType,1>> jacobianConditionNumberField)
{
  LOG(DEBUG) << "computeGradientField (unstructured), functionSpaceType: " << StringUtility::demangle(typeid(FunctionSpaceType).name());

  this->values_->setRepresentationGlobal();
  this->values_->startGhostManipulation();

  // initialize gradient field variable to 0
  gradientField->zeroEntries();

  gradientField->setRepresentationGlobal();
  gradientField->startGhostManipulation();
  gradientField->zeroGhostBuffer();

  // define constants
  const int nDofsPerElement = FunctionSpaceType::nDofsPerElement();
  const int D = FunctionSpaceType::dim();

  const dof_no_t nDofsWithGhosts = this->functionSpace_->nDofsLocalWithGhosts();
  std::vector<int> nAdjacentElements(nDofsWithGhosts, 0);   //< the number of elements that are adjacent to the node

  // count number evaluations for every dof
  // loop over elements
  for (element_no_t elementNoLocal = 0; elementNoLocal < this->functionSpace_->nElementsLocal(); elementNoLocal++)
  {
    // get local dof nos of this element
    std::array<dof_no_t,nDofsPerElement> elementDofs = this->functionSpace_->getElementDofNosLocal(elementNoLocal);

    // loop over dofs in element, where to compute the gradient
    for (int dofIndex = 0; dofIndex < nDofsPerElement; dofIndex++)
    {
      dof_no_t dofNo = elementDofs[dofIndex];

      assert(dofNo < nDofsWithGhosts);

      // increase counter of number of summands for that dof
      nAdjacentElements[dofNo]++;
    }
  }

  // compute gradient value divided by number of evaluations
  // loop over elements
  for (element_no_t elementNoLocal = 0; elementNoLocal < this->functionSpace_->nElementsLocal(); elementNoLocal++)
  {
    // get local dof nos of this element
    std::array<dof_no_t,nDofsPerElement> elementDofs = this->functionSpace_->getElementDofNosLocal(elementNoLocal);

    // compute gradient at every dof, as continuous to current element (gradients have discontinuities between elements at dofs)
    std::array<double,nDofsPerElement> solutionValues;
    this->getElementValues(elementNoLocal, solutionValues);

    // get geometry field (which are the node positions for Lagrange basis and node positions and derivatives for Hermite)
    std::array<Vec3,nDofsPerElement> geometryValues;
    this->functionSpace_->getElementGeometry(elementNoLocal, geometryValues);
    double approximateMeshWidth = MathUtility::computeApproximateMeshWidth<double,nDofsPerElement>(geometryValues);

    std::array<double,D> xi;

    // loop over dofs in element, where to compute the gradient
    for (int dofIndex = 0; dofIndex < nDofsPerElement; dofIndex++)
    {

      dof_no_t dofNo = elementDofs[dofIndex];

      // compute correct gradient from ansatz functions
      // set xi to dofIndex
      for (int i = 0; i < D; i++)
      {
        if (i == 0)
          xi[i] = double(dofIndex % 2);
        else if (i == 1)
          xi[i] = double((dofIndex % 4) / 2);
        else if (i == 2)
          xi[i] = double(dofIndex / 4);
      }

      //VLOG(2) << "element " << elementNoLocal << " dofIndex " << dofIndex << ", xi " << xi << " g:" << geometryValues;

      // compute the 3xD jacobian of the parameter space to world space mapping
      Tensor2<D> jacobianParameterSpace = MathUtility::transformToDxD<D,D>(FunctionSpaceType::computeJacobian(geometryValues, xi));
      double jacobianDeterminant;
      Tensor2<D> inverseJacobianParameterSpace = MathUtility::template computeInverse<double>(jacobianParameterSpace, approximateMeshWidth, jacobianDeterminant);

      // estimate condition value of jacobian
      double conditionNumber = MathUtility::estimateConditionNumber(jacobianParameterSpace, inverseJacobianParameterSpace);
      if (jacobianConditionNumberField != nullptr)
      {
        conditionNumber /= nAdjacentElements[dofNo];
        jacobianConditionNumberField->setValue(dofNo, conditionNumber, ADD_VALUES);
      }
      //VLOG(1) << "conditionNumber: " << conditionNumber;

      // get gradient at dof
      std::array<double,D> gradPhiWorldSpace = this->functionSpace_->interpolateGradientInElement(solutionValues, inverseJacobianParameterSpace, xi);

      // scale value
      gradPhiWorldSpace /= nAdjacentElements[dofNo];

      // add value to gradient field variable
      if (VLOG_IS_ON(2))
      {
        int rankNo = DihuContext::ownRankNoCommWorld();

          if ((rankNo == 0 && dofNo == 150) || (rankNo == 1 && dofNo == 0))
        {
          LOG(DEBUG) << "dofNo " << dofNo << " gradPhiWorldSpace: " << gradPhiWorldSpace << ", nAdjacentElements[dofNo]: " << nAdjacentElements[dofNo]
            << ",geometryValues: " << geometryValues << ", for this dof: " << geometryValues[dofIndex];
          LOG(DEBUG) << "solutionValues: " << solutionValues[dofIndex] << ", all: " << solutionValues;
        }
      }

      gradientField->setValue(dofNo, gradPhiWorldSpace, ADD_VALUES);
    }  // dofIndex
  }  // elementNoLocal

  gradientField->finishGhostManipulation();
}

} // namespace
