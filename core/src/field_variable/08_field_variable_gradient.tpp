#include "field_variable/08_field_variable_vector.h"

#include "control/dihu_context.h"

namespace FieldVariable
{

//#define USE_APPROXIMATE_GRADIENT    // if the method of approximating the gradient by difference quotients between neighouring nodes should be used (this is not so good in general)
// the alternative is to compute the gradient from the ansatz functions (this is more accurate but bad for badly conditioned elements)

//! compute the gradient field
template<typename FunctionSpaceType>
void FieldVariableVector<FunctionSpaceType,1>::
computeGradientField(std::shared_ptr<FieldVariable<FunctionSpaceType, FunctionSpaceType::dim()>> gradientField)
{
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
  std::vector<int> nAdjacentElements(nDofsWithGhosts, 0);   ///< the number of elements that are adjacent to the node

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

#ifdef USE_APPROXIMATE_GRADIENT
    const int nDofsPerNode = this->functionSpace_->nDofsPerNode();
#else
    std::array<double,D> xi;
#endif

    // loop over dofs in element, where to compute the gradient
    for (int dofIndex = 0; dofIndex < nDofsPerElement; dofIndex++)
    {

      dof_no_t dofNo = elementDofs[dofIndex];

#ifndef USE_APPROXIMATE_GRADIENT
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
      Tensor2<D> inverseJacobianParameterSpace = MathUtility::computeInverse<D>(jacobianParameterSpace, jacobianDeterminant);


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

      //VecD<D> test({1.0 + 0.1*(1+rankNo)});
      //LOG(DEBUG) << "test: " << test;
      //gradientField->setValue(dofNo, test, ADD_VALUES);

      /*
      if (fabs(gradPhiWorldSpace[0])+fabs(gradPhiWorldSpace[1]) > fabs(gradPhiWorldSpace[2]) || dofIndex == 123)
      {
        VLOG(2) << "element " << elementNoLocal << " dofIndex " << dofIndex << ", xi " << xi << " geometry:" << geometryValues
          << ", solution: " << solutionValues;
        VLOG(2) << "   jacobianParameterSpace: " << jacobianParameterSpace << ", inverseJacobianParameterSpace: " << inverseJacobianParameterSpace << ", jacobianDeterminant: " << jacobianDeterminant;
        VLOG(2) << "   local dof " << dofNo << ", add value " << MathUtility::normalized<3>(gradPhiWorldSpace) << "/" << nAdjacentElements[dofNo];
        VLOG(2) << "extreml: " << (fabs(gradPhiWorldSpace[0])+fabs(gradPhiWorldSpace[1])) / fabs(gradPhiWorldSpace[2]);

        for (int i = 0; i < 10; i++)
        {
          double f = 0.1;
          xi += Vec3({(rand()%101/50.-1.)*f,(rand()%101/50.-1.)*f,(rand()%101/50.-1.)*f});
          Vec3 a = this->functionSpace_->interpolateGradientInElement(solutionValues, inverseJacobianParameterSpace, xi);
          VLOG(2) << "  xi: " << xi << ", grad: " << MathUtility::normalized<3>(a);
        }
        VLOG(2) << " ";
        VLOG(2) << " ";
      }*/
      gradientField->setValue(dofNo, gradPhiWorldSpace, ADD_VALUES);
#endif


#ifdef USE_APPROXIMATE_GRADIENT
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

        gradientField->setValue(dofNo, gradPhiWorldSpaceD, ADD_VALUES);
      }
#endif
    }  // dofIndex
  }  // elementNoLocal

  gradientField->finishGhostManipulation();
  //gradientField->setRepresentationLocal();

  //LOG(DEBUG) << "gradientField: " << *gradientField;
}

} // namespace
