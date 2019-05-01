#include "field_variable/08_field_variable_vector.h"

#include "control/dihu_context.h"

namespace FieldVariable
{

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
  std::vector<int> nSummands(nDofsWithGhosts, 0);   ///< the number of elements that are adjacent to the node

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
      nSummands[dofNo]++;
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

    std::array<double,D> xi;

    // loop over dofs in element, where to compute the gradient
    for (int dofIndex = 0; dofIndex < nDofsPerElement; dofIndex++)
    {
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

      dof_no_t dofNo = elementDofs[dofIndex];

      /*
      if (fabs(gradPhiWorldSpace[0])+fabs(gradPhiWorldSpace[1]) > fabs(gradPhiWorldSpace[2]) || dofIndex == 123)
      {
        VLOG(2) << "element " << elementNoLocal << " dofIndex " << dofIndex << ", xi " << xi << " geometry:" << geometryValues
          << ", solution: " << solutionValues;
        VLOG(2) << "   jacobianParameterSpace: " << jacobianParameterSpace << ", inverseJacobianParameterSpace: " << inverseJacobianParameterSpace << ", jacobianDeterminant: " << jacobianDeterminant;
        VLOG(2) << "   local dof " << dofNo << ", add value " << MathUtility::normalized<3>(gradPhiWorldSpace) << "/" << nSummands[dofNo];
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

      gradPhiWorldSpace /= nSummands[dofNo];

      // add value to gradient field variable
      if (VLOG_IS_ON(2))
      {
        int rankNo = DihuContext::ownRankNoCommWorld();

          if ((rankNo == 0 && dofNo == 150) || (rankNo == 1 && dofNo == 0))
        {
          LOG(DEBUG) << "dofNo " << dofNo << " gradPhiWorldSpace: " << gradPhiWorldSpace << ", nSummands[dofNo]: " << nSummands[dofNo]
            << ",geometryValues: " << geometryValues << ", for this dof: " << geometryValues[dofIndex];
          LOG(DEBUG) << "solutionValues: " << solutionValues[dofIndex] << ", all: " << solutionValues;
        }
      }

      gradientField->setValue(dofNo, gradPhiWorldSpace, ADD_VALUES);
      //VecD<D> test({1.0 + 0.1*(1+rankNo)});
      //LOG(DEBUG) << "test: " << test;
      //gradientField->setValue(dofNo, test, ADD_VALUES);

    }  // dofIndex
  }  // elementNoLocal

  gradientField->finishGhostManipulation();
  //gradientField->setRepresentationLocal();

  //LOG(DEBUG) << "gradientField: " << *gradientField;
}

} // namespace
