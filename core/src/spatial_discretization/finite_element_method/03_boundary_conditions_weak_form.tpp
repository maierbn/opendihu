#include "spatial_discretization/finite_element_method/03_boundary_conditions.h"

#include <Python.h>
#include <memory>
#include <vector>
#include <petscsys.h>

#include "quadrature/tensor_product.h"
#include "utility/vector_operators.h"

namespace SpatialDiscretization
{

template<typename BasisOnMeshType,typename QuadratureType,typename Term,typename Dummy>
void BoundaryConditions<BasisOnMeshType,QuadratureType,Term,Dummy>::
applyBoundaryConditions()
{
  applyBoundaryConditionsWeakForm();
}

template<typename BasisOnMeshType,typename QuadratureType,typename Term,typename Dummy>
void BoundaryConditions<BasisOnMeshType,QuadratureType,Term,Dummy>::
applyBoundaryConditionsWeakForm()
{
  std::shared_ptr<BasisOnMeshType> mesh = this->data_.mesh();

  // add weak form of Dirichlet BC to rhs
  const int nDofsPerNode = BasisOnMeshType::nDofsPerNode();

  // determine if the BC indices in the config are given for global or local dof nos
  bool inputMeshIsGlobal = PythonUtility::getOptionBool(this->specificSettings_, "inputMeshIsGlobal", true);
  int nDofs = 0;
  if (inputMeshIsGlobal)
  {
    nDofs = mesh->nDofsGlobal();
  }
  else
  {
    nDofs = mesh->nDofsLocalWithoutGhosts();
  }

  LOG(TRACE) << "applyBoundaryConditionsWeakForm";

  FieldVariable::FieldVariable<BasisOnMeshType,1> &rightHandSide = this->data_.rightHandSide();
  std::shared_ptr<PartitionedPetscMat<BasisOnMeshType>> stiffnessMatrix = this->data_.stiffnessMatrix();

  rightHandSide.startVectorManipulation();

  // Boundary conditions are specified for dof numbers, not nodes, such that for Hermite it is possible to prescribe derivatives.
  // However the ordering of the dofs is not known in the config for unstructured meshes. Therefore the ordering is special.
  // For every node there are as many values as dofs, in contiguous order.
  // Example for 2D Hermite, unstructured grid, 2x2 elements:
  //
  // node numbering:
  //  6_7_8
  // 3|_4_|5
  // 0|_1_|2
  //
  // dof numbering:
  //  6_7_8
  // 2|_3_|5
  // 0|_1_|4
  //
  // To specify du/dn = 0 an the left boundary you would set:
  // bc[0*2+1] = 0, bc[3*2+1] = 0, bc[6*2+1] = 0
  //
  // To specifiy u=0 on the bottom, you would set:
  // bc[0] = 0, bc[2] = 0, bc[4] = 0


  // parse all boundary conditions that are given in config
  std::vector<std::pair<int,double>> boundaryConditions;  /// (index, value) pairs

  std::pair<int, double> boundaryCondition = PythonUtility::getOptionDictBegin<int, double>(this->specificSettings_, "DirichletBoundaryCondition");
  for (; !PythonUtility::getOptionDictEnd(this->specificSettings_, "DirichletBoundaryCondition");
          PythonUtility::getOptionDictNext<int, double>(this->specificSettings_, "DirichletBoundaryCondition", boundaryCondition))
  {
    if (boundaryCondition.first < 0)
    {
      boundaryCondition.first += nDofs;
    }
    boundaryConditions.push_back(boundaryCondition);
  }
  if (boundaryConditions.empty())
  {
    LOG(DEBUG) << "no boundary conditions found";
    return;
  }

  // sort all parsed boundary conditions for their element no
  std::sort(boundaryConditions.begin(), boundaryConditions.end(), [](const std::pair<int,double> &item1, const std::pair<int,double> &item2){return item1.first < item2.first;});

  // determine elements with nodes that have prescribed Dirichlet boundary conditions, store them in the vector boundaryConditionElements, which is organized by local elements
  struct ElementWithNodes
  {
    element_no_t elementNoLocal;   ///< local element no
    std::vector<std::pair<int,double>> elementalDofIndex;   ///< the element-local dof index and the value of the boundary condition on this dof
  };
  std::vector<ElementWithNodes> boundaryConditionElements;   ///< nodes grouped by elements on which boundary conditions are specified
  std::vector<dof_no_t> boundaryConditionDofLocalNos;        ///< vector of all local boundary condition dofs
  std::vector<double> boundaryConditionValues;               ///< vector of the local prescribed values, related to boundaryConditionDofLocalNos

  element_no_t lastBoundaryConditionElement = -1;

  // get the first dirichlet boundary condition from the list
  std::vector<std::pair<int, double>>::const_iterator boundaryConditionIter = boundaryConditions.begin();

  // loop over local elements
  for (element_no_t elementNoLocal = 0; elementNoLocal < mesh->nElementsLocal(); elementNoLocal++)
  {
    // loop over the nodes of this element
    int elementalDofIndex = 0;
    for (int nodeIndex = 0; nodeIndex < BasisOnMeshType::nNodesPerElement(); nodeIndex++)
    {
      // get global or local nodeNo of current node (depending on inputMeshIsGlobal)
      global_no_t nodeNo = 0;
      if (inputMeshIsGlobal)
      {
        global_no_t elementNoGlobal = mesh->meshPartition()->getElementNoGlobal(elementNoLocal);
        nodeNo = mesh->getNodeNoGlobal(elementNoGlobal, nodeIndex);
      }
      else
      {
        nodeNo = mesh->getNodeNo(elementNoLocal, nodeIndex);
      }

      // loop over dofs of node
      for (int nodalDofIndex = 0; nodalDofIndex < nDofsPerNode; nodalDofIndex++, elementalDofIndex++)
      {
        // check if a dirichlet boundary condition for the current node is given
        while (boundaryConditionIter->first < nodeNo*nDofsPerNode + nodalDofIndex && boundaryConditionIter != boundaryConditions.end())
        {
          boundaryConditionIter++;
        }

        // if the currently considered boundaryCondition entry from config matches the current nodeNo and nodalDofIndex
        if (boundaryConditionIter->first == nodeNo*nDofsPerNode + nodalDofIndex)
        {
          // if there is not yet an entry in boundaryConditionElements with the current element, create one
          if (elementNoLocal != lastBoundaryConditionElement)
          {
            // add current element
            boundaryConditionElements.emplace_back();
            boundaryConditionElements.back().elementNoLocal = elementNoLocal;

            lastBoundaryConditionElement = elementNoLocal;
          }

          // add current node and boundary condition value to list of boundary conditions for current element
          ElementWithNodes &boundaryConditionElement = boundaryConditionElements.back();

          // check if elementalDofIndex already contains the entry for elementalDofIndex
          bool dofIndexAlreadyContained = false;
          for (int i = 0; i < boundaryConditionElement.elementalDofIndex.size(); i++)
          {
            if (boundaryConditionElement.elementalDofIndex[i].first == elementalDofIndex)
            {
              dofIndexAlreadyContained = true;
              break;
            }
          }

          // add current node and boundary condition value to list of boundary conditions for current element
          if (!dofIndexAlreadyContained)
          {
            boundaryConditionElement.elementalDofIndex.push_back(std::pair<int,double>(elementalDofIndex, boundaryConditionIter->second));
          }

          // also store localDofNo
          dof_no_t dofLocalNo = mesh->getDofNo(elementNoLocal, elementalDofIndex);
          boundaryConditionDofLocalNos.push_back(dofLocalNo);
          boundaryConditionValues.push_back(boundaryConditionIter->second);
        }

        // if there are no more boundary conditions, break the loop
        if (boundaryConditionIter == boundaryConditions.end())
          break;
      }

      // if there are no more boundary conditions, break the loop
      if (boundaryConditionIter == boundaryConditions.end())
        break;
    }

    // if there are no more boundary conditions, break the loop
    if (boundaryConditionIter == boundaryConditions.end())
      break;
  }

  if (VLOG_IS_ON(1))
  {
    VLOG(1) << "dofsLocal of BC: " << boundaryConditionDofLocalNos << " with prescribed values: " << boundaryConditionValues;

    for (auto boundaryConditionElement: boundaryConditionElements)
    {
      VLOG(1) << "elementNo: " << boundaryConditionElement.elementNoLocal << " has (dof,value): " << boundaryConditionElement.elementalDofIndex;
    }
  }

  const int D = BasisOnMeshType::dim();
  typedef Quadrature::TensorProduct<D,QuadratureType> QuadratureDD;
  const int nDofsPerElement = BasisOnMeshType::nDofsPerElement();
  typedef MathUtility::Matrix<nDofsPerElement,nDofsPerElement> EvaluationsType;
  typedef std::array<
            EvaluationsType,
            QuadratureDD::numberEvaluations()
          > EvaluationsArrayType;     // evaluations[nGP^D][nDofs][nDofs]
  // setup arrays used for integration
  std::array<std::array<double,D>, QuadratureDD::numberEvaluations()> samplingPoints = QuadratureDD::samplingPoints();
  EvaluationsArrayType evaluationsArray{};

  // loop over elements that have nodes with prescribed boundary conditions
  for (typename std::vector<ElementWithNodes>::const_iterator iter = boundaryConditionElements.cbegin(); iter != boundaryConditionElements.cend(); iter++)
  {
    element_no_t elementNoLocal = iter->elementNoLocal;
    std::array<dof_no_t,nDofsPerElement> dofNosLocal = mesh->getElementDofNosLocal(elementNoLocal);

    // get geometry field (which are the node positions for Lagrange basis and node positions and derivatives for Hermite)
    std::array<Vec3,BasisOnMeshType::nDofsPerElement()> geometry;
    mesh->getElementGeometry(elementNoLocal, geometry);

    // compute integral
    for (unsigned int samplingPointIndex = 0; samplingPointIndex < samplingPoints.size(); samplingPointIndex++)
    {
      // evaluate function to integrate at samplingPoint
      std::array<double,D> xi = samplingPoints[samplingPointIndex];

      // compute the 3xD jacobian of the parameter space to world space mapping
      auto jacobian = BasisOnMeshType::computeJacobian(geometry, xi);

      // get evaluations of integrand at xi for all (i,j)-dof pairs, integrand is defined in another class
      evaluationsArray[samplingPointIndex]
        = IntegrandStiffnessMatrix<D,EvaluationsType,BasisOnMeshType,Term>::evaluateIntegrand(this->data_, jacobian, xi);

    }  // function evaluations

    // integrate all values for the (i,j) dof pairs at once
    EvaluationsType integratedValues = QuadratureDD::computeIntegral(evaluationsArray);

    // perform integration and add to entry of stiffness matrix
    for (std::vector<std::pair<int,double>>::const_iterator dofsIter1 = iter->elementalDofIndex.begin(); dofsIter1 != iter->elementalDofIndex.end(); dofsIter1++)
    {
      int elementalDofIndex1 = dofsIter1->first;
      double boundaryConditionValue = dofsIter1->second;

      for (int j = 0; j < nDofsPerElement; j++)
      {
        // integrate value and set entry in stiffness matrix
        double integratedValue = integratedValues(elementalDofIndex1, j);
        double value = boundaryConditionValue * integratedValue;

        // the formula for boundary conditions is rhs = -sum_{over BC-dofs} u_{0,i} * integral_Omega ∇phi_i • ∇phi_j dx  for equation j=1,...,N

        VLOG(2) << "  dof pair (" << elementalDofIndex1 << "," << j << ") dofs (" << dofNosLocal[elementalDofIndex1] << "," << dofNosLocal[j]
          << "), integrated value: " << integratedValue << ", boundaryConditionValue: " << boundaryConditionValue << ", value: " << value;

        rightHandSide.setValue(dofNosLocal[j], value, ADD_VALUES);
      }  // j
    }
  }

  // set boundary condition dofs to prescribed values
  rightHandSide.setValues(boundaryConditionDofLocalNos,
                          boundaryConditionValues, INSERT_VALUES);
  rightHandSide.finishVectorManipulation();

  // zero entries in stiffness matrix that correspond to dirichlet dofs
  // set values of row and column of the dofs to zero and diagonal entry to 1
  stiffnessMatrix->zeroRowsColumns(boundaryConditionDofLocalNos.size(), boundaryConditionDofLocalNos.data(), 1.0);
  stiffnessMatrix->assembly(MAT_FINAL_ASSEMBLY);

}

};  // namespace
