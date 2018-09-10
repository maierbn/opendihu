#include "spatial_discretization/finite_element_method/02_boundary_conditions.h"

#include <Python.h>
#include <memory>
#include <vector>
#include <petscsys.h>

#include "quadrature/tensor_product.h"
#include "utility/vector_operators.h"
#include "utility/boundary_conditions.h"

namespace SpatialDiscretization
{

template<typename FunctionSpaceType,typename QuadratureType,typename Term,typename Dummy>
void BoundaryConditions<FunctionSpaceType,QuadratureType,Term,Dummy>::
applyBoundaryConditions()
{
  if (VLOG_IS_ON(4))
  {
    VLOG(4) << "Finite element data before applyBoundaryConditions";
    this->data_.print();
  }

  applyBoundaryConditionsWeakForm();

  if (VLOG_IS_ON(4))
  {
    VLOG(4) << "Finite element data after applyBoundaryConditions";
    this->data_.print();
  }
}

template<typename FunctionSpaceType,typename QuadratureType,typename Term,typename Dummy>
void BoundaryConditions<FunctionSpaceType,QuadratureType,Term,Dummy>::
parseBoundaryConditions()
{
  LOG(TRACE) << "parseBoundaryConditions";

  std::shared_ptr<FunctionSpaceType> functionSpace = this->data_.functionSpace();

  // add weak form of Dirichlet BC to rhs
  const int nDofsPerNode = FunctionSpaceType::nDofsPerNode();

  // determine if the BC indices in the config are given for global or local dof nos
  bool inputMeshIsGlobal = PythonUtility::getOptionBool(this->specificSettings_, "inputMeshIsGlobal", true);

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

  std::vector<std::pair<int,double>> boundaryConditions;  // (index, value)
  ::BoundaryConditions::parseBoundaryConditions(this->specificSettings_, functionSpace, boundaryConditions);

  // sort all parsed boundary conditions for their index no
  auto compareFunction = [](const std::pair<int,double> &item1, const std::pair<int,double> &item2)
  {
    return item1.first < item2.first;
  };
  std::sort(boundaryConditions.begin(), boundaryConditions.end(), compareFunction);

  LOG(DEBUG) << "read in boundary conditions from config: " << boundaryConditions;

  // determine elements with nodes that have prescribed Dirichlet boundary conditions, store them in the vector boundaryConditionElements_,
  // which is organized by local elements
  element_no_t lastBoundaryConditionElement = -1;

  // loop over all local elements
  for (element_no_t elementNoLocal = 0; elementNoLocal < functionSpace->nElementsLocal(); elementNoLocal++)
  {
    VLOG(2) << "elementNoLocal: " << elementNoLocal;

    // loop over the nodes of this element
    int elementalDofIndex = 0;
    for (int nodeIndex = 0; nodeIndex < FunctionSpaceType::nNodesPerElement(); nodeIndex++)
    {
      // get global or local nodeNo of current node (depending on inputMeshIsGlobal)
      global_no_t nodeNo = 0;
      if (inputMeshIsGlobal)
      {
        global_no_t elementNoGlobalNatural = functionSpace->meshPartition()->getElementNoGlobalNatural(elementNoLocal);
        nodeNo = functionSpace->getNodeNoGlobalNatural(elementNoGlobalNatural, nodeIndex);
      }
      else
      {
        nodeNo = functionSpace->getNodeNo(elementNoLocal, nodeIndex);
      }

      // loop over dofs of node and thus over the elemental dofs
      for (int nodalDofIndex = 0; nodalDofIndex < nDofsPerNode; nodalDofIndex++, elementalDofIndex++)
      {
        global_no_t indexForBoundaryCondition = nodeNo*nDofsPerNode + nodalDofIndex;

        // check if a dirichlet boundary condition for the current node is given
        bool boundaryConditionForDofFound = false;

        // loop over all stored boundary conditions
        double boundaryConditionValue;
        for (std::vector<std::pair<int, double>>::const_iterator boundaryConditionIter = boundaryConditions.begin();
             boundaryConditionIter != boundaryConditions.end(); boundaryConditionIter++)
        {
          if (boundaryConditionIter->first == indexForBoundaryCondition)
          {
            boundaryConditionValue = boundaryConditionIter->second;
            boundaryConditionForDofFound = true;
            break;
          }
          if (boundaryConditionIter->first >= indexForBoundaryCondition)
            break;
        }

        VLOG(2) << "boundaryConditionForDofFound: " << boundaryConditionForDofFound;

        // here boundaryConditionIter->first is equal to indexForBoundaryCondition (then the current element/node/dof matches the boundary condition)
        // or boundaryConditionIter->first > indexForBoundaryCondition then the current dofIndex does not have any boundary condition
        // or boundaryConditionIter == boundaryConditions.end() then, too

        // if the currently considered boundaryCondition entry from config matches the current nodeNo and nodalDofIndex
        if (boundaryConditionForDofFound)
        {
          VLOG(2) << "elementNoLocal: " << elementNoLocal << ", lastBoundaryConditionElement: " << lastBoundaryConditionElement;

          // if there is not yet an entry in boundaryConditionElements with the current element, create one
          if (elementNoLocal != lastBoundaryConditionElement)
          {
            // add current element
            boundaryConditionElements_.emplace_back();
            boundaryConditionElements_.back().elementNoLocal = elementNoLocal;

            lastBoundaryConditionElement = elementNoLocal;

            VLOG(2) << "add empty entry for elementNoLocal " << elementNoLocal;
          }

          // add current node and boundary condition value to list of boundary conditions for current element
          ElementWithNodes &boundaryConditionElement = boundaryConditionElements_.back();

          /*
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
          {*/
          VLOG(2) << " add (el-dof, value)" << std::pair<int,double>(elementalDofIndex, boundaryConditionValue) << ", to boundaryConditionElement of element " << boundaryConditionElement.elementNoLocal;
            boundaryConditionElement.elementalDofIndex.push_back(std::pair<int,double>(elementalDofIndex, boundaryConditionValue));
          //}

          // also store localDofNo
          dof_no_t dofLocalNo = functionSpace->getDofNo(elementNoLocal, elementalDofIndex);

          if (dofLocalNo < functionSpace->nDofsLocalWithoutGhosts())
          {
            // check if dofLocalNo is already contained in boundaryConditionNonGhostDofLocalNos_
            bool dofLocalNoIsAlreadyContained = false;
            for (auto dofNo : boundaryConditionNonGhostDofLocalNos_)
            {
              if (dofNo == dofLocalNo)
              {
                dofLocalNoIsAlreadyContained = true;
              }
            }
            if (!dofLocalNoIsAlreadyContained)
            {
              boundaryConditionNonGhostDofLocalNos_.push_back(dofLocalNo);
              boundaryConditionValues_.push_back(boundaryConditionValue);
            }
          }
        }
      }
    }
  }

  if (VLOG_IS_ON(1))
  {
    VLOG(1) << "parsed boundary conditions:";
    VLOG(1) << "  dofsLocal of BC: " << boundaryConditionNonGhostDofLocalNos_ << " with prescribed values: " << boundaryConditionValues_;

    for (auto boundaryConditionElement: boundaryConditionElements_)
    {
      VLOG(1) << "  elementNo: " << boundaryConditionElement.elementNoLocal << " has (dof,value): " << boundaryConditionElement.elementalDofIndex;
    }
  }
}

template<typename FunctionSpaceType,typename QuadratureType,typename Term,typename Dummy>
void BoundaryConditions<FunctionSpaceType,QuadratureType,Term,Dummy>::
applyBoundaryConditionsWeakForm()
{

  LOG(TRACE) << "applyBoundaryConditionsWeakForm";

  // parse boundary conditions from config and store them in boundaryConditionElements_, boundaryConditionNonGhostDofLocalNos_ and boundaryConditionValues_
  parseBoundaryConditions();

  // get abbreviations
  std::shared_ptr<FunctionSpaceType> functionSpace = this->data_.functionSpace();
  FieldVariable::FieldVariable<FunctionSpaceType,1> &rightHandSide = this->data_.rightHandSide();
  std::shared_ptr<PartitionedPetscMat<FunctionSpaceType>> stiffnessMatrix = this->data_.stiffnessMatrix();

  rightHandSide.startGhostManipulation();
  rightHandSide.zeroGhostBuffer();

  const int D = FunctionSpaceType::dim();
  typedef Quadrature::TensorProduct<D,QuadratureType> QuadratureDD;
  const int nDofsPerElement = FunctionSpaceType::nDofsPerElement();
  typedef MathUtility::Matrix<nDofsPerElement,nDofsPerElement> EvaluationsType;
  typedef std::array<
            EvaluationsType,
            QuadratureDD::numberEvaluations()
          > EvaluationsArrayType;     // evaluations[nGP^D][nDofs][nDofs]
  // setup arrays used for integration
  std::array<std::array<double,D>, QuadratureDD::numberEvaluations()> samplingPoints = QuadratureDD::samplingPoints();
  EvaluationsArrayType evaluationsArray{};

  // loop over elements that have nodes with prescribed boundary conditions, only for those the integral term is non-zero
  for (typename std::vector<ElementWithNodes>::const_iterator iter = boundaryConditionElements_.cbegin(); iter != boundaryConditionElements_.cend(); iter++)
  {
    element_no_t elementNoLocal = iter->elementNoLocal;
    std::array<dof_no_t,nDofsPerElement> dofNosLocal = functionSpace->getElementDofNosLocal(elementNoLocal);

    // get geometry field (which are the node positions for Lagrange basis and node positions and derivatives for Hermite)
    std::array<Vec3,FunctionSpaceType::nDofsPerElement()> geometry;
    functionSpace->getElementGeometry(elementNoLocal, geometry);

    // compute integral numerically
    for (unsigned int samplingPointIndex = 0; samplingPointIndex < samplingPoints.size(); samplingPointIndex++)
    {
      // evaluate function to integrate at samplingPoint
      std::array<double,D> xi = samplingPoints[samplingPointIndex];

      // compute the 3xD jacobian of the parameter space to world space mapping
      auto jacobian = FunctionSpaceType::computeJacobian(geometry, xi);

      // get evaluations of integrand at xi for all (i,j)-dof pairs, integrand is defined in another class
      evaluationsArray[samplingPointIndex]
        = IntegrandStiffnessMatrix<D,EvaluationsType,FunctionSpaceType,Term>::evaluateIntegrand(this->data_, jacobian, xi);

    }  // function evaluations

    // integrate all values for the (i,j) dof pairs at once
    EvaluationsType integratedValues = QuadratureDD::computeIntegral(evaluationsArray);

    // perform integration
    for (std::vector<std::pair<int,double>>::const_iterator dofsIter1 = iter->elementalDofIndex.begin(); dofsIter1 != iter->elementalDofIndex.end(); dofsIter1++)
    {
      int elementalDofIndex1 = dofsIter1->first;
      double boundaryConditionValue = dofsIter1->second;

      for (int j = 0; j < nDofsPerElement; j++)
      {
        // integrate value and set entry in stiffness matrix
        double integratedValue = integratedValues(elementalDofIndex1, j);
        double value = boundaryConditionValue * integratedValue;

        // the formula for boundary conditions is rhs = -sum_{over BC-dofs i} u_{0,i} * integral_Omega ∇phi_i • ∇phi_j dx  for equation j=1,...,N

        VLOG(2) << "  dof pair (" << elementalDofIndex1 << "," << j << ") dofs (" << dofNosLocal[elementalDofIndex1] << "," << dofNosLocal[j]
          << "), integrated value: " << integratedValue << ", boundaryConditionValue: " << boundaryConditionValue << ", value: " << value;

        rightHandSide.setValue(dofNosLocal[j], value, ADD_VALUES);
      }  // j
    }
  }

  // Scatter the ghost values to their actual place, with ADD_VALUES,
  // so the ghost value on one rank and the non-ghost value on the other rank are added.
  rightHandSide.finishGhostManipulation();

  // set boundary condition dofs to prescribed values, only non-ghost dofs
  rightHandSide.setValues(boundaryConditionNonGhostDofLocalNos_,
                          boundaryConditionValues_, INSERT_VALUES);

  // zero entries in stiffness matrix that correspond to dirichlet dofs
  // set values of row and column of the dofs to zero and diagonal entry to 1
  stiffnessMatrix->zeroRowsColumns(boundaryConditionNonGhostDofLocalNos_.size(), boundaryConditionNonGhostDofLocalNos_.data(), 1.0);
  stiffnessMatrix->assembly(MAT_FINAL_ASSEMBLY);

}

};  // namespace
