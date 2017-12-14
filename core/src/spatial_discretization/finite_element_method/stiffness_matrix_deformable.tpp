
#include <Python.h>
#include <memory>
#include <vector>
#include <petscksp.h>
#include <petscsys.h>
#include <cmath>

#include "mesh/deformable.h"
#include "basis_function/lagrange.h"

#include "spatial_discretization/spatial_discretization.h"
#include "time_stepping_scheme/discretizable_in_time.h"
#include "control/runnable.h"
#include "control/dihu_context.h"
#include "utility/petsc_utility.h"
#include "utility/math_utility.h"
#include "data_management/finite_elements.h"
#include "equation/laplace.h"
#include "equation/poisson.h"
#include "equation/type_traits.h"
#include "mesh/mesh.h"
#include "integrator/tensor_product.h"
#include "basis_function/tensor_product.h"

namespace SpatialDiscretization
{
  
// 1D stiffness matrix
template<typename BasisFunctionType, typename IntegratorType>
void FiniteElementMethodStiffnessMatrix<Mesh::Deformable<1ul>, BasisFunctionType, IntegratorType>::
setStiffnessMatrix()
{
  LOG(TRACE)<<"setStiffnessMatrix 1D";
 
  std::shared_ptr<Mesh::Deformable<1>> mesh = std::static_pointer_cast<Mesh::Deformable<1>>(this->data_.mesh());
  
  // get settings values
  double prefactor = PythonUtility::getOptionDouble(this->specificSettings_, "prefactor", 1.0);
  
  // fill stiffness matrix
  // M_ij = -int[0,1] dphi_i/dxi * dphi_j/dxi * s^{-2} * J1 * ds
  PetscErrorCode ierr;
  Mat &stiffnessMatrix = this->data_.stiffnessMatrix();
  
  typedef BasisFunction::TensorProduct<1,BasisFunctionType> Basis1D;
  
  int nDofsPerElement = Basis1D::nDofsPerElement();
  int nElements = mesh->nElements();
  std::array<int,1> nElementsArray{nElements};
  
  // node numbers for linear Lagrange basis: 
  // 0 1
  
  // initialize values to zero
  // loop over elements 
  for (element_idx_t elementNo = 0; elementNo < mesh->nElements(); elementNo++)
  {
    //std::array<int,nDofsPerElement>
    auto dof = Basis1D::getElementDofs(elementNo, nElementsArray);
    // for linear Lagrange basis, the previous loop has the same effect as the following lines:
    //int dof0 = elementNo;  
    //int dof1 = elementNo+1;
    
    for (int i=0; i<nDofsPerElement; i++)
    {
      for (int j=0; j<nDofsPerElement; j++)
      {
        ierr = MatSetValue(stiffnessMatrix, dof[i], dof[j], 0, INSERT_VALUES); CHKERRV(ierr);
      }
    }
  }

  // setup arrays used for integration
  std::array<double, IntegratorType::numberEvaluations()> samplingPoints = IntegratorType::samplingPoints();
  std::array<double, IntegratorType::numberEvaluations()> evaluations[nDofsPerElement][nDofsPerElement];
  
  // loop over elements 
  for (element_idx_t elementNo = 0; elementNo < mesh->nElements(); elementNo++)
  {
    // get nodes and element length
    //int dof0 = elementNo;
    //int dof1 = elementNo+1;
    
    auto dof = Basis1D::getElementDofs(elementNo, nElementsArray);
    
    std::array<Vec3,Basis1D::nDofsPerElement()> node;
    for (int dofIndex = 0; dofIndex < nDofsPerElement; dofIndex++)
    {
      // get the global dof no. from element-local dofIndex
      node[dofIndex] = mesh->getNodePosition(dof[dofIndex]);
    }
    //Vec3 node0 = mesh->getNodePosition(dof0);
    //Vec3 node1 = mesh->getNodePosition(dof1);
    
    double elementLength = MathUtility::distance(node[0], node[nDofsPerElement-1]);
    double integralFactor = 1. / elementLength;
    
    // compute integral
    // loop over sampling points (e.g. Gauss points, depending on quadrature scheme)
    for (unsigned int samplingPointIndex=0; samplingPointIndex<samplingPoints.size(); samplingPointIndex++)
    {
      // evaluate function to integrate at samplingPoints[i*D], write value to evaluations[i]
      double xi = samplingPoints[samplingPointIndex];
      
      // loop over pairs of dofs of the current element
      for (int i=0; i<nDofsPerElement; i++)
      {
        for (int j=0; j<nDofsPerElement; j++)
        {
          double integrand = BasisFunctionType::dphi_dxi(i,xi) * BasisFunctionType::dphi_dxi(j,xi) * integralFactor;
          evaluations[i][j][samplingPointIndex] = integrand;
        }
      }
    }
    
    // perform quadrature from stored function evaluations and store result to matrix
    // loop over pairs of dofs of the current element
    for (int i=0; i<nDofsPerElement; i++)
    {
      for (int j=0; j<nDofsPerElement; j++)
      {
        double value = -prefactor*IntegratorType::integrate(evaluations[i][j]);
        ierr = MatSetValue(stiffnessMatrix, dof[i], dof[j], value, ADD_VALUES); CHKERRV(ierr);
      }
    }
  }
}

// 2D stiffness matrix
template<typename BasisFunctionType, typename IntegratorType>
void FiniteElementMethodStiffnessMatrix<Mesh::Deformable<2>, BasisFunctionType, IntegratorType>::
setStiffnessMatrix()
{
  LOG(TRACE)<<"setStiffnessMatrix 2D";
 
  // get settings values
  std::shared_ptr<Mesh::Deformable<2>> mesh = std::static_pointer_cast<Mesh::Deformable<2>>(this->data_.mesh());
  const int nElements0 = mesh->nElements(0);
  const int nElements1 = mesh->nElements(1);
  const int nNodes0 = nElements0 + 1;
  const int nNodes1 = nElements1 + 1;
  const double prefactor = PythonUtility::getOptionDouble(this->specificSettings_, "prefactor", 1.0);
  
  // fill stiffness matrix
  PetscErrorCode ierr;
  Mat &stiffnessMatrix = this->data_.stiffnessMatrix();
  
  typedef Integrator::TensorProduct<2,IntegratorType> Integrator2D;
  typedef BasisFunction::TensorProduct<2,BasisFunctionType> Basis2D;
  
  int nDofsPerElement = Basis2D::nDofsPerElement();
  std::array<int,2> nElementsArray{nElements0, nElements1};
    
  LOG(DEBUG) << "nNodes: " << nNodes0 << "x"<<nNodes1<<", prefactor="<<prefactor;
  
  // node numbers for linear Lagrange basis:
  // 2 3
  // 0 1 
  
  // initialize values to zero
  // loop over elements 
  for (int elementNo = 0; elementNo < mesh->nElements(); elementNo++)
  {
    auto dof = Basis2D::getElementDofs(elementNo, nElementsArray);
    
    for (int i=0; i<nDofsPerElement; i++)
    {
      for (int j=0; j<nDofsPerElement; j++)
      {
        ierr = MatSetValue(stiffnessMatrix, dof[i], dof[j], 0, INSERT_VALUES); CHKERRV(ierr);
      }
    }
  }
  
  // ansatz functions used in integral
  /*
  const auto phi0 = [](double xi){return 1-xi;};
  const auto phi1 = [](double xi){return xi;};
  const auto dphi0_dxi = [](double xi){return -1;};
  const auto dphi1_dxi = [](double xi){return 1;};
  
  const auto dphi0_dxi1 = [&](double xi1, double xi2){return dphi0_dxi(xi1)*phi0(xi2);};
  const auto dphi0_dxi2 = [&](double xi1, double xi2){return phi0(xi1)*dphi0_dxi(xi2);};
  
  const auto dphi1_dxi1 = [&](double xi1, double xi2){return dphi1_dxi(xi1)*phi0(xi2);};
  const auto dphi1_dxi2 = [&](double xi1, double xi2){return phi1(xi1)*dphi0_dxi(xi2);};
  
  const auto dphi2_dxi1 = [&](double xi1, double xi2){return dphi0_dxi(xi1)*phi1(xi2);};
  const auto dphi2_dxi2 = [&](double xi1, double xi2){return phi0(xi1)*dphi1_dxi(xi2);};
  
  const auto dphi3_dxi1 = [&](double xi1, double xi2){return dphi1_dxi(xi1)*phi1(xi2);};
  const auto dphi3_dxi2 = [&](double xi1, double xi2){return phi1(xi1)*dphi1_dxi(xi2);};
   */
  
  // setup arrays used for integration
  std::array<double, Integrator2D::samplingArraySize()> samplingPoints = Integrator2D::samplingPoints();
  std::array<double, Integrator2D::numberEvaluations()> evaluations[nDofsPerElement][nDofsPerElement];
  
  LOG(DEBUG) << "1D integration with " << IntegratorType::numberEvaluations() << " evaluations";
  LOG(DEBUG) << "2D integration with " << Integrator2D::numberEvaluations() << " evaluations";
#ifdef DEBUG
  LOG(DEBUG) << "SAMPLING POINTS: ";
  for  (auto value : samplingPoints)
    LOG(DEBUG) << "   " << value;
#endif
  
  // loop over elements 
  for (int elementNo = 0; elementNo < mesh->nElements(); elementNo++)
  {
    // get indices of element-local dofs
    auto dof = Basis2D::getElementDofs(elementNo, nElementsArray);
    
    // get node positions for dofs
    std::array<Vec3,Basis2D::nDofsPerElement()> node;
    for (int dofIndex = 0; dofIndex < nDofsPerElement; dofIndex++)
    {
      node[dofIndex] = mesh->getNodePosition(dof[dofIndex]);
    }
    
    // compute integral
    for (unsigned int samplingPointIndex=0; samplingPointIndex<samplingPoints.size()/2; samplingPointIndex++)
    {
      // evaluate function to integrate at samplingPoints[i*2], write value to evaluations[i]
      double xi1 = samplingPoints[samplingPointIndex*2+0];
      double xi2 = samplingPoints[samplingPointIndex*2+1];
      
      // compute the 3x2 jacobian of the parameter space to world space mapping
      //Vec3 jacobianColumn0 = (1-xi2) * (node[1]-node[0]) + xi2 * (node[3]-node[2]);
      //Vec3 jacobianColumn1 = (1-xi2) * (node[2]-node[0]) + xi2 * (node[3]-node[1]);
      auto jacobian = Basis2D::computeJacobian(node, std::array<double,2>({xi1,xi2}));
      Vec3 jacobianColumn0 = jacobian[0];
      Vec3 jacobianColumn1 = jacobian[1];
      
      Vec3 &zeta1 = jacobianColumn0;
      Vec3 &zetah = jacobianColumn1;
      
      double integrationFactor = MathUtility::compute2DIntegrationFactor(zeta1, zetah);
      
      double l1 = MathUtility::length(zeta1);
      double lh = MathUtility::length(zetah);
      double beta = acos((zeta1[0]*zetah[0] + zeta1[1]*zetah[1] + zeta1[2]*zetah[2]) / (l1 * lh));
      double alpha = MathUtility::sqr(MathUtility::PI) / (4.*beta);
      
      Vec3 zeta2 = cos(alpha) * zeta1 + sin(alpha) * zetah;
      double l2 = MathUtility::length(zeta2);
      double l1squared = MathUtility::sqr(l1);
      double l2squared = MathUtility::sqr(l2);
      
      // compute the 3x3 transformation matrix T = J^{-1}J^{-T} and the absolute of the determinant of the jacobian
      std::array<double,4> transformationMatrix = {
        MathUtility::sqr(cos(alpha))/l2squared + 1./l1squared,
        sin(alpha)*cos(alpha)/l2squared,
        sin(alpha)*cos(alpha)/l2squared,
        MathUtility::sqr(sin(alpha))/l2squared
      };

#ifdef DEBUG        
      VLOG(3) << "transformationMatrix:";
      std::stringstream s;
      for(int i=0; i<2; i++)
      {
        for(int j=0; j<2; j++)
        {
          s << transformationMatrix[i*2+j] << " ";
        }
        s << std::endl;
      }
      VLOG(3) << s.str(); 
#endif          
      
      // initialize gradient vectors of ansatz function phi_i, for node i of current element
      std::array<Vec2,Basis2D::nDofsPerElement()> gradPhi = Basis2D::getGradPhi(std::array<double,2>({xi1,xi2}));
      
      //gradPhi[0] = {dphi0_dxi1(xi1,xi2), dphi0_dxi2(xi1,xi2)};
      //gradPhi[1] = {dphi1_dxi1(xi1,xi2), dphi1_dxi2(xi1,xi2)};
      //gradPhi[2] = {dphi2_dxi1(xi1,xi2), dphi2_dxi2(xi1,xi2)};
      //gradPhi[3] = {dphi3_dxi1(xi1,xi2), dphi3_dxi2(xi1,xi2)};
      
      //LOG(DEBUG) << "XI: " << xi1<<","<<xi2<<","<<xi3;
      //LOG(DEBUG) << "gradPhi0:" << gradPhi[0][0] << "," << gradPhi[0][1] << ","<< gradPhi[0][2];
      
      for (int i=0; i<nDofsPerElement; i++)
      {
        for (int j=0; j<nDofsPerElement; j++)
        {
          double integrand = MathUtility::applyTransformation(transformationMatrix, gradPhi[i], gradPhi[j]) * integrationFactor;
          evaluations[i][j][samplingPointIndex] = integrand;
        }
      }
    }  // function evaluations
    
    // perform integration and add to entry of stiffness matrix
    for (int i=0; i<nDofsPerElement; i++)
    {
      for (int j=0; j<nDofsPerElement; j++)
      {
        double value = -prefactor * Integrator2D::integrate(evaluations[i][j]);
        ierr = MatSetValue(stiffnessMatrix, dof[i], dof[j], value, ADD_VALUES); CHKERRV(ierr);
      }  // j
    }  // i
  }  // elementNo
}
  
// 3D stiffness matrix
template<typename BasisFunctionType, typename IntegratorType>
void FiniteElementMethodStiffnessMatrix<Mesh::Deformable<3>, BasisFunctionType, IntegratorType>::
setStiffnessMatrix()
{
  LOG(TRACE)<<"setStiffnessMatrix 3D";
 
  // get settings values
  std::shared_ptr<Mesh::Deformable<3>> mesh = std::static_pointer_cast<Mesh::Deformable<3>>(this->data_.mesh());
  const int nElements0 = mesh->nElements(0);
  const int nElements1 = mesh->nElements(1);
  const int nElements2 = mesh->nElements(2);
  const int nNodes0 = nElements0 + 1;
  const int nNodes1 = nElements1 + 1;
  const int nNodes2 = nElements2 + 1;
  const double prefactor = PythonUtility::getOptionDouble(this->specificSettings_, "prefactor", 1.0);
  
  // fill stiffness matrix
  PetscErrorCode ierr;
  Mat &stiffnessMatrix = this->data_.stiffnessMatrix();
  
  typedef Integrator::TensorProduct<3,IntegratorType> Integrator3D;
  typedef BasisFunction::TensorProduct<3,BasisFunctionType> Basis3D;
  
  int nDofsPerElement = Basis3D::nDofsPerElement();
  std::array<int,3> nElementsArray{nElements0, nElements1, nElements2};
    
  LOG(DEBUG) << "nNodes: " << nNodes0 << "x"<<nNodes1<<"x"<<nNodes2<<", prefactor="<<prefactor;
  
  // node numbers: 
  // 6 7  (top)
  // 4 5
  
  // 2 3  (bottom)
  // 0 1 
  
  // initialize values to zero
  // loop over elements 
  for (int elementNo = 0; elementNo < mesh->nElements(); elementNo++)
  {
    auto dof = Basis3D::getElementDofs(elementNo, nElementsArray);

    for (int i=0; i<nDofsPerElement; i++)
    {
      for (int j=0; j<nDofsPerElement; j++)
      {
        ierr = MatSetValue(stiffnessMatrix, dof[i], dof[j], 0, INSERT_VALUES); CHKERRV(ierr);
      }
    }
  }
  
  // ansatz functions used in integral
  /*
  const auto phi0 = [](double xi){return 1-xi;};
  const auto phi1 = [](double xi){return xi;};
  const auto dphi0_dxi = [](double xi){return -1;};
  const auto dphi1_dxi = [](double xi){return 1;};
  
  const auto dphi0_dxi1 = [&](double xi1, double xi2, double xi3){return dphi0_dxi(xi1)*phi0(xi2)*phi0(xi3);};
  const auto dphi0_dxi2 = [&](double xi1, double xi2, double xi3){return phi0(xi1)*dphi0_dxi(xi2)*phi0(xi3);};
  const auto dphi0_dxi3 = [&](double xi1, double xi2, double xi3){return phi0(xi1)*phi0(xi2)*dphi0_dxi(xi3);};
  
  const auto dphi1_dxi1 = [&](double xi1, double xi2, double xi3){return dphi1_dxi(xi1)*phi0(xi2)*phi0(xi3);};
  const auto dphi1_dxi2 = [&](double xi1, double xi2, double xi3){return phi1(xi1)*dphi0_dxi(xi2)*phi0(xi3);};
  const auto dphi1_dxi3 = [&](double xi1, double xi2, double xi3){return phi1(xi1)*phi0(xi2)*dphi0_dxi(xi3);};
  
  const auto dphi2_dxi1 = [&](double xi1, double xi2, double xi3){return dphi0_dxi(xi1)*phi1(xi2)*phi0(xi3);};
  const auto dphi2_dxi2 = [&](double xi1, double xi2, double xi3){return phi0(xi1)*dphi1_dxi(xi2)*phi0(xi3);};
  const auto dphi2_dxi3 = [&](double xi1, double xi2, double xi3){return phi0(xi1)*phi1(xi2)*dphi0_dxi(xi3);};
  
  const auto dphi3_dxi1 = [&](double xi1, double xi2, double xi3){return dphi1_dxi(xi1)*phi1(xi2)*phi0(xi3);};
  const auto dphi3_dxi2 = [&](double xi1, double xi2, double xi3){return phi1(xi1)*dphi1_dxi(xi2)*phi0(xi3);};
  const auto dphi3_dxi3 = [&](double xi1, double xi2, double xi3){return phi1(xi1)*phi1(xi2)*dphi0_dxi(xi3);};
  
  const auto dphi4_dxi1 = [&](double xi1, double xi2, double xi3){return dphi0_dxi(xi1)*phi0(xi2)*phi1(xi3);};
  const auto dphi4_dxi2 = [&](double xi1, double xi2, double xi3){return phi0(xi1)*dphi0_dxi(xi2)*phi1(xi3);};
  const auto dphi4_dxi3 = [&](double xi1, double xi2, double xi3){return phi0(xi1)*phi0(xi2)*dphi1_dxi(xi3);};
  
  const auto dphi5_dxi1 = [&](double xi1, double xi2, double xi3){return dphi1_dxi(xi1)*phi0(xi2)*phi1(xi3);};
  const auto dphi5_dxi2 = [&](double xi1, double xi2, double xi3){return phi1(xi1)*dphi0_dxi(xi2)*phi1(xi3);};
  const auto dphi5_dxi3 = [&](double xi1, double xi2, double xi3){return phi1(xi1)*phi0(xi2)*dphi1_dxi(xi3);};
  
  const auto dphi6_dxi1 = [&](double xi1, double xi2, double xi3){return dphi0_dxi(xi1)*phi1(xi2)*phi1(xi3);};
  const auto dphi6_dxi2 = [&](double xi1, double xi2, double xi3){return phi0(xi1)*dphi1_dxi(xi2)*phi1(xi3);};
  const auto dphi6_dxi3 = [&](double xi1, double xi2, double xi3){return phi0(xi1)*phi1(xi2)*dphi1_dxi(xi3);};
  
  const auto dphi7_dxi1 = [&](double xi1, double xi2, double xi3){return dphi1_dxi(xi1)*phi1(xi2)*phi1(xi3);};
  const auto dphi7_dxi2 = [&](double xi1, double xi2, double xi3){return phi1(xi1)*dphi1_dxi(xi2)*phi1(xi3);};
  const auto dphi7_dxi3 = [&](double xi1, double xi2, double xi3){return phi1(xi1)*phi1(xi2)*dphi1_dxi(xi3);};
  */
  
  // setup arrays used for integration
  std::array<double, Integrator3D::samplingArraySize()> samplingPoints = Integrator3D::samplingPoints();
  std::array<double, Integrator3D::numberEvaluations()> evaluations[nDofsPerElement][nDofsPerElement];
  
  LOG(DEBUG) << "1D integration with " << IntegratorType::numberEvaluations() << " evaluations";
  LOG(DEBUG) << "3D integration with " << Integrator3D::numberEvaluations() << " evaluations";
#ifdef DEBUG
  LOG(DEBUG) << "SAMPLING POINTS: ";
  for  (auto value : samplingPoints)
    LOG(DEBUG) << "   " << value;
#endif
  
  // loop over elements 
  for (int elementNo = 0; elementNo < mesh->nElements(); elementNo++)
  {
    // get global dof index of element
    auto dof = Basis3D::getElementDofs(elementNo, nElementsArray);
    
    // get node position of dofs
    std::array<Vec3,Basis3D::nDofsPerElement()> node;
    for (int dofIndex = 0; dofIndex < nDofsPerElement; dofIndex++)
    {
      node[dofIndex] = mesh->getNodePosition(dof[dofIndex]);
    }
    
    // compute integral
    for (unsigned int samplingPointIndex=0; samplingPointIndex<samplingPoints.size()/3; samplingPointIndex++)
    {
      // evaluate function to integrate at samplingPoints[i*3], write value to evaluations[i]
      double xi1 = samplingPoints[samplingPointIndex*3+0];
      double xi2 = samplingPoints[samplingPointIndex*3+1];
      double xi3 = samplingPoints[samplingPointIndex*3+2];
      
      // compute the 3x3 jacobian of the parameter space to world space mapping
      auto jacobian = Basis3D::computeJacobian(node, std::array<double,3>({xi1,xi2,xi3}));
      Vec3 jacobianColumn0 = jacobian[0];
      Vec3 jacobianColumn1 = jacobian[1];
      Vec3 jacobianColumn2 = jacobian[2];
        
      // compute the 3x3 transformation matrix T = J^{-1}J^{-T} and the absolute of the determinant of the jacobian
      double determinant;
      std::array<double, 9> transformationMatrix 
      = MathUtility::computeTransformationMatrixAndDeterminant(jacobianColumn0, jacobianColumn1, jacobianColumn2, determinant);

#ifdef DEBUG        
      VLOG(3) << "transformationMatrix:";
      std::stringstream s;
      for(int i=0; i<3; i++)
      {
        for(int j=0; j<3; j++)
        {
          s << transformationMatrix[i*3+j] << " ";
        }
        s << std::endl;
      }
      VLOG(3) << s.str(); 
#endif          
      
      // initialize gradient vectors of ansatz function phi_i, for node i of current element
      std::array<Vec3,Basis3D::nDofsPerElement()> gradPhi = Basis3D::getGradPhi(std::array<double,3>({xi1,xi2,xi3}));
    
      // initialize gradient vectors of ansatz function phi_i, for node i of current element
      /*Vec3 gradPhi[8];
      gradPhi[0] = {dphi0_dxi1(xi1,xi2,xi3), dphi0_dxi2(xi1,xi2,xi3), dphi0_dxi3(xi1,xi2,xi3)};
      gradPhi[1] = {dphi1_dxi1(xi1,xi2,xi3), dphi1_dxi2(xi1,xi2,xi3), dphi1_dxi3(xi1,xi2,xi3)};
      gradPhi[2] = {dphi2_dxi1(xi1,xi2,xi3), dphi2_dxi2(xi1,xi2,xi3), dphi2_dxi3(xi1,xi2,xi3)};
      gradPhi[3] = {dphi3_dxi1(xi1,xi2,xi3), dphi3_dxi2(xi1,xi2,xi3), dphi3_dxi3(xi1,xi2,xi3)};
      gradPhi[4] = {dphi4_dxi1(xi1,xi2,xi3), dphi4_dxi2(xi1,xi2,xi3), dphi4_dxi3(xi1,xi2,xi3)};
      gradPhi[5] = {dphi5_dxi1(xi1,xi2,xi3), dphi5_dxi2(xi1,xi2,xi3), dphi5_dxi3(xi1,xi2,xi3)};
      gradPhi[6] = {dphi6_dxi1(xi1,xi2,xi3), dphi6_dxi2(xi1,xi2,xi3), dphi6_dxi3(xi1,xi2,xi3)};
      gradPhi[7] = {dphi7_dxi1(xi1,xi2,xi3), dphi7_dxi2(xi1,xi2,xi3), dphi7_dxi3(xi1,xi2,xi3)};
      */
      
      //LOG(DEBUG) << "XI: " << xi1<<","<<xi2<<","<<xi3;
      //LOG(DEBUG) << "gradPhi0:" << gradPhi[0][0] << "," << gradPhi[0][1] << ","<< gradPhi[0][2];
      
      for (int i=0; i<nDofsPerElement; i++)
      {
        for (int j=0; j<nDofsPerElement; j++)
        {
          double integrand = MathUtility::applyTransformation(transformationMatrix, gradPhi[i], gradPhi[j]) * fabs(determinant);
          evaluations[i][j][samplingPointIndex] = integrand;
        }
      }
    }  // function evaluations
    
    // perform integration and add to entry of stiffness matrix
    for (int i=0; i<nDofsPerElement; i++)
    {
      for (int j=0; j<nDofsPerElement; j++)
      {
        double value = -prefactor * Integrator3D::integrate(evaluations[i][j]);
        ierr = MatSetValue(stiffnessMatrix, dof[i], dof[j], value, ADD_VALUES); CHKERRV(ierr);
      }  // j
    }  // i
  }  // elementNo
}

// 1D rhs
template<typename BasisFunctionType, typename IntegratorType>
void FiniteElementMethodBaseRhs<Mesh::Deformable<1>, BasisFunctionType, IntegratorType>::
transferRhsToWeakForm()
{
  LOG(TRACE)<<"transferRhsToWeakForm (1D)";
 
  std::shared_ptr<Mesh::Deformable<1>> mesh = std::static_pointer_cast<Mesh::Deformable<1>>(this->data_.mesh());
  
  // get settings values
  int nElements = mesh->nElements();
  PetscErrorCode ierr;
  Vec &rightHandSide = this->data_.rightHandSide();
  
  typedef BasisFunction::TensorProduct<1,BasisFunctionType> Basis1D;
  
  int nDofsPerElement = Basis1D::nDofsPerElement();
  std::array<int,1> nElementsArray{nElements};
  
  // get all entries
  std::vector<double> rhsValues;
  PetscUtility::getVectorEntries(rightHandSide, rhsValues);
  
  // initialize values to zero
  VecZeroEntries(rightHandSide);

  // setup arrays used for integration
  std::array<double, IntegratorType::numberEvaluations()> samplingPoints = IntegratorType::samplingPoints();
  std::array<double, IntegratorType::numberEvaluations()> evaluations[nDofsPerElement][nDofsPerElement];
  
  // node numbers for linear Lagrange basis: 
  // 0 1
  
  // loop over elements 
  for (element_idx_t elementNo = 0; elementNo < mesh->nElements(); elementNo++)
  {
    // get nodes and element length
    auto dof = Basis1D::getElementDofs(elementNo, nElementsArray);
    
    std::array<Vec3,Basis1D::nDofsPerElement()> node;
    for (int dofIndex = 0; dofIndex < nDofsPerElement; dofIndex++)
    {
      // get the global dof no. from element-local dofIndex
      node[dofIndex] = mesh->getNodePosition(dof[dofIndex]);
    }
    
    double elementLength = MathUtility::distance(node[0], node[nDofsPerElement-1]);
    double integralFactor = elementLength;
    
    // compute integral
    // loop over sampling points (e.g. Gauss points, depending on quadrature scheme)
    for (unsigned int samplingPointIndex=0; samplingPointIndex<samplingPoints.size(); samplingPointIndex++)
    {
      // evaluate function to integrate at samplingPoints[i*D], write value to evaluations[i]
      double xi = samplingPoints[samplingPointIndex];
      
      // loop over pairs of dofs of the current element
      for (int i=0; i<nDofsPerElement; i++)
      {
        for (int j=0; j<nDofsPerElement; j++)
        {
          double integrand = BasisFunctionType::phi(i,xi) * BasisFunctionType::phi(j,xi) * integralFactor;
          evaluations[i][j][samplingPointIndex] = integrand;
        }
      }
    }
        
    // perform quadrature from stored function evaluations and store result to matrix
    // loop over pairs of dofs of the current element
    for (int i=0; i<nDofsPerElement; i++)
    {
      for (int j=0; j<nDofsPerElement; j++)
      {
        double value = IntegratorType::integrate(evaluations[i][j]);
        ierr = VecSetValue(rightHandSide, dof[i], value*rhsValues[dof[j]], ADD_VALUES); CHKERRV(ierr);
      }
    }
  }
  
  VecAssemblyBegin(rightHandSide);
  VecAssemblyEnd(rightHandSide);
}

// 2D rhs
template<typename BasisFunctionType, typename IntegratorType>
void FiniteElementMethodBaseRhs<Mesh::Deformable<2>, BasisFunctionType, IntegratorType>::
transferRhsToWeakForm()
{
  LOG(TRACE)<<"transferRhsToWeakForm (2D)";
 
  std::shared_ptr<Mesh::Deformable<2>> mesh = std::static_pointer_cast<Mesh::Deformable<2>>(this->data_.mesh());
  const int nElements0 = mesh->nElements(0);
  const int nElements1 = mesh->nElements(1);
  const int nNodes0 = nElements0 + 1;
  const int nNodes1 = nElements1 + 1;
  
  int nElements = mesh->nElements();
  PetscErrorCode ierr;
  Vec &rightHandSide = this->data_.rightHandSide();
  
  typedef Integrator::TensorProduct<2,IntegratorType> Integrator2D;
  typedef BasisFunction::TensorProduct<2,BasisFunctionType> Basis2D;
  
  int nDofsPerElement = Basis2D::nDofsPerElement();
  std::array<int,2> nElementsArray{nElements0, nElements1};
    
  // get all entries
  std::vector<double> rhsValues;
  PetscUtility::getVectorEntries(rightHandSide, rhsValues);
  
  // initialize values to zero
  VecZeroEntries(rightHandSide);
  
  // setup arrays used for integration
  std::array<double, Integrator2D::samplingArraySize()> samplingPoints = Integrator2D::samplingPoints();
  std::array<double, Integrator2D::numberEvaluations()> evaluations[nDofsPerElement][nDofsPerElement];
  
  // node numbers for linear Lagrange basis:
  // 2 3
  // 0 1 
  
  // loop over elements
  for (int elementNo = 0; elementNo < mesh->nElements(); elementNo++)
  {
    auto dof = Basis2D::getElementDofs(elementNo, nElementsArray);
    
    Vec3 node[nDofsPerElement];
    for (int dofIndex = 0; dofIndex < nDofsPerElement; dofIndex++)
    {
      // get the global dof no. from element-local dofIndex
      node[dofIndex] = mesh->getNodePosition(dof[dofIndex]);
    }
    
    // compute integral
    for (unsigned int samplingPointIndex=0; samplingPointIndex<samplingPoints.size()/2; samplingPointIndex++)
    {
      // evaluate function to integrate at samplingPoints[i*2], write value to evaluations[i]
      double xi1 = samplingPoints[samplingPointIndex*2+0];
      double xi2 = samplingPoints[samplingPointIndex*2+1];
              
      // compute the 3x2 jacobian of the parameter space to world space mapping
      //Vec3 jacobianColumn0 = (1-xi2) * (node[1]-node[0]) + xi2 * (node[3]-node[2]);
      //Vec3 jacobianColumn1 = (1-xi2) * (node[2]-node[0]) + xi2 * (node[3]-node[1]);
      auto jacobian = Basis2D::computeJacobian(node, std::array<double,2>({xi1,xi2}));
      Vec3 jacobianColumn0 = jacobian[0];
      Vec3 jacobianColumn1 = jacobian[1];
      
      double integrationFactor = MathUtility::compute2DIntegrationFactor(jacobianColumn0, jacobianColumn1);
      
      for (int i=0; i<nDofsPerElement; i++)
      {
        for (int j=0; j<nDofsPerElement; j++)
        {
          double integrand = Basis2D::phi(i,std::array<double,2>({xi1,xi2})) * Basis2D::phi(j,std::array<double,2>({xi1,xi2}))
            * integrationFactor;
          evaluations[i][j][samplingPointIndex] = integrand;
        }
      }
    }  // function evaluations
    
    // perform integration and add to entry of stiffness matrix
    for (int i=0; i<nDofsPerElement; i++)
    {
      for (int j=0; j<nDofsPerElement; j++)
      {
        double value = Integrator2D::integrate(evaluations[i][j]);
        ierr = VecSetValue(rightHandSide, dof[i], value*rhsValues[dof[j]], ADD_VALUES); CHKERRV(ierr);
      }  // j
    }  // i
  }  // elementNo
  
  VecAssemblyBegin(rightHandSide);
  VecAssemblyEnd(rightHandSide);
}

// 3D rhs
template<typename BasisFunctionType, typename IntegratorType>
void FiniteElementMethodBaseRhs<Mesh::Deformable<3>, BasisFunctionType, IntegratorType>::
transferRhsToWeakForm()
{
  LOG(DEBUG)<<"transferRhsToWeakForm (3D)";
  std::shared_ptr<Mesh::Deformable<3>> mesh = std::static_pointer_cast<Mesh::Deformable<3>>(this->data_.mesh());

  // get settings values
  const int nElements0 = mesh->nElements(0);
  const int nElements1 = mesh->nElements(1);
  const int nElements2 = mesh->nElements(2);
  const int nNodes0 = nElements0 + 1;
  const int nNodes1 = nElements1 + 1;
  const int nNodes2 = nElements2 + 1;
  
  int nElements = mesh->nElements();
  PetscErrorCode ierr;
  Vec &rightHandSide = this->data_.rightHandSide();
  
  typedef Integrator::TensorProduct<3,IntegratorType> Integrator3D;
  typedef BasisFunction::TensorProduct<3,BasisFunctionType> Basis3D;
  
  int nDofsPerElement = Basis3D::nDofsPerElement();
  std::array<int,3> nElementsArray{nElements0, nElements1, nElements2};
    
  // get all entries
  std::vector<double> rhsValues;
  PetscUtility::getVectorEntries(rightHandSide, rhsValues);
  
  // initialize values to zero
  VecZeroEntries(rightHandSide);
  
  // node numbers: 
  // 6 7  (top)
  // 4 5
  
  // 2 3  (bottom)
  // 0 1  
  
  // setup arrays used for integration
  std::array<double, Integrator3D::samplingArraySize()> samplingPoints = Integrator3D::samplingPoints();
  std::array<double, Integrator3D::numberEvaluations()> evaluations[nDofsPerElement][nDofsPerElement];
  
  LOG(DEBUG) << "1D integration with " << IntegratorType::numberEvaluations() << " evaluations";
  LOG(DEBUG) << "3D integration with " << Integrator3D::numberEvaluations() << " evaluations";
#ifdef DEBUG
  LOG(DEBUG) << "SAMPLING POINTS: ";
  for  (auto value : samplingPoints)
    LOG(DEBUG) << "   " << value;
#endif
  
  // loop over elements 
  for (int elementNo = 0; elementNo < mesh->nElements(); elementNo++)
  {
    auto dof = Basis3D::getElementDofs(elementNo, nElementsArray);
    
    Vec3 node[nDofsPerElement];
    for (int dofIndex = 0; dofIndex < nDofsPerElement; dofIndex++)
    {
      // get the global dof no. from element-local dofIndex
      node[dofIndex] = mesh->getNodePosition(dof[dofIndex]);
    }
    
    // compute integral
    for (unsigned int samplingPointIndex=0; samplingPointIndex<samplingPoints.size()/3; samplingPointIndex++)
    {
      // evaluate function to integrate at samplingPoints[i*3], write value to evaluations[i]
      double xi1 = samplingPoints[samplingPointIndex*3+0];
      double xi2 = samplingPoints[samplingPointIndex*3+1];
      double xi3 = samplingPoints[samplingPointIndex*3+2];
      
      // compute the 3x3 jacobian of the parameter space to world space mapping
      auto jacobian = Basis3D::computeJacobian(node, std::array<double,3>({xi1,xi2,xi3}));
      Vec3 jacobianColumn0 = jacobian[0];
      Vec3 jacobianColumn1 = jacobian[1];
      Vec3 jacobianColumn2 = jacobian[2];
        
      // compute determinant of the jacobian
      double determinant = MathUtility::computeDeterminant(jacobianColumn0, jacobianColumn1, jacobianColumn2);
      
      for (int i=0; i<nDofsPerElement; i++)
      {
        for (int j=0; j<nDofsPerElement; j++)
        {
          double integrand = Basis3D::phi(i,std::array<double,3>({xi1,xi2,xi3})) * Basis3D::phi(j,std::array<double,3>({xi1,xi2,xi3})) * fabs(determinant);
          evaluations[i][j][samplingPointIndex] = integrand;
        }
      }
    }  // function evaluations
    
    // perform integration and add to entry of stiffness matrix
    for (int i=0; i<nDofsPerElement; i++)
    {
      for (int j=0; j<nDofsPerElement; j++)
      {
        double value = Integrator3D::integrate(evaluations[i][j]);
        ierr = VecSetValue(rightHandSide, dof[i], value*rhsValues[dof[j]], ADD_VALUES); CHKERRV(ierr);
      }  // j
    }  // i
  }  // elementNo
  
  VecAssemblyBegin(rightHandSide);
  VecAssemblyEnd(rightHandSide);
}

// 1D discretization matrix
template<typename BasisFunctionType, typename IntegratorType>
void FiniteElementMethodBaseTimeStepping<Mesh::Deformable<1>, BasisFunctionType, IntegratorType>::
createRhsDiscretizationMatrix()
{
  // check if matrix discretization matrix exists
  if (!this->data_.discretizationMatrixInitialized())
  {
    this->data_.initializeDiscretizationMatrix();
    
    // set entries of matrix
    LOG(TRACE)<<"createRhsDiscretizationMatrix 1D";
    
    // dmatrix * f_strong = rhs_weak
    // row of dmatrix: contributions to a single entry in rhs_weak
    
    std::shared_ptr<Mesh::Deformable<1>> mesh = std::static_pointer_cast<Mesh::Deformable<1>>(this->data_.mesh());
    
    // get settings values
    int nElements = mesh->nElements();
    PetscErrorCode ierr;
    Mat &dmatrix = this->data_.discretizationMatrix();
    
    typedef BasisFunction::TensorProduct<1,BasisFunctionType> Basis1D;
    
    int nDofsPerElement = Basis1D::nDofsPerElement();
    std::array<int,1> nElementsArray{nElements};
  
    // node numbers: 
    // 0 1
      
    // initialize values to zero
    // loop over elements 
    for (element_idx_t elementNo = 0; elementNo < mesh->nElements(); elementNo++)
    {
      //std::array<int,nDofsPerElement>
      auto dof = Basis1D::getElementDofs(elementNo, nElementsArray);
      // for linear Lagrange basis, the previous loop has the same effect as the following lines:
      //int dof0 = elementNo;  
      //int dof1 = elementNo+1;
      
      for (int i=0; i<nDofsPerElement; i++)
      {
        for (int j=0; j<nDofsPerElement; j++)
        {
          ierr = MatSetValue(dmatrix, dof[i], dof[j], 0, INSERT_VALUES); CHKERRV(ierr);
        }
      }
    }

    // setup arrays used for integration
    std::array<double, IntegratorType::numberEvaluations()> samplingPoints = IntegratorType::samplingPoints();
    std::array<double, IntegratorType::numberEvaluations()> evaluations[nDofsPerElement][nDofsPerElement];
    
    // loop over elements 
    for (element_idx_t elementNo = 0; elementNo < mesh->nElements(); elementNo++)
    {
      auto dof = Basis1D::getElementDofs(elementNo, nElementsArray);
      
      Vec3 node[nDofsPerElement];
      for (int dofIndex = 0; dofIndex < nDofsPerElement; dofIndex++)
      {
        // get the global dof no. from element-local dofIndex
        node[dofIndex] = mesh->getNodePosition(dof[dofIndex]);
      }
      
      double elementLength = MathUtility::distance(node[0], node[nDofsPerElement-1]);
      double integralFactor = elementLength;
      
      // compute integral
      // loop over sampling points (e.g. Gauss points, depending on quadrature scheme)
      for (unsigned int samplingPointIndex=0; samplingPointIndex<samplingPoints.size(); samplingPointIndex++)
      {
        // evaluate function to integrate at samplingPoints[i*D], write value to evaluations[i]
        double xi = samplingPoints[samplingPointIndex];
        
        // loop over pairs of dofs of the current element
        for (int i=0; i<nDofsPerElement; i++)
        {
          for (int j=0; j<nDofsPerElement; j++)
          {
            double integrand = BasisFunctionType::phi(i,xi) * BasisFunctionType::phi(j,xi) * integralFactor;
            evaluations[i][j][samplingPointIndex] = integrand;
          }
        }
      }
      
      // perform quadrature from stored function evaluations and store result to matrix
      // loop over pairs of dofs of the current element
      for (int i=0; i<nDofsPerElement; i++)
      {
        for (int j=0; j<nDofsPerElement; j++)
        {
          double value = IntegratorType::integrate(evaluations[i][j]);
          ierr = MatSetValue(dmatrix, dof[i], dof[j], value, ADD_VALUES); CHKERRV(ierr);
        }
      }
    }
    
    ierr = MatAssemblyBegin(dmatrix, MAT_FINAL_ASSEMBLY); CHKERRV(ierr);
    ierr = MatAssemblyEnd(dmatrix, MAT_FINAL_ASSEMBLY); CHKERRV(ierr);
  }
}

// 2D discretization matrix
template<typename BasisFunctionType, typename IntegratorType>
void FiniteElementMethodBaseTimeStepping<Mesh::Deformable<2>, BasisFunctionType, IntegratorType>::
createRhsDiscretizationMatrix()
{
  // check if matrix discretization matrix exists
  if (!this->data_.discretizationMatrixInitialized())
  {
    this->data_.initializeDiscretizationMatrix();
    
    // set entries of matrix
    LOG(DEBUG)<<"createRhsDiscretizationMatrix 2D";
        
    // dmatrix * f_strong = rhs_weak
    // row of dmatrix: contributions to a single entry in rhs_weak
    
    std::shared_ptr<Mesh::Deformable<2>> mesh = std::static_pointer_cast<Mesh::Deformable<2>>(this->data_.mesh());
    const int nElements0 = mesh->nElements(0);
    const int nElements1 = mesh->nElements(1);
    
    PetscErrorCode ierr;
    Mat &dmatrix = this->data_.discretizationMatrix();
      
    typedef Integrator::TensorProduct<2,IntegratorType> Integrator2D;
    typedef BasisFunction::TensorProduct<2,BasisFunctionType> Basis2D;
      
    int nDofsPerElement = Basis2D::nDofsPerElement();
    std::array<int,2> nElementsArray{nElements0, nElements1};
      
    // node numbers for linear Lagrange basis:
    // 2 3
    // 0 1 
    
    // initialize values to zero
    // loop over elements 
    for (int elementNo = 0; elementNo < mesh->nElements(); elementNo++)
    {
      auto dof = Basis2D::getElementDofs(elementNo, nElementsArray);
      
      for (int i=0; i<nDofsPerElement; i++)
      {
        for (int j=0; j<nDofsPerElement; j++)
        {
          ierr = MatSetValue(dmatrix, dof[i], dof[j], 0, INSERT_VALUES); CHKERRV(ierr);
        }
      }
    }
      
    // setup arrays used for integration
    std::array<double, Integrator2D::samplingArraySize()> samplingPoints = Integrator2D::samplingPoints();
    std::array<double, Integrator2D::numberEvaluations()> evaluations[nDofsPerElement][nDofsPerElement];
    
    // loop over elements 
    for (int elementNo = 0; elementNo < mesh->nElements(); elementNo++)
    {
      // get indices of element-local dofs
      auto dof = Basis2D::getElementDofs(elementNo, nElementsArray);
      
      // get node positions for dofs
      std::array<Vec3,Basis2D::nDofsPerElement()> node;
      for (int dofIndex = 0; dofIndex < nDofsPerElement; dofIndex++)
      {
        node[dofIndex] = mesh->getNodePosition(dof[dofIndex]);
      }
      
      // compute integral
      for (unsigned int samplingPointIndex=0; samplingPointIndex<samplingPoints.size()/2; samplingPointIndex++)
      {
        // evaluate function to integrate at samplingPoints[i*2], write value to evaluations[i]
        double xi1 = samplingPoints[samplingPointIndex*2+0];
        double xi2 = samplingPoints[samplingPointIndex*2+1];
        
        // compute the 3x2 jacobian of the parameter space to world space mapping
        //Vec3 jacobianColumn0 = (1-xi2) * (node[1]-node[0]) + xi2 * (node[3]-node[2]);
        //Vec3 jacobianColumn1 = (1-xi2) * (node[2]-node[0]) + xi2 * (node[3]-node[1]);
        auto jacobian = Basis2D::computeJacobian(node, std::array<double,2>({xi1,xi2}));
        Vec3 jacobianColumn0 = jacobian[0];
        Vec3 jacobianColumn1 = jacobian[1];
        
        double integrationFactor = MathUtility::compute2DIntegrationFactor(jacobianColumn0, jacobianColumn1);
            
        for (int i=0; i<nDofsPerElement; i++)
        {
          for (int j=0; j<nDofsPerElement; j++)
          {
            double integrand = Basis2D::phi(i,std::array<double,2>({xi1,xi2})) * Basis2D::phi(j,std::array<double,2>({xi1,xi2}))
              * integrationFactor;
            evaluations[i][j][samplingPointIndex] = integrand;
          }
        }
      }  // function evaluations
      
      // perform integration and add to entry of stiffness matrix
      for (int i=0; i<nDofsPerElement; i++)
      {
        for (int j=0; j<nDofsPerElement; j++)
        {
          double value = Integrator2D::integrate(evaluations[i][j]);
          ierr = MatSetValue(dmatrix, dof[i], dof[j], value, ADD_VALUES); CHKERRV(ierr);
        }  // j
      }  // i
    }  // elementNo
    
    ierr = MatAssemblyBegin(dmatrix, MAT_FINAL_ASSEMBLY); CHKERRV(ierr);
    ierr = MatAssemblyEnd(dmatrix, MAT_FINAL_ASSEMBLY); CHKERRV(ierr);
  }
}
  
// 3D discretization matrix
template<typename BasisFunctionType, typename IntegratorType>
void FiniteElementMethodBaseTimeStepping<Mesh::Deformable<3>, BasisFunctionType, IntegratorType>::
createRhsDiscretizationMatrix()
{
  // check if matrix discretization matrix exists
  if (!this->data_.discretizationMatrixInitialized())
  {
    this->data_.initializeDiscretizationMatrix();
    
    // set entries of matrix
    LOG(DEBUG)<<"createRhsDiscretizationMatrix 3D";
      
    // get settings values
    std::shared_ptr<Mesh::Deformable<3>> mesh = std::static_pointer_cast<Mesh::Deformable<3>>(this->data_.mesh());
    const int nElements0 = mesh->nElements(0);
    const int nElements1 = mesh->nElements(1);
    const int nElements2 = mesh->nElements(2);
    
    PetscErrorCode ierr;
    Mat &dmatrix = this->data_.discretizationMatrix();
        
    typedef Integrator::TensorProduct<3,IntegratorType> Integrator3D;
    typedef BasisFunction::TensorProduct<3,BasisFunctionType> Basis3D;
      
    int nDofsPerElement = Basis3D::nDofsPerElement();
    std::array<int,3> nElementsArray{nElements0, nElements1, nElements2};
      
    // dmatrix * f_strong = rhs_weak
    // row of dmatrix: contributions to a single entry in rhs_weak
    
    // node numbers: 
    // 6 7  (top)
    // 4 5
    
    // 2 3  (bottom)
    // 0 1 
      
    // initialize values to zero
    // loop over elements 
    
    for (int elementNo = 0; elementNo < mesh->nElements(); elementNo++)
    {
      auto dof = Basis3D::getElementDofs(elementNo, nElementsArray);
    
      for (int i=0; i<nDofsPerElement; i++)
      {
        for (int j=0; j<nDofsPerElement; j++)
        {
          ierr = MatSetValue(dmatrix, dof[i], dof[j], 0, INSERT_VALUES); CHKERRV(ierr);
        }
      }
    }
    
    // setup arrays used for integration
    std::array<double, Integrator3D::samplingArraySize()> samplingPoints = Integrator3D::samplingPoints();
    std::array<double, Integrator3D::numberEvaluations()> evaluations[nDofsPerElement][nDofsPerElement];
    
    LOG(DEBUG) << "1D integration with " << IntegratorType::numberEvaluations() << " evaluations";
    LOG(DEBUG) << "3D integration with " << Integrator3D::numberEvaluations() << " evaluations";
#ifdef DEBUG
    LOG(DEBUG) << "SAMPLING POINTS: ";
    for  (auto value : samplingPoints)
      LOG(DEBUG) << "   " << value;
#endif
        
    // loop over elements 
    for (int elementNo = 0; elementNo < mesh->nElements(); elementNo++)
    {
      // get global dof index of element
      auto dof = Basis3D::getElementDofs(elementNo, nElementsArray);
      
      // get node position of dofs
      std::array<Vec3,Basis3D::nDofsPerElement()> node;
      for (int dofIndex = 0; dofIndex < nDofsPerElement; dofIndex++)
      {
        node[dofIndex] = mesh->getNodePosition(dof[dofIndex]);
      }
      
      // compute integral
      for (unsigned int samplingPointIndex=0; samplingPointIndex<samplingPoints.size()/3; samplingPointIndex++)
      {
        // evaluate function to integrate at samplingPoints[i*3], write value to evaluations[i]
        double xi1 = samplingPoints[samplingPointIndex*3+0];
        double xi2 = samplingPoints[samplingPointIndex*3+1];
        double xi3 = samplingPoints[samplingPointIndex*3+2];
        
        // compute the 3x3 jacobian of the parameter space to world space mapping
        auto jacobian = Basis3D::computeJacobian(node, std::array<double,3>({xi1,xi2,xi3}));
        Vec3 jacobianColumn0 = jacobian[0];
        Vec3 jacobianColumn1 = jacobian[1];
        Vec3 jacobianColumn2 = jacobian[2];
          
        // compute determinant of the jacobian
        double determinant = MathUtility::computeDeterminant(jacobianColumn0, jacobianColumn1, jacobianColumn2);
          
        for (int i=0; i<nDofsPerElement; i++)
        {
          for (int j=0; j<nDofsPerElement; j++)
          {
            double integrand = Basis3D::phi(i,std::array<double,3>({xi1,xi2,xi3})) * Basis3D::phi(j,std::array<double,3>({xi1,xi2,xi3})) * fabs(determinant);
            evaluations[i][j][samplingPointIndex] = integrand;
          }
        }
      }  // function evaluations
      
      // perform integration and add to entry of stiffness matrix
      for (int i=0; i<nDofsPerElement; i++)
      {
        for (int j=0; j<nDofsPerElement; j++)
        {
          double value = Integrator3D::integrate(evaluations[i][j]);
          ierr = MatSetValue(dmatrix, dof[i], dof[j], value, ADD_VALUES); CHKERRV(ierr);
        }  // j
      }  // i
    }  // elementNo
    
    ierr = MatAssemblyBegin(dmatrix, MAT_FINAL_ASSEMBLY); CHKERRV(ierr);
    ierr = MatAssemblyEnd(dmatrix, MAT_FINAL_ASSEMBLY); CHKERRV(ierr);
  }
}
 

};    // namespace