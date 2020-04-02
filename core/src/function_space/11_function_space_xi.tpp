#include "function_space/11_function_space_xi.h"

#include <Python.h>  // has to be the first included header
#include <array>
#include "easylogging++.h"
#include "function_space/00_function_space_base_dim.h"
#include "utility/math_utility.h"
#include <chrono>

namespace FunctionSpace
{

//const double POINT_IN_ELEMENT_EPSILON = 1e-4;    // (1e-5 is too small)
const int N_NEWTON_ITERATIONS = 16;    // (7) (4) (5) number of newton iterations to find out if a point is inside an element
const double RESIDUUM_NORM_TOLERANCE = 1e-4;  // 1e-4^2 is 1e-8, usually it takes 2-3 iterations to reach 1e-2
 
// general implementation
template<typename MeshType,typename BasisFunctionType,typename DummyForTraits>
bool FunctionSpacePointInElement<MeshType,BasisFunctionType,DummyForTraits>::
pointIsInElement(Vec3 point, element_no_t elementNo, std::array<double,MeshType::dim()> &xi, double &residual, double xiTolerance)
{
  // This method computes the xi coordinates in the element-local coordinate system [0,1]^D of the point p and 
  // then checks if the point is inside the element with given xiTolerance (then returns true).
  // This is accomplished by a simple inversion of the mapping Phi(xi) = point using a Newton scheme.
  // However, there are a lot of tweaks to ensure convergence also for almost singular mappings, as well as ensuring good performance.
  // The algorithm is basically two steps: First, do the normal Newton scheme for maximum N_NEWTON_ITERATIONS iterations until the inversion converged by RESIDUUM_NORM_TOLERANCE.
  // In 99% of all cases this is sufficient. In some cases the inversion is still not solved. For those cases, repeat the Newton solve
  // with special reset operations, at the end use the best xi value that was found during the search, not the one of the last iteration.
  
  // timing measurements are disabled, they showed that 'computeApproximateXiForPoint' makes sense and is faster than just initializing the initial guess to 0
#if 0 
  static double durationApproximation = 0.0;
  static double durationNewton = 0.0;
  static int nMeasurements = 0;
 
  auto tStart = std::chrono::steady_clock::now();
#endif

  // define constants
  const int D = MeshType::dim();
  const int nDofsPerElement = FunctionSpaceFunction<MeshType,BasisFunctionType>::nDofsPerElement();
  
  VLOG(2) << "pointIsInElement(" << point << " element " << elementNo << ")";

  // for 3D mesh and linear Lagrange basis function compute approximate xi by heuristic, else set to 0.5
  this->computeApproximateXiForPoint(point, elementNo, xi);
   
#if 0 
  auto tEnd = std::chrono::steady_clock::now();
  durationApproximation += std::chrono::duration_cast<std::chrono::duration<double> >(tEnd - tStart).count();
  nMeasurements++;
  
  tStart = std::chrono::steady_clock::now();
#endif

  // get geometry field (which are the node positions for Lagrange basis and node positions and derivatives for Hermite)
  std::array<Vec3,nDofsPerElement> geometryValues;
  this->getElementGeometry(elementNo, geometryValues);

  VLOG(2) << "point " << point << ", geometryValues: " << geometryValues;

  if (this->pointIsOutsideBoundingBox(point, geometryValues))
  {
    VLOG(2) << "point " << point << " is outside bounding box";
    return false;
  }

  if (this->pointIsNodePosition(point, geometryValues, xi))
  {
    VLOG(2) << "point " << point << " is a node position of the element";
    return true;
  }

  std::array<double,MeshType::dim()> xiPrevious = xi;

  // compute initial residuum
  Vec3 residuum = point - this->template interpolateValueInElement<3>(geometryValues, xi);
  double residuumNormSquared = MathUtility::normSquared<3>(residuum);
  double residuumNormSquaredPrevious = residuumNormSquared;
  
  if (VLOG_IS_ON(2))
  {
    VLOG(2) << " xi0 = " << xi << ", residuum: " << residuum << " (norm: " << sqrt(residuumNormSquared) << ")";
  }

  // initialize the increment for xi if the Newton step fails
  std::array<std::array<double,MeshType::dim()>,6> xiStep{};
  VLOG(2) << " xiStep uninit: " << xiStep;
  if (MeshType::dim() == 3)
  {
    xiStep[0][0] = 0.01;
    xiStep[1][1] = 0.01;
    xiStep[2][2] = 0.01;
    xiStep[3][0] = -0.01;
    xiStep[4][1] = -0.01;
    xiStep[5][2] = -0.01;
  }
  else if (MeshType::dim() == 2)
  {
    xiStep[0][0] = 0.01;
    xiStep[1][1] = 0.01;
    xiStep[2][0] = -0.01;
    xiStep[3][1] = -0.01;
    xiStep[4][0] = 0.01;
    xiStep[4][1] = 0.01;
    xiStep[5][0] = -0.01;
    xiStep[5][1] = -0.01;
  }
  else if (MeshType::dim() == 1)
  {
    xiStep[0][0] = 0.01;
    xiStep[1][0] = -0.01;
    xiStep[2][0] = 0.01;
    xiStep[3][0] = -0.01;
    xiStep[4][0] = 0.01;
    xiStep[5][0] = -0.01;
  }
  VLOG(2) << " xiStep: " << xiStep;
  int nIterations = 0;

  // while the residuum norm is above the tolerance
  for (int iterationNo = 0; iterationNo < N_NEWTON_ITERATIONS && residuumNormSquared > MathUtility::sqr(RESIDUUM_NORM_TOLERANCE); iterationNo++, nIterations++)
  { 
    // perform Newton step
    // Phi(xi) = point
    // Phi(xi) = Phi(xi0) + J*(xi-xi0)  => xi = xi0 + Jinv*(point - Phi(xi0))  
    Tensor2<D> inverseJacobian = this->getInverseJacobian(geometryValues, elementNo, xi);
    xi += inverseJacobian * MathUtility::transformToD<D,3>(residuum);
    
    // compute residuum
    residuum = point - this->template interpolateValueInElement<3>(geometryValues, xi);
    residuumNormSquared = MathUtility::normSquared<3>(residuum);
    
    // if the residuum value jumped more than 1000 in one step, discard current step and restart with a slightly different xi value
    if (residuumNormSquared - residuumNormSquaredPrevious > 1000)
    {
      VLOG(2) << "jump in norm " << residuumNormSquaredPrevious << " -> " << residuumNormSquared << ", xi: " << xi;

      xi = xiPrevious + xiStep[iterationNo%6];
      residuum = point - this->template interpolateValueInElement<3>(geometryValues, xi);
      residuumNormSquared = MathUtility::normSquared<3>(residuum);

      VLOG(2) << "reset to xi=" << xi << ", norm: " << residuumNormSquared;
    }
    
    residuumNormSquaredPrevious = residuumNormSquared;
    xiPrevious = xi;

    if (VLOG_IS_ON(3))   // extra if because of sqrt
    {
      VLOG(3) << " xi_" << iterationNo << " = "  << xi << ", residuum: " << residuum << " (norm: " << sqrt(residuumNormSquared) << ")";
    }
  }
  residual = residuumNormSquared;


  // check if point is inside the element by looking at the value of xi
  double epsilon = xiTolerance;
  bool pointIsInElement = true;
  for (int i = 0; i < D; i++)
  {
    if (!(0.0-epsilon <= xi[i] && xi[i] <= 1.0+epsilon))
    {
      pointIsInElement = false;
    }
  }

  // if the norm is not good, i.e. the solution was not found in the previous loop
  if (residuumNormSquared > MathUtility::sqr(RESIDUUM_NORM_TOLERANCE))
  {
    // if current value for residuum is bad, restart completely
    if (residuumNormSquared > 1)
    {
      // for 3D mesh and linear Lagrange basis function compute approximate xi by heuristic, else set to 0.5
      this->computeApproximateXiForPoint(point, elementNo, xi);
          
      // compute initial residuum
      residuum = point - this->template interpolateValueInElement<3>(geometryValues, xi);
      residuumNormSquared = MathUtility::normSquared<3>(residuum);
      residuumNormSquaredPrevious = residuumNormSquared;
    }
    
    // iterate further and note the best found xi on the way
    double initialResidual = residuumNormSquared;
    
    // variables for the best value of xi with the lowest residual norm that was found
    std::array<double,MeshType::dim()> bestXi = xi;
    double bestResidual = residuumNormSquared;
    VLOG(2) << "point not found yet, restart with Xi=" << xi << ", bestResidual: " << bestResidual;


    // Again, do a lot of iterations to get closer to the correct xi value. This occurs rarely.
    // Note, increasing thet number of iterations further has no significant effect. For the cases where this loop is required, the Newton scheme kind of fails,
    // the residual drops step by step and then in one single step increases sharply.
    for (int iterationNo = 0; iterationNo < 2*N_NEWTON_ITERATIONS && residuumNormSquared > MathUtility::sqr(RESIDUUM_NORM_TOLERANCE); iterationNo++, nIterations++)
    {
      // perform Newton step
      Tensor2<D> inverseJacobian = this->getInverseJacobian(geometryValues, elementNo, xi);
      xi += inverseJacobian * MathUtility::transformToD<D,3>(residuum);

      // compute residuum
      residuum = point - this->template interpolateValueInElement<3>(geometryValues, xi);
      residuumNormSquared = MathUtility::normSquared<3>(residuum);
      
      // if the residuum value jumped more than 1000 in one step, discard current step and restart with a slightly different xi value
      if (residuumNormSquared - residuumNormSquaredPrevious > 1000)
      {
        VLOG(2) << "jump in norm " << residuumNormSquaredPrevious << " -> " << residuumNormSquared << ", xi: " << xi;

        xi = xiPrevious + xiStep[iterationNo%6];
        residuum = point - this->template interpolateValueInElement<3>(geometryValues, xi);
        residuumNormSquared = MathUtility::normSquared<3>(residuum);

        VLOG(2) << "reset to xi=" << xi << ", norm: " << residuumNormSquared;
      }

      if (residuumNormSquared < bestResidual)
      {
        bestResidual = residuumNormSquared;
        bestXi = xi;
      }

      residuumNormSquaredPrevious = residuumNormSquared;
      xiPrevious = xi;

      if (VLOG_IS_ON(3))
      {
        VLOG(3) << " xi_" << iterationNo << " = "  << xi << ", residuum: " << residuum << " (norm: " << sqrt(residuumNormSquared) << ")";
      }
    }

    // Use the best values that were found during the iterations. This is usually not the value from the last iteration.
    xi = bestXi;
    residual = bestResidual;
    residuumNormSquared = bestResidual;

    VLOG(2) << "pointIsInElement point: " << point << ", elementNo: " << elementNo
      << ", residual improved from " << initialResidual << " to " << residual << ", xi: " << xi;

    // check if point is inside the element by looking at the value of xi
    pointIsInElement = true;
    for (int i = 0; i < D; i++)
    {
      if (!(0.0-epsilon <= xi[i] && xi[i] <= 1.0+epsilon))
      {
        pointIsInElement = false;
      }
    }
  }
    
  // if the correct value was still not found, emit a warning, but only in debug mode
#ifndef NDEBUG
  if (pointIsInElement && residuumNormSquared > MathUtility::sqr(RESIDUUM_NORM_TOLERANCE))
  {
    LOG(WARNING) << "pointIsInElement failed after " << 3*N_NEWTON_ITERATIONS << " iterations, point: " 
      << point << ", elementNo: " << elementNo << ", geometryValues: " << geometryValues << ", found xi: " 
      << xi << ", but residual: " << sqrt(residuumNormSquared) << " > " << RESIDUUM_NORM_TOLERANCE 
      << ", (" << residuumNormSquared << " > " << MathUtility::sqr(RESIDUUM_NORM_TOLERANCE) << ")";
  }
#endif

  // measure number of iterations
#if 0
  static std::vector<int> nIterationsList;
  nIterationsList.push_back(nIterations);

  static int counter = 0;
  counter++;

  if (counter % 100000 == 0)
  {
    // compute mean and max
    long long int sum = 0;
    int maximum = 0;
    for (int value : nIterationsList)
    {
      sum += value;
      maximum = std::max(maximum, value);
    }
    LOG(INFO) << "nIterations mean: " << (sum/nIterationsList.size()) << ", max: " << maximum;
    // nIterations mean: 5, max: 96 for N_NEWTON_ITERATIONS = 32, RESIDUUM_NORM_TOLERANCE = 1e-4
  }
#endif

  // measure runtime
#if 0  
  tEnd = std::chrono::steady_clock::now();
  durationNewton += std::chrono::duration_cast<std::chrono::duration<double> >(tEnd - tStart).count();
  
  if (nMeasurements%10000 == 0)
  {
    LOG(INFO) << "duration initial approximation: " << durationApproximation / nMeasurements;
    LOG(INFO) << "duration Newton: " << durationNewton / nMeasurements;
    LOG(INFO) << "total: " << (durationNewton+durationApproximation) / nMeasurements;
  }
#endif  
  VLOG(2) << "  " << (pointIsInElement? "inside" : "outside");
  
  return pointIsInElement;
}

// regular fixed 1D
template<typename BasisFunctionType>
bool FunctionSpacePointInElement<Mesh::StructuredRegularFixedOfDimension<1>, BasisFunctionType>::
pointIsInElement(Vec3 point, element_no_t elementNo, std::array<double,1> &xi, double &residual, double xiTolerance)
{
  const int nDofsPerElement = FunctionSpaceBaseDim<1,BasisFunctionType>::nDofsPerElement();  //=2
  std::array<Vec3, nDofsPerElement> geometryValues;
  
  this->getElementGeometry(elementNo, geometryValues);
  const double elementLength = this->meshWidth() * FunctionSpaceBaseDim<1,BasisFunctionType>::averageNNodesPerElement();
  
  xi[0] = MathUtility::norm<3>(point - geometryValues[0]) / elementLength;
  
  return -xiTolerance <= xi[0] && xi[0] <= 1.0+xiTolerance;
}

// regular fixed 2D
template<typename BasisFunctionType>
bool FunctionSpacePointInElement<Mesh::StructuredRegularFixedOfDimension<2>, BasisFunctionType>::
pointIsInElement(Vec3 point, element_no_t elementNo, std::array<double,2> &xi, double &residual, double xiTolerance)
{
  const int nDofsPerElement = FunctionSpaceBaseDim<2,BasisFunctionType>::nDofsPerElement();
  std::array<Vec3, nDofsPerElement> geometryValues;
  
  this->getElementGeometry(elementNo, geometryValues);
  const double elementLength = this->meshWidth() * FunctionSpaceBaseDim<1,BasisFunctionType>::averageNNodesPerElement();
  
  const double xi1 = (point[0] - geometryValues[0][0]) / elementLength;
  const double xi2 = (point[1] - geometryValues[0][1]) / elementLength;
  
  return (-xiTolerance <= xi1 && xi1 <= 1.0+xiTolerance) && (-xiTolerance <= xi2 && xi2 <= 1.0+xiTolerance);
}

// regular fixed 3D
template<typename BasisFunctionType>
bool FunctionSpacePointInElement<Mesh::StructuredRegularFixedOfDimension<3>, BasisFunctionType>::
pointIsInElement(Vec3 point, element_no_t elementNo, std::array<double,3> &xi, double &residual, double xiTolerance)
{
  const int nDofsPerElement = FunctionSpaceBaseDim<3,BasisFunctionType>::nDofsPerElement();
  std::array<Vec3, nDofsPerElement> geometryValues;
  
  this->getElementGeometry(elementNo, geometryValues);
  const double elementLength = this->meshWidth() * FunctionSpaceBaseDim<1,BasisFunctionType>::averageNNodesPerElement();
  
  const double xi1 = (point[0] - geometryValues[0][0]) / elementLength;
  const double xi2 = (point[1] - geometryValues[0][1]) / elementLength;
  const double xi3 = (point[2] - geometryValues[0][2]) / elementLength;
  
  return (-xiTolerance <= xi1 && xi1 <= 1.0+xiTolerance) && (-xiTolerance <= xi2 && xi2 <= 1.0+xiTolerance) && (-xiTolerance <= xi3 && xi3 <= 1.0+xiTolerance);
}

// 1D deformable meshes and linear shape function
template<typename MeshType>
bool FunctionSpacePointInElement<MeshType, BasisFunction::LagrangeOfOrder<1>, Mesh::isDeformableWithDim<1,MeshType>>::
pointIsInElement(Vec3 point, element_no_t elementNo, std::array<double,1> &xi, double &residual, double xiTolerance)
{
  //const int nDofsPerElement = FunctionSpaceBaseDim<1,BasisFunction::LagrangeOfOrder<1>>::nDofsPerElement();  //=2
  const int nDofsPerElement = 2;
  std::array<Vec3, nDofsPerElement> geometryValues;
  
  // get needed variables
  this->getElementGeometry(elementNo, geometryValues);
  
  const double xp1 = point[0];
  
  const double x11 = geometryValues[0][0];
  const double x21 = geometryValues[1][0];
  
  // compute analytic solution for xi
  const double xi1 = (x11 - xp1)/(x11 - x21);
  xi[0] = xi1;
  
  return -xiTolerance <= xi1 && xi1 <= 1.0+xiTolerance;
}
/*
// 2D deformable meshes and linear shape function
template<typename MeshType>
bool FunctionSpacePointInElement<MeshType, BasisFunction::LagrangeOfOrder<1>, Mesh::isDeformableWithDim<2,MeshType>>::
pointIsInElement(Vec3 point, element_no_t elementNo, std::array<double,2> &xi, double xiTolerance)
{
  //const int nDofsPerElement = FunctionSpaceBaseDim<2,BasisFunction::LagrangeOfOrder<1>>::nDofsPerElement();  //=4
  const int nDofsPerElement = 4;
  std::array<Vec3, nDofsPerElement> geometryValues;
  
  // get needed variables
  this->getElementGeometry(elementNo, geometryValues);
  
  MathUtility::quadrilateralGetPointCoordinates(geometryValues, point, xi);

  VLOG(2) << "pointIsInElement (2D linear) el no. " << elementNo << ", point " << point << ", xiTolerance: " << xiTolerance << " -> xi: " << xi;

  const double xi1 = xi[0];
  const double xi2 = xi[1];

  return (-xiTolerance <= xi1 && xi1 <= 1.0+xiTolerance) && (-xiTolerance <= xi2 && xi2 <= 1.0+xiTolerance);
}*/

// 3D deformable meshes and linear shape function

//! check if the point lies inside the 4-point tetrahedron, if yes, return true and set xi to the value of the point
template<typename MeshType>
bool ComputeXiApproximation<MeshType, BasisFunction::LagrangeOfOrder<1>, Mesh::isDeformableWithDim<3,MeshType>>::
pointIsInTetrahedron(Vec3 point, std::array<Vec3,4> tetrahedron, std::array<bool,3> correctOrientation, std::array<double,3> &xi)
{
  const double xp1 = point[0];
  const double xp2 = point[1];
  const double xp3 = point[2];
  
  const double x11 = tetrahedron[1][0];
  const double x12 = tetrahedron[1][1];
  const double x13 = tetrahedron[1][2];
  
  const double x21 = tetrahedron[2][0];
  const double x22 = tetrahedron[2][1];
  const double x23 = tetrahedron[2][2];
  
  const double x31 = tetrahedron[3][0];
  const double x32 = tetrahedron[3][1];
  const double x33 = tetrahedron[3][2];
  
  const double x41 = tetrahedron[0][0];
  const double x42 = tetrahedron[0][1];
  const double x43 = tetrahedron[0][2];
  
  const double det = (x11 - x41)*(x22 - x42)*(x33 - x43) - (x11 - x41)*(x23 - x43)*(x32 - x42) - (x12 - x42)*(x21 - x41)*(x33 - x43) + (x12 - x42)*(x23 - x43)*(x31 - x41) + (x13 - x43)*(x21 - x41)*(x32 - x42) - (x13 - x43)*(x22 - x42)*(x31 - x41);
  xi[0] = 1./det * ((-x41 + xp1)*((x22 - x42)*(x33 - x43) - (x23 - x43)*(x32 - x42)) + (-x42 + xp2)*(-(x21 - x41)*(x33 - x43) + (x23 - x43)*(x31 - x41)) + (-x43 + xp3)*((x21 - x41)*(x32 - x42) - (x22 - x42)*(x31 - x41)));
  xi[1] = 1./det * ((-x41 + xp1)*(-(x12 - x42)*(x33 - x43) + (x13 - x43)*(x32 - x42)) + (-x42 + xp2)*((x11 - x41)*(x33 - x43) - (x13 - x43)*(x31 - x41)) + (-x43 + xp3)*(-(x11 - x41)*(x32 - x42) + (x12 - x42)*(x31 - x41)));
  xi[2] = 1./det * ((-x41 + xp1)*((x12 - x42)*(x23 - x43) - (x13 - x43)*(x22 - x42)) + (-x42 + xp2)*(-(x11 - x41)*(x23 - x43) + (x13 - x43)*(x21 - x41)) + (-x43 + xp3)*((x11 - x41)*(x22 - x42) - (x12 - x42)*(x21 - x41)));

  const double epsilon = 1e-12;
  bool pointIsInQuarterSpace = (xi[0] >= -epsilon && xi[1] >= -epsilon && xi[2] >= -epsilon);
  
  if (!correctOrientation[0])
    xi[0] = 1. - xi[0];
  if (!correctOrientation[1])
    xi[1] = 1. - xi[1];
  if (!correctOrientation[2])
    xi[2] = 1. - xi[2];
   
  pointIsInQuarterSpace = (xi[0] >= -epsilon && xi[1] >= -epsilon && xi[2] >= -epsilon);
  
  VLOG(3) << "   xi: " << xi << "(in:" << pointIsInQuarterSpace << ")";
  return pointIsInQuarterSpace;
}

// 3D deformable meshes and linear shape function
template<typename MeshType>
bool ComputeXiApproximation<MeshType, BasisFunction::LagrangeOfOrder<1>, Mesh::isDeformableWithDim<3,MeshType>>::
pointIsInElementQuick(Vec3 point, const std::array<Vec3, FunctionSpaceBaseDim<3,BasisFunction::LagrangeOfOrder<1>>::nDofsPerElement()> &geometryValues) const
{
  VLOG(3) << "pointIsInElementQuick, point " << point << ", element " << geometryValues;

  // this is a heuristic to check whether point is inside the element
  // for every quadrilateral it is checked if the pyramid with the quadrilateral as base and the point as top is oriented such that the point is on the correct side of the element
  // Note that this heuristic is wrong sometimes! Therefore it does not get used anywhere.
  
  const double epsilon = 1e-12;
  
  const Vec3 p0 = (-geometryValues[0]+point);
  const Vec3 p1 = (-geometryValues[1]+point);
  const Vec3 p2 = (-geometryValues[2]+point);
  const Vec3 p3 = (-geometryValues[3]+point);
  const Vec3 p4 = (-geometryValues[4]+point);
  const Vec3 p5 = (-geometryValues[5]+point);
  const Vec3 p6 = (-geometryValues[6]+point);
  const Vec3 p7 = (-geometryValues[7]+point);
  
  if (MathUtility::normSquared<3>(p0) < epsilon || MathUtility::normSquared<3>(p1) < epsilon || MathUtility::normSquared<3>(p2) < epsilon || MathUtility::normSquared<3>(p3) < epsilon
    || MathUtility::normSquared<3>(p4) < epsilon || MathUtility::normSquared<3>(p5) < epsilon || MathUtility::normSquared<3>(p6) < epsilon || MathUtility::normSquared<3>(p7) < epsilon)
  {
    return true;
  }
 
  // bottom [geometryValues[0],geometryValues[1],geometryValues[3],geometryValues[2]]
  const Vec3 temp30 = (-geometryValues[3]+geometryValues[0]);
  const Vec3 temp01 = (-geometryValues[0]+geometryValues[1]);
  const Vec3 temp12 = (-geometryValues[1]+geometryValues[2]);
  const Vec3 temp32 = (-geometryValues[3]+geometryValues[2]);
  const Vec3 temp20 = (-geometryValues[2]+geometryValues[0]);
  
  const bool v0 = MathUtility::dot(MathUtility::cross(temp30, temp01), p0) >= -epsilon;
  const bool v1 = MathUtility::dot(MathUtility::cross(temp01, temp12), p1) >= -epsilon;
  const bool v2 = MathUtility::dot(MathUtility::cross(temp12, temp32), p3) >= -epsilon;
  const bool v3 = MathUtility::dot(MathUtility::cross(temp32, temp20), p2) >= -epsilon;
  
  // top [geometryValues[4],geometryValues[6],geometryValues[7],geometryValues[5]]
  const Vec3 temp74 = (-geometryValues[7]+geometryValues[4]);
  const Vec3 temp46 = (-geometryValues[4]+geometryValues[6]);
  const Vec3 temp65 = (-geometryValues[6]+geometryValues[5]);
  const Vec3 temp75 = (-geometryValues[7]+geometryValues[5]);
  const Vec3 temp54 = (-geometryValues[5]+geometryValues[4]);
  
  const bool v4 = MathUtility::dot(MathUtility::cross(temp74, temp46), p4) >= -epsilon;
  const bool v5 = MathUtility::dot(MathUtility::cross(temp46, temp65), p6) >= -epsilon;
  const bool v6 = MathUtility::dot(MathUtility::cross(temp65, temp75), p7) >= -epsilon;
  const bool v7 = MathUtility::dot(MathUtility::cross(temp75, temp54), p5) >= -epsilon;
  
  // right [geometryValues[1],geometryValues[5],geometryValues[7],geometryValues[3]]
  const Vec3 temp71 = (-geometryValues[7]+geometryValues[1]);
  const Vec3 temp15 = (-geometryValues[1]+geometryValues[5]);
  const Vec3 temp53 = (-geometryValues[5]+geometryValues[3]);
  const Vec3 temp73 = (-geometryValues[7]+geometryValues[3]);
  const Vec3 temp31 = (-geometryValues[3]+geometryValues[1]);
  
  const bool v8 = MathUtility::dot(MathUtility::cross(temp71, temp15), p1) >= -epsilon;
  const bool v9 = MathUtility::dot(MathUtility::cross(temp15, temp53), p5) >= -epsilon;
  const bool v10 = MathUtility::dot(MathUtility::cross(temp53, temp73), p7) >= -epsilon;
  const bool v11 = MathUtility::dot(MathUtility::cross(temp73, temp31), p3) >= -epsilon;
  
  // left [geometryValues[0],geometryValues[2],geometryValues[6],geometryValues[4]]
  const Vec3 temp60 = (-geometryValues[6]+geometryValues[0]);
  const Vec3 temp02 = (-geometryValues[0]+geometryValues[2]);
  const Vec3 temp24 = (-geometryValues[2]+geometryValues[4]);
  const Vec3 temp64 = (-geometryValues[6]+geometryValues[4]);
  const Vec3 temp40 = (-geometryValues[4]+geometryValues[0]);
  
  const bool v12 = MathUtility::dot(MathUtility::cross(temp60, temp02), p0) >= -epsilon;
  const bool v13 = MathUtility::dot(MathUtility::cross(temp02, temp24), p2) >= -epsilon;
  const bool v14 = MathUtility::dot(MathUtility::cross(temp24, temp64), p6) >= -epsilon;
  const bool v15 = MathUtility::dot(MathUtility::cross(temp64, temp40), p4) >= -epsilon;
  
  // front [geometryValues[0],geometryValues[4],geometryValues[5],geometryValues[1]]
  const Vec3 temp50 = (-geometryValues[5]+geometryValues[0]);
  const Vec3 temp04 = (-geometryValues[0]+geometryValues[4]);
  const Vec3 temp41 = (-geometryValues[4]+geometryValues[1]);
  const Vec3 temp51 = (-geometryValues[5]+geometryValues[1]);
  const Vec3 temp10 = (-geometryValues[1]+geometryValues[0]);
  
  const bool v16 = MathUtility::dot(MathUtility::cross(temp50, temp04), p0) >= -epsilon;
  const bool v17 = MathUtility::dot(MathUtility::cross(temp04, temp41), p4) >= -epsilon;
  const bool v18 = MathUtility::dot(MathUtility::cross(temp41, temp51), p5) >= -epsilon;
  const bool v19 = MathUtility::dot(MathUtility::cross(temp51, temp10), p1) >= -epsilon;
  
  // bconst Vec3 tempck [geometryValues[2],geometryValues[3],geometryValues[7],geometryValues[6]]
  const Vec3 temp72 = (-geometryValues[7]+geometryValues[2]);
  const Vec3 temp23 = (-geometryValues[2]+geometryValues[3]);
  const Vec3 temp36 = (-geometryValues[3]+geometryValues[6]);
  const Vec3 temp76 = (-geometryValues[7]+geometryValues[6]);
  const Vec3 temp62 = (-geometryValues[6]+geometryValues[2]);
  
  const bool v20 = MathUtility::dot(MathUtility::cross(temp72, temp23), p2) >= -epsilon;
  const bool v21 = MathUtility::dot(MathUtility::cross(temp23, temp36), p3) >= -epsilon;
  const bool v22 = MathUtility::dot(MathUtility::cross(temp36, temp76), p7) >= -epsilon;
  const bool v23 = MathUtility::dot(MathUtility::cross(temp76, temp62), p6) >= -epsilon;
  
  
  bool isInside = v0 && v1 && v2 && v3 && v4 && v5 && v6 && v7 && v8 && v9 && v10 && v11 && v12 && v13 && v14 
    && v15 && v16 && v17 && v18 && v19 && v20 && v21 && v22 && v23;
    
  VLOG(3) << (isInside? "inside" : "outside");

  if (!isInside)
  {
    VLOG(1) << " isInside: " << v0 << "," <<  v1 << "," <<  v2 << "," <<  v3 << ", " <<  v4 << "," <<  v5 << "," <<  v6 << "," <<  v7 << ", " <<  v8 << "," <<  v9 << "," <<  v10 << "," <<  v11 << ", " <<  v12 << "," <<  v13 << "," <<  v14 
      << "," <<  v15 << ", " <<  v16 << "," <<  v17 << "," <<  v18 << "," <<  v19 << ", " <<  v20 << "," <<  v21 << "," <<  v22 << "," <<  v23;
  }
  else
  {
    VLOG(3) << " isInside: " << v0 << "," <<  v1 << "," <<  v2 << "," <<  v3 << "," <<  v4 << "," <<  v5 << "," <<  v6 << "," <<  v7 << "," <<  v8 << "," <<  v9 << "," <<  v10 << "," <<  v11 << "," <<  v12 << "," <<  v13 << "," <<  v14 
      << "," <<  v15 << "," <<  v16 << "," <<  v17 << "," <<  v18 << "," <<  v19 << "," <<  v20 << "," <<  v21 << "," <<  v22 << "," <<  v23;
  }

  return isInside;
}

template<typename MeshType,typename BasisFunctionType,typename DummyForTraits>
void ComputeXiApproximation<MeshType,BasisFunctionType,DummyForTraits>::
computeApproximateXiForPoint(Vec3 point, element_no_t elementNo, std::array<double,MeshType::dim()> &xi)
{
  xi.fill(0.5);
}

template<typename MeshType>
void ComputeXiApproximation<MeshType, BasisFunction::LagrangeOfOrder<1>, Mesh::isDeformableWithDim<3,MeshType>>::
computeApproximateXiForPoint(Vec3 point, element_no_t elementNo, std::array<double,3> &xi)
{
  // This computes a fast approximation to xi, which can then be refined by a newton scheme
 
  //const int nDofsPerElement = FunctionSpaceBaseDim<3,BasisFunction::LagrangeOfOrder<1>>::nDofsPerElement();  //=8
  std::array<Vec3, 8> geometryValues;
  
  this->getElementGeometry(elementNo, geometryValues);
  
  VLOG(3) << "computeApproximateXiForPoint, point " << point << ", element " << elementNo << geometryValues;
  
  xi.fill(0.0);
  std::array<double,3> xiSum{};
  int nSummands = 0;
  
  // p0
  if (pointIsInTetrahedron(point, std::array<Vec3,4>{geometryValues[0], geometryValues[1], geometryValues[2], geometryValues[4]}, std::array<bool,3>{true,true,true}, xiSum))
  {
    xi += xiSum;
    nSummands++;
  }
  else 
  {
    VLOG(3) << "p0 out";
  }
#if 0
  // p1
  if (pointIsInTetrahedron(point, std::array<Vec3,4>{geometryValues[1], geometryValues[0], geometryValues[5], geometryValues[3]}, std::array<bool,3>{false,true,true}, xiSum))
  {
    std::swap(xi[1],xi[2]);
    xi += xiSum;
    nSummands++;
  }
  else 
  {
    VLOG(3) << "p1 out";
  }
  
  //p2
  if (pointIsInTetrahedron(point, std::array<Vec3,4>{geometryValues[2], geometryValues[3], geometryValues[6], geometryValues[0]}, std::array<bool,3>{true,true,false}, xiSum))
  {
    std::swap(xi[1],xi[2]);
    xi += xiSum;
    nSummands++;
  }
  else 
  {
    VLOG(3) << "p2 out";
  }
#endif    
  // p3
  if (pointIsInTetrahedron(point, std::array<Vec3,4>{geometryValues[3], geometryValues[2], geometryValues[1], geometryValues[7]}, std::array<bool,3>{false,false,true}, xiSum))
  {
    xi += xiSum;
    nSummands++;
  }
  else 
  {
    VLOG(3) << "p3 out";
  }
#if   0
  // p4
  if (pointIsInTetrahedron(point, std::array<Vec3,4>{geometryValues[4], geometryValues[5], geometryValues[0], geometryValues[6]}, std::array<bool,3>{true,false,true}, xiSum))
  {
    std::swap(xi[1],xi[2]);
    xi += xiSum;
    nSummands++;
  }
  else 
  {
    VLOG(3) << "p4 out";
  }
#endif  
  // p5
  if (pointIsInTetrahedron(point, std::array<Vec3,4>{geometryValues[5], geometryValues[4], geometryValues[7], geometryValues[1]}, std::array<bool,3>{false,true,false}, xiSum))
  {
    xi += xiSum;
    nSummands++;
  }
  else 
  {
    VLOG(3) << "p5 out";
  }
  
  // p6
  if (pointIsInTetrahedron(point, std::array<Vec3,4>{geometryValues[6], geometryValues[7], geometryValues[4], geometryValues[2]}, std::array<bool,3>{true,false,false}, xiSum))
  {
    xi += xiSum;
    nSummands++;
  }
  else 
  {
    VLOG(3) << "p6 out";
  }
#if 0  
  // p7
  if (pointIsInTetrahedron(point, std::array<Vec3,4>{geometryValues[7], geometryValues[6], geometryValues[3], geometryValues[5]}, std::array<bool,3>{false,false,false}, xiSum))
  {
    std::swap(xi[1],xi[2]);
    xi += xiSum;
    nSummands++;
  }
  else 
  {
    VLOG(3) << "p7 out";
  }
#endif    

  if (nSummands != 0)
    xi /= nSummands;
}

template<typename MeshType,typename BasisFunctionType>
bool FunctionSpaceXi<MeshType,BasisFunctionType>::
pointIsOutsideBoundingBox(Vec3 point, const std::array<Vec3,FunctionSpaceFunction<MeshType,BasisFunctionType>::nDofsPerElement()> &geometryValues) const
{
  double xmin = geometryValues[0][0];
  double xmax = geometryValues[0][0];
  double ymin = geometryValues[0][1];
  double ymax = geometryValues[0][1];
  double zmin = geometryValues[0][2];
  double zmax = geometryValues[0][2];

  // loop over all geometryValues and find out bounding box values
  for (int i = 1; i < geometryValues.size(); i++)
  {
    xmin = std::min(xmin, geometryValues[i][0]);
    xmax = std::max(xmax, geometryValues[i][0]);
    ymin = std::min(ymin, geometryValues[i][1]);
    ymax = std::max(ymax, geometryValues[i][1]);
    zmin = std::min(zmin, geometryValues[i][2]);
    zmax = std::max(zmax, geometryValues[i][2]);
  }

  // tolerance has to be at least 1e-9 or higher
  const double eps = 1e-9;
  if (point[0] < xmin-eps || point[0] > xmax+eps || point[1] < ymin-eps || point[1] > ymax+eps || point[2] < zmin-eps || point[2] > zmax+eps)
  {
    VLOG(2) << "point " << point << " is outside bounding box [" << xmin << "," << xmax << "]x[" << ymin << "," << ymax << "]x[" << zmin << "," << zmax << "], diff: " 
    "[" << point[0]-xmin << "," << xmax-point[0] << "]x[" << point[1]-ymin << "," << ymax-point[1] << "]x[" << point[2]-zmin << "," << zmax-point[2] << "], eps=" << eps;
    return true;
  }
  return false;
}

template<typename MeshType,typename BasisFunctionType>
bool FunctionSpaceXi<MeshType,BasisFunctionType>::
pointIsNodePosition(Vec3 point, const std::array<Vec3,FunctionSpaceFunction<MeshType,BasisFunctionType>::nDofsPerElement()> &geometryValues, std::array<double,MeshType::dim()> &xi) const
{
  // loop over all geometryValues
  int nNodesPerDimensionX = FunctionSpaceBaseDim<1,BasisFunctionType>::nNodesPerElement();
  int nNodesPerDimensionY = FunctionSpaceBaseDim<1,BasisFunctionType>::nNodesPerElement();
  int nNodesPerDimensionZ = FunctionSpaceBaseDim<1,BasisFunctionType>::nNodesPerElement();

  const int D = MeshType::dim();
  
  if (D <= 2)
    nNodesPerDimensionZ = 1;
  
  if (D == 1)
    nNodesPerDimensionY = 1;

  int index = 0;
  for (int k = 0; k < nNodesPerDimensionZ; k++)
  {
    for (int j = 0; j < nNodesPerDimensionY; j++)
    {
      for (int i = 0; i < nNodesPerDimensionX; i++, index++)
      {
        if (MathUtility::template equals<3>((const Vec3)geometryValues[index], (const Vec3)point, 1e-6))
        {
          xi[0] = (double)i / (nNodesPerDimensionX-1);
          if (D >= 2)
            xi[1] = (double)j / (nNodesPerDimensionY-1);
          if (D == 3)
            xi[2] = (double)k / (nNodesPerDimensionZ-1);
          return true;
        }
        else 
        {
          VLOG(2) << "point " << point << " is not equal to node at " << index << " (" << i << "," << j << "," << k << ") " << geometryValues[index];
        }
      }
    }
  }
  return false;
}

} // namespace
