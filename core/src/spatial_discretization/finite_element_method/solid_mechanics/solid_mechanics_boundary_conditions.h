#pragma once

#include <petscvec.h>
#include <memory>

namespace SpatialDiscretization
{

/** This class stores and handles boundary conditions for the nonlinear solid mechanics problems. There are Dirichlet BC, traction and body force boundary, in reference and current configuration each.
  */
template<typename FunctionSpaceType,typename Term>
class SolidMechanicsBoundaryConditions
{
public:

  //! copy all values that are not constrained by dirichlet BC from the input to the output vector, it is no problem if the output vector is bigger than needed
  void reduceVector(Vec &input, Vec &output, const int nLocalUnknownsInputVector);

  //! reverse operation to reduceVector, adds values of Dirichlet BC
  void expandVector(Vec &input, Vec &output, const int nLocalUnknownsOutputVector);

  //! set entries in displacements to Dirichlet BC values
  void applyDirichletBoundaryConditionsInDisplacements(Data::FiniteElements<FunctionSpaceType,Term> &data);

  //! set entries in f to the entry in rhs for which Dirichlet BC are set
  void applyDirichletBoundaryConditionsInNonlinearFunction(Vec &f, Data::FiniteElements<FunctionSpaceType,Term> &data);

  //! set rows and columns in stiffness matrix to 0 for which boundary conditions are specified
  void applyDirichletBoundaryConditionsInStiffnessMatrix(std::shared_ptr<PartitionedPetscMat<FunctionSpaceType>> matrix, Data::FiniteElements<FunctionSpaceType,Term> &data);

  //! This transforms a 2D mesh input vector to a 3D mesh output vector by inserting 0's. It can only be called for 2D problems.
  void expandVectorTo3D(Vec &input, Vec &output, const int nLocalUnknowns3D);

protected:

  //! extract the submatrix that only contains entries for dofs that are not constraint by Dirichlet BCs
  void reduceMatrix(std::shared_ptr<PartitionedPetscMat<FunctionSpaceType>> input, std::shared_ptr<PartitionedPetscMat<FunctionSpaceType>> output, const int nLocalUnknowns);

  //! initialize Dirichlet boundary conditions
  void initializeBoundaryConditions(bool &externalVirtualWorkIsConstant, const int nLocalUnknowns, PythonConfig specificSettings, Data::FiniteElements<FunctionSpaceType,Term> &data);

  //! print boundary conditions
  void printBoundaryConditions();

  std::vector<dof_no_t> dirichletIndices_;  ///< the indices of unknowns (not dofs) for which the displacement is fixed
  std::vector<double> dirichletValues_;     ///< the to dirichletIndices corresponding fixed values for the displacement
  std::vector<double> zeros_;           ///< a vector of 0s, number of dirichlet values

  struct TractionBoundaryCondition
  {
    element_no_t elementGlobalNo;

    Mesh::face_t face;
    std::vector<std::pair<dof_no_t, VecD<FunctionSpaceType::dim()>>> dofVectors;  //<element-local dof no, value>

    // parse values from python config, e.g. {"element": 1, "face": "0+", "dofVectors:", {0: [tmax,0,0], 1: [tmax,0,0], 2: [tmax,0,0], 3: [tmax,0,0]}}
    TractionBoundaryCondition(PythonConfig specificSettings, std::shared_ptr<typename FunctionSpaceType::HighOrderFunctionSpace> mesh);
  };

  std::vector<TractionBoundaryCondition> tractionReferenceConfiguration_;   //< tractions for elements
  std::vector<TractionBoundaryCondition> tractionCurrentConfiguration_;

  std::vector<std::pair<element_no_t, VecD<FunctionSpaceType::dim()>>> bodyForceReferenceConfiguration_;  //< <element global no, vector>
  std::vector<std::pair<element_no_t, VecD<FunctionSpaceType::dim()>>> bodyForceCurrentConfiguration_;    //< <element global no, vector>

};

};  // namespace

#include "spatial_discretization/finite_element_method/solid_mechanics/solid_mechanics_boundary_conditions.tpp"
