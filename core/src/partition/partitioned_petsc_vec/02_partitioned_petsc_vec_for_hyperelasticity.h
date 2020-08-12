#pragma once

#include <memory>

#include "control/types.h"
#include "partition/rank_subset.h"
#include "function_space/function_space.h"
#include "mesh/type_traits.h"
#include "partition/partitioned_petsc_vec/01_partitioned_petsc_vec_with_dirichlet_bc.h"
#include "spatial_discretization/dirichlet_boundary_conditions/01_dirichlet_boundary_conditions.h"

/** This vector is an extension of PartitionedPetscVecWithDirichletBc (read the description there).
 *  In its standard version (nComponents == 3), it is a combined vector of 3 displacements components and 1 pressure component.
 *  If nComponents == 6, it also has velocity components, so 6 components in total (3 displacements, 3 velocities) + 1 pressure.
 *
 *  for the incompressible formulation, where the pressure is an unkown:
 *    for nComponents = 3:     for nComponents = 6:
 *
 *      (ux)                        (ux)
 *      (uy)                        (uy)
 *      (uz)                        (uz)
 *      ( p)                        (vx)
 *                                  (vy)
 *                                  (vz)
 *                                  ( p)
 *
 *  for the compressible formulation, where the pressure is given by a constitutive equation:
 *    for nComponents = 3:     for nComponents = 6:
 *
 *      (ux)                        (ux)
 *      (uy)                        (uy)
 *      (uz)                        (uz)
 *                                  (vx)
 *                                  (vy)
 *                                  (vz)
 *
 *  The template parameter nComponents is the number of components without pressure, i.e. 3 for displacements only or 6 for displacements and velocities.
 */
template<typename DisplacementsFunctionSpaceType, typename PressureFunctionSpaceType, typename Term, int nComponents=3, typename DummyForTraits=Term>
class PartitionedPetscVecForHyperelasticity
{};

/** Partition specialization for incompressible formulation, which has the (p) component at the end
 */
template<typename DisplacementsFunctionSpaceType, typename PressureFunctionSpaceType, typename Term, int nComponents>
class PartitionedPetscVecForHyperelasticity<DisplacementsFunctionSpaceType,PressureFunctionSpaceType,Term,nComponents,std::enable_if_t<Term::isIncompressible,Term>>:
  public PartitionedPetscVecWithDirichletBc<DisplacementsFunctionSpaceType,nComponents+1,nComponents>  // <DisplacementsFunctionSpaceType,4,3>
{
public:
 
  //! constructor
  PartitionedPetscVecForHyperelasticity(
    std::shared_ptr<Partition::MeshPartition<DisplacementsFunctionSpaceType>> meshPartitionDisplacements,
    std::shared_ptr<Partition::MeshPartition<PressureFunctionSpaceType>> meshPartitionPressure,
    std::shared_ptr<SpatialDiscretization::DirichletBoundaryConditions<DisplacementsFunctionSpaceType,nComponents>> dirichletBoundaryConditions,
    std::string name
  );

  //! write the vector values to a file in natural ordering, such that the output for different numbers of ranks can be compared
  void dumpGlobalNatural(std::string filename);

  //! the meshPartition of the pressure function space, the other meshPartition, the one of the displacements function space can be obtained by meshPartition()
  std::shared_ptr<Partition::MeshPartition<PressureFunctionSpaceType>> meshPartitionPressure();

  //! get a string representation of the whole vector in global natural ordering
  //! @param horizontal if the string is for console output (less newlines) or for file output
  std::string getString(bool horizontal = true, std::string vectorName = "") const;

  //! get the number of displacement dofs ux,uy,uz without Dirichlet Bc constraint dofs
  dof_no_t nDisplacementDofsWithoutBcLocal();

  //! get the number of velocity dofs vx,vy,vz without Dirichlet Bc constraint dofs (if nComponents == 6), or 0 if there are no velocities (if nComponents == 3)
  dof_no_t nVelocityDofsWithoutBcLocal();

  //! check if the vector contains nan or inf values and output error
  virtual bool containsNanOrInf();

  //! get a Petsc index set (IS) of all global non-bc indices comprising ux,uy,uz.
  //! This is a sequence of consecutive numbers (e.g. 5,6,7,...) and can be used to access submatrices of the global system matrix.
  IS displacementDofsGlobal();

  //! get a Petsc index set (IS) of all global non-bc indices comprising vx,vy,vz, only if nComponents == 6
  IS velocityDofsGlobal();

  //! get a Petsc index set (IS) of all global non-bc indices p
  IS pressureDofsGlobal();

protected:

  //! do the additional initalization after the initialization of the parent class, that accounts also for the pressure component
  void initializeForPressure();

  std::shared_ptr<Partition::MeshPartition<PressureFunctionSpaceType>> meshPartitionPressure_;   //< the mesh partition for the pressure function space, which has a lower degree of ansatz function than the displacement function space
  int componentNoPressure_ = nComponents;     //< which of the component is dedicated for the pressure, this is the last component
};

/** Partition specialization for compressible formulation, without the (p) component at the end
 */
template<typename DisplacementsFunctionSpaceType, typename PressureFunctionSpaceType, typename Term, int nComponents>
class PartitionedPetscVecForHyperelasticity<DisplacementsFunctionSpaceType,PressureFunctionSpaceType,Term,nComponents,std::enable_if_t<!Term::isIncompressible,Term>>:
  public PartitionedPetscVecWithDirichletBc<DisplacementsFunctionSpaceType,nComponents,nComponents>  // <DisplacementsFunctionSpaceType,3,3>
{
public:

  //! constructor
  PartitionedPetscVecForHyperelasticity(
    std::shared_ptr<Partition::MeshPartition<DisplacementsFunctionSpaceType>> meshPartitionDisplacements,
    std::shared_ptr<Partition::MeshPartition<PressureFunctionSpaceType>> meshPartitionPressure,
    std::shared_ptr<SpatialDiscretization::DirichletBoundaryConditions<DisplacementsFunctionSpaceType,nComponents>> dirichletBoundaryConditions,
    std::string name
  );

  //! write the vector values to a file in natural ordering, such that the output for different numbers of ranks can be compared
  void dumpGlobalNatural(std::string filename);

  //! the meshPartition of the pressure function space, the other meshPartition, the one of the displacements function space can be obtained by meshPartition()
  std::shared_ptr<Partition::MeshPartition<PressureFunctionSpaceType>> meshPartitionPressure();

  //! get a string representation of the whole vector in global natural ordering
  //! @param horizontal if the string is for console output (less newlines) or for file output
  std::string getString(bool horizontal = true, std::string vectorName = "") const;

  //! get the number of displacement dofs ux,uy,uz without Dirichlet Bc constraint dofs
  dof_no_t nDisplacementDofsWithoutBcLocal();

  //! get the number of velocity dofs vx,vy,vz without Dirichlet Bc constraint dofs (if nComponents == 6), or 0 if there are no velocities (if nComponents == 3)
  dof_no_t nVelocityDofsWithoutBcLocal();

  //! check if the vector contains nan or inf values and output error
  virtual bool containsNanOrInf();

  //! get a Petsc index set (IS) of all global non-bc indices comprising ux,uy,uz.
  //! This is a sequence of consecutive numbers (e.g. 5,6,7,...) and can be used to access submatrices of the global system matrix.
  IS displacementDofsGlobal();

  //! get a Petsc index set (IS) of all global non-bc indices comprising vx,vy,vz, only if nComponents == 6
  IS velocityDofsGlobal();

  //! get a Petsc index set (IS) of all global non-bc indices p
  IS pressureDofsGlobal();

protected:

  std::shared_ptr<Partition::MeshPartition<PressureFunctionSpaceType>> meshPartitionPressure_;   //< the mesh partition for the pressure function space, which has a lower degree of ansatz function than the displacement function space
};

template<typename DisplacementsFunctionSpaceType, typename PressureFunctionSpaceType, typename Term>
std::ostream &operator<<(std::ostream &stream, const PartitionedPetscVecForHyperelasticity<DisplacementsFunctionSpaceType,PressureFunctionSpaceType,Term> &vector);

#include "partition/partitioned_petsc_vec/02_partitioned_petsc_vec_for_hyperelasticity_compressible.tpp"
#include "partition/partitioned_petsc_vec/02_partitioned_petsc_vec_for_hyperelasticity_incompressible.tpp"
