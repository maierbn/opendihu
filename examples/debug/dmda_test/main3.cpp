#include <iostream>
#include <cstdlib>
#include <assert.h>

#include <petscdmda.h>  
  
  
int main(int argc, char *argv[])
{
  // initialize
  MPI_Init(&argc, &argv);
  PetscInitialize(&argc, &argv, NULL, NULL);
  
  // get own rank and size of MPI_COMM_WORLD
  PetscMPIInt ownRankNo, size;
  MPI_Comm_rank(PETSC_COMM_WORLD, &ownRankNo);
  MPI_Comm_size(PETSC_COMM_WORLD,&size);
  
  const PetscInt M = 5;
  
  DM da;
  PetscErrorCode ierr;
  
  // create 1d decomposition
  ierr = DMDACreate1d(PETSC_COMM_WORLD, DM_BOUNDARY_NONE, M, 1, 1,
                      NULL, &da); CHKERRQ(ierr);

  // get global coordinates of local partition
  PetscInt x, m;
  ierr = DMDAGetCorners(da, &x, NULL, NULL, &m, NULL, NULL); CHKERRQ(ierr);
  std::cout << ownRankNo << "/" << size << ": own range is from x=" << x << ", with size m=" << m << std::endl;

  // get local sizes on the ranks
  const PetscInt *lxData;
  ierr = DMDAGetOwnershipRanges(da, &lxData, NULL, NULL); CHKERRQ(ierr);

  std::cout << ownRankNo << "/" << size << ": lx: " << lxData << std::endl;
  assert(lxData != NULL);
  
  PetscFinalize();
  MPI_Finalize();
  
  return EXIT_SUCCESS;
}
