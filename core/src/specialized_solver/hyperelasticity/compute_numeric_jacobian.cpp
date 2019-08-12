#include "specialized_solver/hyperelasticity/compute_numeric_jacobian.h"

#include <petsc/private/snesimpl.h>
#include <petscdm.h>
#include <iomanip>

#include "easylogging++.h"

PetscErrorCode  SNESComputeJacobianDefaultNested(SNES snes,Vec x1,Mat J,Mat B,void *ctx)
{
  Vec               j1a,j2a,x2;
  PetscErrorCode    ierr;
  PetscInt          i,N,start,end,j,value,root,n_nested,n_nested0, n_nested1,i_row_no,i_col_no;
  PetscScalar       dx,*y,wscale;
  const PetscScalar *xx;
  PetscReal         amax,epsilon = PETSC_SQRT_MACHINE_EPSILON;
  PetscReal         dx_min = 1.e-16,dx_par = 1.e-1,unorm;
  MPI_Comm          comm;
  PetscBool         assembled,use_wp = PETSC_TRUE,flg;
  const char        *list[2] = {"ds","wp"};
  PetscMPIInt       size;
  const PetscInt    *ranges;
  Vec               *nested_x1, *nested_x2, current_x1, current_x2, *nested_j1a, *nested_j2a;
  Mat               **nested_B, current_B;

  const bool debug = false;

  // determine if matrix is nested
  MatType type;
  MatGetType(B,&type);
  if (std::string(type) != std::string(MATNEST))
    return SNESComputeJacobianDefault(snes, x1, J, B, ctx);

  /* Since this Jacobian will possibly have "extra" nonzero locations just turn off errors for these locations */
  MatSetOption(B,MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
  PetscOptionsGetReal(((PetscObject)snes)->options,((PetscObject)snes)->prefix,"-snes_test_err",&epsilon,0);

  PetscObjectGetComm((PetscObject)x1,&comm);
  MPI_Comm_size(comm,&size);
  MatAssembled(B,&assembled);
  if (assembled) {
    MatZeroEntries(B);
  }
  if (!snes->nvwork) {
    if (snes->dm) {
      DMGetGlobalVector(snes->dm,&j1a);
      DMGetGlobalVector(snes->dm,&j2a);
      DMGetGlobalVector(snes->dm,&x2);
    } else {
      snes->nvwork = 3;
      VecDuplicateVecs(x1,snes->nvwork,&snes->vwork);
      PetscLogObjectParents(snes,snes->nvwork,snes->vwork);
      j1a = snes->vwork[0]; j2a = snes->vwork[1]; x2 = snes->vwork[2];
    }
  }

  //MatNestGetSubMats(Mat A,PetscInt *M,PetscInt *N,Mat ***mat)

  MatNestGetSubMats(B,&n_nested,&n_nested0,&nested_B);

  if (n_nested != n_nested0)
  {
    LOG(FATAL) << n_nested << "!=" << n_nested0;
  }

  for (i_row_no=0; i_row_no<n_nested; i_row_no++)
  {
    for (i_col_no=0; i_col_no<n_nested; i_col_no++)
    {
      // set values in B
      current_B = nested_B[i_row_no][i_col_no];
      /* Since this Jacobian will possibly have "extra" nonzero locations just turn off errors for these locations */
      MatSetOption(current_B,MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
    }
  }


  VecNestGetSubVecs(x1,&n_nested1,&nested_x1);

  if (n_nested != n_nested1)
  {
    LOG(FATAL) << n_nested << "!=" << n_nested1;
  }

  nested_x2 = new Vec[n_nested];
  nested_j1a = new Vec[n_nested];
  nested_j2a = new Vec[n_nested];
  for(i=0; i<n_nested; i++)
  {
    VecDuplicate(nested_x1[i],&nested_x2[i]);
    VecDuplicate(nested_x1[i],&nested_j1a[i]);
    VecDuplicate(nested_x1[i],&nested_j2a[i]);
  }

  PetscOptionsBegin(PetscObjectComm((PetscObject)snes),((PetscObject)snes)->prefix,"Differencing options","SNES");
  PetscOptionsEList("-mat_fd_type","Algorithm to compute difference parameter","SNESComputeJacobianDefault",list,2,"wp",&value,&flg);
  PetscOptionsEnd();
  if (flg && !value) use_wp = PETSC_FALSE;

  if (use_wp) {
    // compute norm of x1
    PetscReal unorm_part = 0;
    unorm = 0;
    for (i = 0; i < n_nested; i++)
    {
      VecNorm(nested_x1[i],NORM_2,&unorm_part);
      unorm_part = unorm_part*unorm_part;
      unorm += unorm_part;
    }
    unorm = sqrt(unorm);
  }

  // create nested vec
  VecCreateNest(comm, n_nested, NULL, nested_x2, &x2);
  VecCreateNest(comm, n_nested, NULL, nested_j1a, &j1a);
  VecCreateNest(comm, n_nested, NULL, nested_j2a, &j2a);

  // compute reference value
  SNESComputeFunction(snes,x1,j1a);

  // iterate over block columns of the resulting matrix or equivalently over block rows of the variable vector
  for (i_col_no=0; i_col_no<n_nested; i_col_no++)
  {
    // col no.: to what variable we derive
    // row no.: what part of the output will be put into the jacobian

    // part of the variable vector that will be modified
    current_x1 = nested_x1[i_col_no];

    if (debug)
    {
      std::stringstream s;
      s<<"j1a F(";
      for(int k=0; k<n_nested;k++)
      {
        PetscInt x1size;
        VecGetLocalSize(nested_x1[k],&x1size);
        VecGetOwnershipRange(nested_x1[k],&start,&end);
        for(int l=start; l<end; l++)
        {
          PetscReal value;
          VecGetValues(nested_x1[k],1,&l,&value);
          s<<std::setprecision(10) << value;
          if (l != x1size-1)
            s<<",";
        }
        if (k != n_nested-1)
          s<<"|";
      }
      s << ") " <<std::endl<<" = ";
      for(int k=0; k<n_nested;k++)
      {
        PetscInt j1asize;
        VecGetLocalSize(nested_j1a[k],&j1asize);
        VecGetOwnershipRange(nested_j1a[k],&start,&end);
        for(int l=start; l<end; l++)
        {
          PetscReal value;
          VecGetValues(nested_j1a[k],1,&l,&value);
          s<<std::setprecision(10) << value;
          if (l != j1asize-1)
            s<<",";
        }
        if (k != n_nested-1)
          s<<"|";
      }
      s << std::endl;
      LOG(DEBUG) << s.str();
    }

    /* Compute Jacobian approximation, 1 column at a time.
        x1 = current iterate, j2a = F(x1)
        x2 = perturbed iterate, j2a = F(x2)
    */
    VecGetSize(current_x1,&N);
    // loop over columns of current matrix part
    for (i=0; i<N; i++) {

      // reset x2 as copy of x1
      for (j=0; j<n_nested; j++)
      {
        VecCopy(nested_x1[j],nested_x2[j]);
      }

      // now modify nested_x2[i_col_no,i] at (block row i_col_no, row i), this corresponds to entry (block column i_col_no, col i) in the jacobian
      current_x2 = nested_x2[i_col_no];
      VecGetOwnershipRange(current_x2,&start,&end);
      if (i>= start && i<end) {
        VecGetArrayRead(current_x1,&xx);
        if (use_wp) dx = PetscSqrtReal(1.0 + unorm);
        else        dx = xx[i-start];
        VecRestoreArrayRead(current_x1,&xx);
        if (PetscAbsScalar(dx) < dx_min) dx = (PetscRealPart(dx) < 0. ? -1. : 1.) * dx_par;
        dx    *= epsilon;
        wscale = 1.0/dx;
        VecSetValues(current_x2,1,&i,&dx,ADD_VALUES);
      } else {
        wscale = 0.0;
      }
      VecAssemblyBegin(current_x2);
      VecAssemblyEnd(current_x2);

      // compute function at modified value x2
      SNESComputeFunction(snes,x2,j2a);

      if (debug)
      {
        std::stringstream s;
        s<<"j2a F(";
        for(int k=0; k<n_nested;k++)
        {
          PetscInt x2size;
          VecGetLocalSize(nested_x2[k],&x2size);
          VecGetOwnershipRange(nested_x2[k],&start,&end);
          for(int l=start; l<end; l++)
          {
            PetscReal value;
            VecGetValues(nested_x2[k],1,&l,&value);
            s<<std::setprecision(10) << value;
            if (l != x2size-1)
              s<<",";
          }
          if (k != n_nested-1)
            s<<"|";
        }
        s << ") " <<std::endl<<" = ";
        for(int k=0; k<n_nested;k++)
        {
          PetscInt j2asize;
          VecGetLocalSize(nested_j2a[k],&j2asize);
          VecGetOwnershipRange(nested_j2a[k],&start,&end);
          for(int l=start; l<end; l++)
          {
            PetscReal value;
            VecGetValues(nested_j2a[k],1,&l,&value);
            s<<std::setprecision(10) << value;
            if (l != j2asize-1)
              s<<",";
          }
          if (k != n_nested-1)
            s<<"|";
        }
        s << std::endl;
        LOG(DEBUG) << s.str();
      }

      // compute j2a -= j1a for nested vec
      //VecAXPY(j2a,-1.0,j1a);
      for (j=0; j<n_nested; j++)
      {
        VecAXPY(nested_j2a[j],-1.0,nested_j1a[j]);
      }

      /* Communicate scale=1/dx_i to all processors */
      VecGetOwnershipRanges(current_x2,&ranges);

      // wscale is only set on that rank where x2 was modified. This is the rank that owns dof i
      root = size;  // n ranks
      for (j=size-1; j>-1; j--) {
        root--;
        if (i>=ranges[j]) break;
      }
      MPI_Bcast(&wscale,1,MPIU_SCALAR,root,comm);

      // scale whole nested j2a with wscale=1/dx, then j2a will be the column of the jacobian
      // and determine the maximum entry of the jacobian in this column, amax
      PetscReal amax_part = 0;
      amax = 0;
      for (j=0; j<n_nested; j++)
      {
        VecScale(nested_j2a[j],wscale);
        VecNorm(nested_j2a[j],NORM_INFINITY,&amax_part);
        amax = std::max(amax_part, amax);
      }
      amax *= 1.e-14;

      // loop over block columns of jacobian
      for (i_row_no=0; i_row_no<n_nested; i_row_no++)
      {
        // set values in B
        current_B = nested_B[i_row_no][i_col_no];
        /*PetscInt m, n;
        MatGetLocalSize(current_B, &m, &n);*/

        VecGetOwnershipRange(nested_j2a[i_row_no],&start,&end);

        VecGetArray(nested_j2a[i_row_no],&y);
        // loop over rows in block row
        for (j=start; j<end; j++) {
          if (PetscAbsScalar(y[j-start]) > amax || j == i) {
            PetscReal value = y[j-start];

            if (debug)
            {
              LOG(DEBUG) << "B[" << i_row_no << "," << i_col_no << "] at " << j << "," << i << ", value " << std::setprecision(10) << value << std::endl;
            }

            MatSetValues(current_B,1,&j,1,&i,&value,INSERT_VALUES);
          }
        }
        VecRestoreArray(nested_j2a[i_row_no],&y);
      }
    }
  }
  if (snes->dm) {
    DMRestoreGlobalVector(snes->dm,&j1a);
    DMRestoreGlobalVector(snes->dm,&j2a);
    DMRestoreGlobalVector(snes->dm,&x2);
  }
  MatAssemblyBegin(B,MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(B,MAT_FINAL_ASSEMBLY);
  if (B != J) {
    MatAssemblyBegin(J,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(J,MAT_FINAL_ASSEMBLY);
  }
  delete [] nested_x2;
  delete [] nested_j1a;
  delete [] nested_j2a;

  return(0);
}
