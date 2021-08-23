c$Id: paropm.F,v 1.1 2006/11/20 20:43:46 rlt Exp $
      subroutine paropm(pneq,u,v)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Compute product of 'b' with 'v'

c      Inputs:
c         pneq     - Number of partition equations
c         u(*)     - Iteration vectors

c      Outputs:
c         v(*)     - B * u
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

#     include  "include/finclude/petsc.h"
#     include  "include/finclude/petscvec.h"
#     include  "include/finclude/petscmat.h"
#     include  "pfeapa.h"
#     include  "pfeapb.h"
#     include  "pfeapc.h"

#     include  "fdata.h"
#     include  "sdata.h"
#     include  "pointer.h"
#     include  "comblk.h"

      integer   pneq, i,j,k
      real*8    u(*),v(*)

c     Petsc values

      PetscScalar    sol_array(1)
      PetscOffset    i_sol
      PetscTruth     chk
      PetscInt       starti,endi
      PetscErrorCode ierr

      save

c     Set global PETSc vectors

      call VecValid(xvec, chk, ierr)
      if(chk .eq. PETSC_FALSE) then
        call VecCreate        (PETSC_COMM_WORLD, xvec, ierr)
        call VecSetSizes      (xvec, numpeq, PETSC_DECIDE, ierr)
        call VecSetFromOptions(xvec, ierr)
        call VecDuplicate     (xvec, yvec, ierr)
      endif
      call VecValid(yvec, chk, ierr)
      if(chk .eq. PETSC_FALSE) then
        call VecDuplicate     (xvec, yvec, ierr)
      endif

c     Multiply u(*) by M

      j = 0
      do i = 0,numpn*ndf
        k = mr(np(245)+i) - 1
        if(k.ge.0) then
          j = j + 1
          call VecSetValue     (xvec,k,u(j),INSERT_VALUES, ierr)
        endif
      end do ! i

      call VecAssemblyBegin(xvec, ierr)
      call VecAssemblyEnd  (xvec, ierr)

      if(fl(1)) then
        call MatMult           (Mmat, xvec, yvec, ierr)
      else
        call VecPointwiseMult  (yvec, Mdiag, xvec, ierr)
      endif

c     Set my values in v(1)

      call VecGetOwnershipRange(yvec, starti, endi, ierr)
      call VecGetArray         (yvec, sol_array, i_sol, ierr)

c     do i = 1,pneq
      do i = 1,j
        v(i) = sol_array(i_sol+i)
      end do ! i

      call VecRestoreArray     (yvec, sol_array, i_sol, ierr)

      end
