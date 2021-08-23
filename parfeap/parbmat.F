c$Id: parbmat.F,v 1.1 2006/11/21 16:44:39 rlt Exp $
      subroutine parbmat(phib,mb,cmass,lmass,t, neq)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose: Compute base matrix coupling term

c     Inputs:

c     Outputs:

c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

#     include   "include/finclude/petsc.h"
#     include   "include/finclude/petscvec.h"
#     include   "include/finclude/petscmat.h"
#     include   "pfeapa.h"
#     include   "pfeapb.h"
#     include   "pfeapc.h"

#     include   "sdata.h"

c     Petsc values

      PetscScalar    sol_array(1)
      PetscOffset    i_sol
      PetscInt       starti,endi
      PetscTruth     chk
      PetscErrorCode ierr

#     include   "fdata.h"
#     include   "pointer.h"
#     include   "comblk.h"

      integer    neq, i,j,k
      real*8     t(*),cmass(*),lmass(*),phib(*),mb(*)

      save

c     Parallel solution

c     Set PETSc vectors

      call VecValid(xvec, chk, ierr)
      if(chk .eq. PETSC_FALSE) then
        call VecCreate        (PETSC_COMM_WORLD, xvec, ierr)
        call VecSetSizes      (xvec, numpeq, PETSC_DECIDE, ierr)
        call VecSetFromOptions(xvec, ierr)
      endif
      call VecValid(yvec, chk, ierr)
      if(chk .eq. PETSC_FALSE) then
        call VecDuplicate     (xvec, yvec, ierr)
      endif

c     Multiply phib(*) by M

      j = 0
      do i = 0,ndf*numpn-1
        k = mr(np(245)+i) - 1
        if(k.ge.0) then
          j = j + 1 
          call VecSetValue(xvec,k,phib(j),INSERT_VALUES, ierr)
        endif
      end do ! i
      call VecAssemblyBegin(xvec, ierr)
      call VecAssemblyEnd  (xvec, ierr)

      if(fl(1)) then
        call MatMult( Mmat, xvec, yvec, ierr )
      else
        call VecPointwiseMult( yvec, Mdiag, xvec, ierr )
      endif

c     Set my values in t(*)

      call VecGetOwnershipRange(yvec, starti, endi, ierr)
      call VecGetArray         (yvec, sol_array, i_sol, ierr)

      do i = 1,numpeq
        t(i) = sol_array(i_sol+i) - mb(i)
      end do ! i

      call VecRestoreArray     (yvec, sol_array, i_sol, ierr)

      end
