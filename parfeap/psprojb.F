c$Id: psprojb.F,v 1.1 2006/11/21 16:44:39 rlt Exp $
      subroutine psprojb(v,t,h,pneq,nv)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Compute subspace projection of 'b' to form 'h'

c      Inputs:
c         v(neq,*) - Set of iteration vectors
c         pneq      - Number of partition equations
c         nv        - Size of projected matrix

c      Scratch:
c         t(vneq)   - Working vector

c      Outputs:
c         h(*)     - Projected matrix V_trans * B * V
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

#     include  "include/finclude/petsc.h"
#     include  "include/finclude/petscvec.h"
#     include  "include/finclude/petscmat.h"
#     include  "pfeapa.h"
#     include  "pfeapc.h"

#     include  "pfeapb.h"
#     include  "fdata.h"
#     include  "sdata.h"
#     include  "pointer.h"
#     include  "comblk.h"

      integer   pneq,nv, i,j,k,l
      real*8    dot,v(vneq,*),t(*),h(*)

c     Petsc values

      PetscScalar    sol_array(1)
      PetscOffset    i_sol
      PetscInt       starti,endi
      PetscTruth     chk
      PetscErrorCode ierr

      save

c     Compute 'z' and 'b' projection to form 'h'

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

c     Multiply v(*,j) by M

      do j = 1,nv
        l = 0
        do i = 0,numpn*ndf
          k = mr(np(245)+i) - 1
          if(k.ge.0) then
            l = l + 1
            call VecSetValue(xvec,k,v(l,j),INSERT_VALUES, ierr)
          endif
        end do ! i

        call VecAssemblyBegin(xvec, ierr)
        call VecAssemblyEnd  (xvec, ierr)

        if(fl(1)) then
          call MatMult( Mmat, xvec, yvec, ierr )
        else
          call VecPointwiseMult( yvec, Mdiag, xvec, ierr )
        endif

c       Set my values in v(1,j)

        call VecGetOwnershipRange(yvec, starti, endi, ierr)
        call VecGetArray         (yvec, sol_array, i_sol, ierr)

        do i = 1,numpeq
          t(i) = sol_array(i_sol+i)
        end do ! i

        call VecRestoreArray     (yvec, sol_array, i_sol, ierr)

c       Project 'z' and 'v' vectors to form 'h'

        k = j*(j+1)/2
        do i = j,nv
          h(k) = dot(t,v(1,i),numpeq)
          k = k + i
        end do ! i
        do i = 1,numpeq
          v(i,j) = t(i)
        end do ! i
      end do ! j

      end
