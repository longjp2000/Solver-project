c$Id: ppmodin.F,v 1.1 2006/11/21 16:44:39 rlt Exp $
      subroutine ppmodin(phi, y, u,v,tu,nv,tneq,mtyp)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Initial Conditions for modal solution of
c               linear transient problems

c      Inputs:
c         phi(vneq,*) - Mass orthonormal eigenvectors
c         u(*)        - Initial displacements
c         v(*)        - Initial velocities
c         nv          - Number of eigenpairs (converged)
c         tneq        - Number of d.o.f. in model
c         mtyp        - Type of mass (true for consistent,false for lumped)

c      Scratch:
c         tu(*,2)     - Temporary vector

c      Outputs:
c         y(nv,3)     - Eigensolution at time t_0
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

#     include  "include/finclude/petsc.h"
#     include  "include/finclude/petscvec.h"
#     include  "include/finclude/petscmat.h"
#     include  "pfeapa.h"
#     include  "pfeapc.h"

#     include  "pfeapb.h"

#     include  "sdata.h"
#     include  "pointer.h"
#     include  "comblk.h"

      logical   mtyp
      integer   i,j,k,l,n, nv,tneq
      real*8    phi(vneq,*), y(nv,3), u(*),v(*), tu(tneq,*)

c     Petsc values

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

      do j = 1,2
        l = 0
        do i = 0,numpn*ndf-1
          k = mr(np(245)+i) - 1
          if(k.ge.0) then
            l = l + 1
            call VecSetValue(xvec,k,tu(l,j),INSERT_VALUES, ierr)
          endif
        end do ! i

        call VecAssemblyBegin(xvec, ierr)
        call VecAssemblyEnd  (xvec, ierr)

        if(mtyp) then
          call MatMult( Mmat, xvec, yvec, ierr )
        else
          call VecPointwiseMult( yvec, Mdiag, xvec, ierr )
        endif

        do n = 1,nv
          l = 0
          do i = 0,numpn*ndf-1
            k = mr(np(245)+i) - 1
            if(k.ge.0) then
              l = l + 1
              call VecSetValue(xvec,k,phi(l,n),INSERT_VALUES, ierr)
            endif
          end do ! i

          call VecAssemblyBegin(xvec, ierr)
          call VecAssemblyEnd  (xvec, ierr)

          call VecDot(xvec,yvec, y(n,j), ierr)
        end do ! n
      end do ! j

      end
