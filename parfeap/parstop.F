c$Id: parstop.F,v 1.1 2006/11/21 16:44:39 rlt Exp $
      subroutine parstop()

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Close any open parallel array and delete memory use
c               PETSc version

c      Inputs:
c         none

c      Outputs:
c         none
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

#     include   "include/finclude/petsc.h"
#     include   "include/finclude/petscmat.h"
#     include   "include/finclude/petscvec.h"
#     include   "pfeapc.h"

      PetscErrorCode ierr
      PetscTruth     chk

      save

c     Close out PETSc matrices

      call MatValid(Kmat, chk, ierr)
      if(chk .eq. PETSC_TRUE) then
        call MatDestroy(Kmat, ierr)
      endif
      call MatValid(Pmat, chk, ierr)
      if(chk .eq. PETSC_TRUE) then
        call MatDestroy(Pmat, ierr)
      endif
      call MatValid(Mmat, chk, ierr)
      if(chk .eq. PETSC_TRUE) then
        call MatDestroy(Mmat, ierr)
      endif

c     Close out PETSc vectors

      call VecValid(Mdiag, chk, ierr)
      if(chk .eq. PETSC_TRUE) then
        call VecDestroy(Mdiag, ierr)
      endif
      call VecValid(Msqrt, chk, ierr)
      if(chk .eq. PETSC_TRUE) then
        call VecDestroy(Msqrt, ierr)
      endif
      call VecValid(xvec, chk, ierr)
      if(chk .eq. PETSC_TRUE) then
        call VecDestroy(xvec, ierr)
      endif
      call VecValid(yvec, chk, ierr)
      if(chk .eq. PETSC_TRUE) then
        call VecDestroy(yvec, ierr)
      endif
      call VecValid(zvec, chk, ierr)
      if(chk .eq. PETSC_TRUE) then
        call VecDestroy(zvec, ierr)
      endif
      call VecValid(rhs , chk, ierr)
      if(chk .eq. PETSC_TRUE) then
        call VecDestroy(rhs , ierr)
      endif
      call VecValid(sol , chk, ierr)
      if(chk .eq. PETSC_TRUE) then
        call VecDestroy(sol , ierr)
      endif

c     Stop PETSc

      call PetscFinalize(ierr)

      end
