c$Id: upremas.F,v 1.1 2006/11/21 16:44:39 rlt Exp $
      subroutine upremas(fl)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose:  Mass interface for PETSc

c     Inputs:
c        fl(1) - Form Consistent mass if true
c        fl(2) - Form Lumped mass if true

c     Outputs:
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

#     include   "include/finclude/petsc.h"
#     include   "include/finclude/petscvec.h"
#     include   "include/finclude/petscmat.h"
#     include   "pfeapa.h"
#     include   "pfeapc.h"

#     include   "cdata.h"
#     include   "comblk.h"
#     include   "endata.h"
#     include   "iofile.h"
#     include   "pointer.h"
#     include   "sdata.h"
#     include   "pfeapb.h"

      PetscErrorCode ierr
      PetscTruth     chk

      logical    fl(2)

      save

c     Preform setup

      if(fl(1)) then   ! Consistent mass allocate
        call MatValid(Mmat, chk, ierr)
        if(chk .eq. PETSC_FALSE) then
          if(pfeap_blk) then
            call MatCreateMPIBAIJ(PETSC_COMM_WORLD, ndf,numpeq,numpeq,
     &                            PETSC_DETERMINE, PETSC_DETERMINE,
     &                            PETSC_NULL_INTEGER, mr(np(246)),
     &                            PETSC_NULL_INTEGER, mr(np(247)),
     &                            Mmat, ierr )
          else
            call MatCreateMPIAIJ (PETSC_COMM_WORLD, numpeq, numpeq,
     &                            PETSC_DETERMINE, PETSC_DETERMINE,
     &                            PETSC_NULL_INTEGER, mr(np(246)),
     &                            PETSC_NULL_INTEGER, mr(np(247)),
     &                            Mmat, ierr )
          endif
        endif
        call MatSetFromOptions(Mmat, ierr)
      elseif(fl(2)) then   ! Lumped mass allocate
        call VecValid(Mdiag, chk, ierr)
        if(chk .eq. PETSC_FALSE) then
          call VecCreate        (PETSC_COMM_WORLD, Mdiag, ierr)
          call VecSetSizes      (Mdiag, numpeq, PETSC_DECIDE, ierr)
          call VecSetFromOptions(Mdiag, ierr)
        endif
      elseif(.not.fl(1) .and. .not.fl(2)) then
        write(iow,*) ' ERROR DID NOT ALLOCATE MASS MATRIX'
        write(  *,*) ' ERROR DID NOT ALLOCATE MASS MATRIX'
        call plstop()
      endif
      if(ierr .ne. 0) then
        write(iow,*) 'Error on MatCreate'
        write(  *,*) 'Error on MatCreate'
        call plstop()
      endif

      if(fl(1)) then
        call MatZeroEntries   (Mmat,ierr)
      else
        call VecZeroEntries   (Mdiag,ierr)
      endif

c     Vector for matrix multiply

      call VecValid(xvec, chk, ierr)
      if(chk .eq. PETSC_FALSE) then
        call VecCreate        (PETSC_COMM_WORLD, xvec, ierr)
        call VecSetSizes      (xvec, numpeq, PETSC_DECIDE, ierr)
        call VecSetFromOptions(xvec, ierr)
      endif

      end
