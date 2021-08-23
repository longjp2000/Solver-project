c$Id: usolve.F,v 1.2 2006/11/25 02:12:48 rlt Exp $
      subroutine usolve(flags,b)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Revise send/receive (pmaxr,pmaxs use)            24/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose:  Solver interface for PETSc

c     Inputs:
c        flags(1) - Allocation and/or initialization phase
c        flags(2) - Perform factorization for direct solutions
c        flags(3) - Coefficient array unsymmetric
c        flags(4) - Solve equations
c        flags(5) - Purge storage of pointers
c        b(*)     - RHS vector

c     Outputs:

c        flags(5) - True if error occurs (for factor/solve only)
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

#     include   "include/finclude/petsc.h"
#     include   "include/finclude/petscpc.h"
#     include   "include/finclude/petscvec.h"
#     include   "include/finclude/petscmat.h"
#     include   "include/finclude/petscksp.h"
#     include   "include/finclude/petscviewer.h"
#     include   "pfeapa.h"
#     include   "pfeapb.h"
#     include   "pfeapc.h"

#     include   "cdata.h"
#     include   "comblk.h"
#     include   "compas.h"
#     include   "debugs.h"
#     include   "endata.h"
#     include   "iofile.h"
#     include   "pointer.h"
#     include   "print.h"
#     include   "sdata.h"

      PetscScalar    sol_array(1)
      PetscOffset    i_sol
      PetscInt       starti,endi, iters
      PetscErrorCode ierr
      PetscTruth     chk
      PetscViewer    K_view, B_view
      KSP            kspsol
      KSPConvergedReason reason
      PC             pc
      PetscReal      p_rtol,p_atol,p_dtol
      PetscInt       p_maxit

      real*8         info(MAT_INFO_SIZE), nza,nzr,mal

      logical    flags(5), setvar, palloc
      character  creason(7)*29,nreason(10)*26, pctype*20 

      external   FeapUserPCApply

      integer    i, j, k, pmaxr, pmaxs
      real*8     b(*)

      save

      data       creason / ' Unknown -1',
     &                     ' Satisfied residual tolerance',
     &                     ' Satisfied absolute tolerance',
     &                     ' Satisfied iteration limit',
     &                     ' Unused - 5',
     &                     ' Unused - 6',
     &                     ' Satisfied step length'/

      data       nreason / ' Unknown ',
     &                     ' Solution null',
     &                     ' Iterations exceeded',
     &                     ' Solution diverged',
     &                     ' Breakdown?',
     &                     ' Breakdown BiCG?',
     &                     ' Matrix nonsymmetric',
     &                     ' Indefinite Preconditioner',
     &                     ' NaN in solution',
     &                     ' Indefinite Matrix'/

c     Presolve setups

      if(flags(1)) then

        call MatValid(Kmat, chk, ierr)
        if(chk .eq. PETSC_FALSE) then
          if(pfeap_blk) then
            call MatCreateMPIBAIJ(PETSC_COMM_WORLD,nsbk, numpeq, numpeq,
     &                            PETSC_DETERMINE, PETSC_DETERMINE,
     &                            PETSC_NULL_INTEGER, mr(np(246)),
     &                            PETSC_NULL_INTEGER, mr(np(247)),
     &                            Kmat, ierr )
          else
            call MatCreateMPIAIJ (PETSC_COMM_WORLD, numpeq, numpeq,
     &                            PETSC_DETERMINE, PETSC_DETERMINE,
     &                            PETSC_NULL_INTEGER, mr(np(246)),
     &                            PETSC_NULL_INTEGER, mr(np(247)),
     &                            Kmat, ierr )
          endif
          if(ierr .ne. 0) then
            flags(5) = .true.
            write(*,*) 'Error on MatCreate'
          endif
          call MatSetFromOptions(Kmat, ierr)
        endif
        call MatZeroEntries     (Kmat,ierr)

c     Solution steps for assembled equations

      else

c       Factor equations

        if(flags(2)) then

          call MatAssemblyBegin (Kmat, MAT_FINAL_ASSEMBLY, ierr)
          call MatAssemblyEnd   (Kmat, MAT_FINAL_ASSEMBLY, ierr)

c         Check on preallocation of matrix

          if(debug) then
            call MatGetInfo(Kmat, MAT_LOCAL, info, ierr)
            nza = info(MAT_INFO_NZ_ALLOCATED)
            nzr = info(MAT_INFO_NZ_USED)
            mal = info(MAT_INFO_MALLOCS)
            write(iow,2000) ' Matrix Information',rank,nza,nzr,mal
            write(  *,2000) ' Matrix Information',rank,nza,nzr,mal
2000  format(/a/' Rank[',i3,']'/
     &          ' Matrix non-zeros allocated =',1p,1d21.12/
     &          ' Matrix non-zeros used      =',1p,1d21.12/
     &          ' Number of mallocs          =',1p,1d21.12)

          endif

c         Stiffness matrix viewer

          if(pfeap_view) then
            call PetscViewerASCIIOpen(PETSC_COMM_WORLD,"stiff.m",K_view,
     &                                ierr)
            call PetscViewerSetFormat(K_view,PETSC_VIEWER_ASCII_MATLAB,
     &                                ierr)
            call MatView             (Kmat,K_view,ierr)
            call PetscViewerDestroy  (K_view,ierr)
          endif

        endif

c       Perform assembly of RHS

        if(flags(4)) then

          call VecValid(rhs, chk, ierr)
          if(chk .eq. PETSC_FALSE) then
            call VecCreate        (PETSC_COMM_WORLD, rhs, ierr)
            call VecSetSizes      (rhs, numpeq, PETSC_DECIDE, ierr)
            call VecSetFromOptions(rhs, ierr)
          endif
          call VecValid(sol, chk, ierr)
          if(chk .eq. PETSC_FALSE) then
            call VecDuplicate     (rhs, sol, ierr)
          endif

          j = 0
          do i = 0, ndf*numpn-1
            k = mr(np(245) + i)
            if( k .gt. 0 ) then
              j = j + 1
              call VecSetValue  (rhs, k-1, b(j), INSERT_VALUES, ierr)
            endif
          end do

          call VecAssemblyBegin (rhs, ierr)
          call VecAssemblyEnd   (rhs, ierr)

c         RHS Vector viewer

          if(pfeap_view) then
            call PetscViewerASCIIOpen(PETSC_COMM_WORLD,"rhs.m",B_view,
     &                                ierr)
            call PetscViewerSetFormat(B_view,PETSC_VIEWER_ASCII_MATLAB,
     &                                ierr)
            call VecView             (rhs,B_view,ierr)
            call PetscViewerDestroy  (B_view,ierr)
          endif

c         Create solution handle solution

          if(pfeap_tang) then
            if(pfeap_dstr) then
              call KSPDestroy       (kspsol,           ierr)
            endif
            pfeap_dstr = .true.
            call KSPCreate        (PETSC_COMM_WORLD, kspsol   ,ierr)
            call KSPSetOperators  (kspsol, Kmat, Kmat,
     &                             DIFFERENT_NONZERO_PATTERN  ,ierr)

c           Set preconditioner type

            call KSPGetPC             (kspsol, pc ,      ierr)
            call PetscOptionsGetString(PETSC_NULL_CHARACTER,'-pc_type',
     &           pctype,chk,ierr)

            if (chk .eq. PETSC_TRUE) then

              if(pctype.eq.'shell') then

                call PCShellSetApply(pc,FeapUserPCApply,ierr)
                call FeapUserPCSetUp(Kmat,kspsol,rhs,sol,ierr)

              else if (pctype.eq.'prometheus') then

c               New style for petsc-2.3.1
                call PCSetCoordinates(pc,ndm,hr(np(43)),ierr)
c               old style for petsc-2.3.0
c               call PCPrometheusSetCoordinates(pc,hr(np(43)),ierr)
              endif

            endif
            pfeap_tang = .false.

c           Set tolerance values for convergence checks

            call KSPSetFromOptions(kspsol,           ierr) ! From command line
          endif
          call KSPGetTolerances (kspsol,p_rtol,p_atol,p_dtol,p_maxit,
     &                           ierr)
          if(pmaxit.le.0) then
            pmaxit = p_maxit
          endif
          call KSPSetTolerances (kspsol,  itol,  atol,  dtol,pmaxit,
     &                           ierr)                   ! From command data

c         Solve equations and destroy handle

          call KSPSolve         (kspsol, rhs, sol, ierr)
          call KSPGetConvergedReason(kspsol,reason,ierr)
          if(reason.gt.0 .and. echo) then
            call KSPGetIterationNumber(kspsol,iters,ierr)
            write(  *,*) ' CONVERGENCE: ',creason(reason),
     &                   ' Iterations = ',iters
            write(iow,*) ' CONVERGENCE: ',creason(reason),
     &                   ' Iterations = ',iters
          elseif(reason.lt.0) then
            write(  *,*) ' NO CONVERGENCE REASON: ',nreason(-reason)
            write(iow,*) ' NO CONVERGENCE REASON: ',nreason(-reason)
          endif
c         call KSPDestroy       (kspsol,           ierr)

c         Set local values in b(*)

          call VecGetOwnershipRange(sol, starti, endi ,ierr)
          call VecGetArray         (sol, sol_array, i_sol, ierr)

          do i = 1,numpeq
            b(i) = sol_array(i_sol + i)
          end do ! i

          call VecRestoreArray(sol, sol_array, i_sol, ierr)

c         Compute max message buffer sizes for ghost communication

          pmaxr = 0
          do i = 0, ntasks-1
            pmaxr = max( pmaxr, mr(np(248)+i) )
          end do
          pmaxr = pmaxs*ndf

          pmaxs = 0
          do i = 0, ntasks-1
            pmaxs = max( pmaxs, mr(np(250)+i) )
          end do
          pmaxs = pmaxs*ndf

c         Perform communication and set of ghost values in solution

          setvar = palloc( 120, 'TEMP0', pmaxr+pmaxs, 2)
          call psetb(b,mr(np(248)),mr(np(249)),mr(np(250)),mr(np(251)),
     &                 mr(np( 31)),ndf, hr(np(120)),hr(np(120)+pmaxr),
     &                 pmaxr,pmaxs)
          setvar = palloc( 120, 'TEMP0',    0, 2)

c         Compute 'energy' for convergence test

          call VecDot(rhs,sol, aengy,ierr)

        endif

c       Purge storage in 'factor'

        if(flags(5)) then
c         write(*,*) 'Rank: ',rank,' USOLVE flag(5) true'
        endif

      endif

      end
