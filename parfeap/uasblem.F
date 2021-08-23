c$Id: uasblem.F,v 1.1 2006/11/21 16:44:39 rlt Exp $
      subroutine uasblem(s,p,ld,ns,afl,bfl, isw)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose: User function to assemble current mass matrix/vector

c     Inputs:
c       s(ns,ns)  - element consistent mass (matrix)
c       p(ns)     - element lumped mass     (vector)
c       ld(ns,2)  - local/global active equation numbers
c                   1 = rows; 2 = columns.
c       ns        - size of arrays
c       afl       - Assemble s(ns,ns) into global storage
c       bfl       - Assemble p(ns)    into global storage
c       isw       - Switch: (1) Send values to Mmat; (2) Assemble

c     Outputs:
c                 - Assembled Mass matrix/vector
c-----[--.----+----.----+----.-----------------------------------------]
      implicit    none

#     include    "include/finclude/petsc.h"
#     include    "include/finclude/petscvec.h"
#     include    "include/finclude/petscmat.h"
#     include    "include/finclude/petscviewer.h"

#     include    "pfeapc.h"
#     include    "pfeapd.h"

#     include    "debugs.h"
#     include    "compas.h"
#     include    "fdata.h"
#     include    "rdat1.h"
#     include    "pfeapb.h"

      PetscInt    ierr
      PetscViewer M_view

      logical     afl,bfl
      integer     isw, i, ns, n, m, j
      integer     ld(ns,*)
      real*8      s(ns,ns),p(ns),val

c     Assemble consistent mass matrix

      if(isw.eq.1) then
        if(afl) then
          do m = 1,nar
            if(ld(m,5).eq.0) then       ! Assemble active rows only
              i = ld(m,1) - 1
              do n = ld(m,3),ld(m,4)
                if(ld(n,5).eq.0) then   ! Assemble active cols only
                  j = ld(n,2) - 1
                  val = s(m,n)
                  if(val.ne.0.0d0) then
                    call MatSetValue( Mmat, i, j, val, ADD_VALUES, ierr)
                  endif
                endif
              end do ! m
            endif
          end do ! n

        endif

c       Assemble lumped mass vector

        if(bfl) then
          val = 0.0d0
          do n = 1, nar
            if(ld(n,5).eq.0) then
              i = ld(n,1) - 1
              call VecSetValue( Mdiag, i, p(n), ADD_VALUES, ierr)
            else
              i = ld(n,1) - 1
              call VecSetValue( Mdiag, i, val , ADD_VALUES, ierr)
            endif
          end do ! n
        endif

c     Assemble mass for parallel form

      else
        if(fl(1)) then
          call MatAssemblyBegin (Mmat, MAT_FINAL_ASSEMBLY, ierr)
          call MatAssemblyEnd   (Mmat, MAT_FINAL_ASSEMBLY, ierr)
        else
          call VecAssemblyBegin (Mdiag, ierr)
          call VecAssemblyEnd   (Mdiag, ierr)
        endif

        if(pfeap_view) then
          call PetscViewerASCIIOpen(PETSC_COMM_WORLD,"mass.m",M_view,
     &                              ierr)
          call PetscViewerSetFormat(M_view,PETSC_VIEWER_ASCII_MATLAB,
     &                              ierr)
          if(fl(1)) then
            call MatView           (Mmat,M_view,ierr)
          else
            call VecView           (Mdiag,M_view,ierr)
          endif
          call PetscViewerDestroy  (M_view , ierr)
        endif

      endif

      end
