c$Id: uasble.F,v 1.1 2006/11/21 16:44:39 rlt Exp $
      subroutine uasble(s,p,ld,ns,afl,bfl,b)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose: User to assemble current matrix/vector
c              User is responsible for creating the storage locations
c              for global arrays.

c     Inputs:
c       s(ns,ns,2)- element matrix
c       p(ns)     - element vector
c       ld(ns,2)  - local/global active equation numbers
c                   1 = rows; 2 = columns.
c       ns        - size of arrays
c       afl       - Assemble s(ns,ns) into global storage
c       bfl       - Assemble p(ns)    into global storage

c     Outputs:
c       b(*)      - Assembled RHS vector
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

#     include   "include/finclude/petsc.h"
#     include   "include/finclude/petscvec.h"
#     include   "include/finclude/petscmat.h"
#     include   "pfeapc.h"
#     include   "pfeapd.h"

#     include   "compas.h"
#     include   "rdat1.h"
#     include   "pfeapb.h"

      PetscInt   ierr

      logical    afl,bfl
      integer    i, ns, n, m, j
      integer    ld(ns,*)
      real*8     s(ns,ns,2),p(ns),b(*),val

c     Assemble matrix

      if(afl) then
        do m = 1,nar
          if(ld(m,1).gt.0) then
            i = ld(m,1) - 1
            do n = ld(m,3),ld(m,4)
              j = ld(n,2) - 1
              val = s(m,n,1)
              call MatSetValue( Kmat, i, j, val, ADD_VALUES, ierr )
            end do ! n
          endif
        end do ! m

c       Set flag indicating new tangent formed

        pfeap_tang = .true.

      endif

c     Assemble vector

      if(bfl) then
        do n = 1, nar
          if(ld(n,6).gt.0) then
            b(ld(n,6))  = b(ld(n,6))  + p(n)
          endif
        end do ! n
        if(compre) then  ! Compute residual norm for convergence check
          do n = 1,nar
            rnorm1 = rnorm1 + abs(p(n))
            rnormn = rnormn + 1.d0
          end do ! i
        endif
      endif

      end
