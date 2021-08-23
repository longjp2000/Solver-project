c$Id: smodify.f,v 1.1 2006/11/21 16:44:39 rlt Exp $
      subroutine smodify(ld, p, s, nst, afl)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Modify tangent for imposed displacements, set RHS = 0
c
c      Inputs:
c        ld(*)   - Boundary condition indicator
c        p(*)    - Unmodified element residual
c        s(*,*)  - Element tangent array
c        afl     - Tangent flag

c      Outputs:
c        p(*)    - Modified residual
c        s(*,*)  - Modified tangent array
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include   'pfeapb.h'

      logical    afl
      integer    n,m, nst
      integer    ld(nst,*)
      real*8     p(nst), s(nst,*)

      if(afl) then               ! Modify tangent and residual
        do n = 1,nst
          if(ld(n,5).ne.0) then  ! inactive equation
            do m = 1,nst
              s(n,m) = 0.0d0
              s(m,n) = 0.0d0
            end do ! m
            s(n,n) = 1.0d0
            p(n)   = 0.0d0
          endif
        end do ! n
      else                       ! Modify residual
        do n = 1,nst
          if(ld(n,5).ne.0) then  ! inactive equation
            p(n)   = 0.0d0
          endif
        end do ! n
      endif

      end
