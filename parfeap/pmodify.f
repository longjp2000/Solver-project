c$Id: pmodify.f,v 1.1 2006/11/21 16:44:39 rlt Exp $
      subroutine pmodify(ld, p, s, ub, nst)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Modify RHS for imposed displacements
c
c      Inputs:
c        ld(*,2) - Equations for rows and columns
c        p(*)    - Unmodified element residual
c        s(*,*)  - Element tangent array

c      Outputs:
c        p(*)    - Modified residual for effects of non-zero displacement
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include   'pfeapb.h'

      integer    n,m, nst
      integer    ld(nst,*)
      real*8     p(nst), s(nst,*), ub(nst)

      if(pfeap_blk) then
        do n = 1,nst
          if(ld(n,5).eq.0) then  ! active equation
            do m = 1,nst
              if(ub(m).ne.0.0d0) then
                 p(n) = p(n) - s(n,m)*ub(m)
              endif
            end do ! m
          endif
        end do ! n
        do n = 1,nst
          if(ld(n,5).ne.0) then  ! inactive equation
            do m = 1,nst
              s(n,m) = 0.0d0
              s(m,n) = 0.0d0
            end do ! m
            s(n,n) = 1.0d0
            p(n)   = ub(n)
          endif
        end do ! n
      else
        do n = 1,nst
          if(ld(n,1).gt.0) then  ! active equation
            do m = 1,nst
              if(ld(m,2).le.0) then
                 p(n) = p(n) - s(n,m)*ub(m)
              endif
            end do ! m
          endif
        end do ! n
      endif

      end
