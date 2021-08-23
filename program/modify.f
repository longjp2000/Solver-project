c$Id: modify.f,v 1.1 2006/11/20 20:32:37 rlt Exp $
      subroutine modify(ld,p,s,dul,nst)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Modify element residual for effects of specified
c               boundary values.

c               p(i) = p(i) - s(i,j)*dul(j)

c      Inputs:
c         ld(*)  - Array with negative entries where boundary
c                  solution to be imposed
c         p(*)   - Unmodified residual from element
c         s(*,*) - Element tangent array
c         dul(*) - Value of specified solution increments
c         nst    - Dimension of element arrays

c      Outputs:
c         p(*)   - Residual modified for effect of increments
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      integer   nst,i,j
      integer   ld(nst)
      real*8    p(nst),s(nst,nst),dul(nst)

      save

c     Loop over columns and search for boundary terms

      do j = 1,nst
        if(ld(j).lt.0) then

c         Loop over rows to modify active equations

          do i = 1,nst
            if(ld(i).gt.0) then
              p(i) = p(i) - s(i,j)*dul(j)
            endif
          end do ! i

        endif
      end do ! j

      end
