c$Id: cmprot.f,v 1.1 2006/11/20 20:32:37 rlt Exp $
      subroutine cmprot(mo,lam, numnp, lmr)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Set rotational array for finite rotation transformations

c      Inputs:

c      Outputs:
c-----[--.----+----.----+----.-----------------------------------------]
      implicit none

      integer  numnp, lmr, n,i
      integer  mo(numnp,2)
      real*8   lam(54,numnp)

      lmr = 0
      do n = 1,numnp
        if(mo(n,1).ne.0) then
          lmr = lmr + 1
          mo(n,2) = lmr
          if(lmr.ne.n) then
            do i = 1,54
              lam(i,lmr) = lam(i,n)
            end do ! i
          endif
        endif
      end do ! n

      end
