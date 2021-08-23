c$Id: pscalc.f,v 1.1 2006/11/20 20:32:37 rlt Exp $
      subroutine pscalc(ir,jc,ad,al,au,da,neq, cfr)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose: Scaling of sparse matrix to have unit diagonals
c              Diagonals in separate array

c     Inputs:

c     Outputs:
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      logical    cfr
      integer    neq, n, j, ir(*),jc(*)
      real*8     ad(*),al(*),au(*),da(*)

      save

c     Normalize diagonal

      do n = 1,neq
        if(ad(n).ne.0.0d0) then
          da(n) = 1.d0/sqrt(abs(ad(n)))
          ad(n) = sign(1.d0,ad(n))
        else
          da(n) = 1.d0
        endif
      end do ! n

c     Scale AU

      do n = 2,neq
        do j = jc(n-1)+1,jc(n)
          au(j) = au(j)*da(n)*da(ir(j))
        end do ! j
      end do ! n

c     Scale AL

      if(cfr) then
        do n = 1,neq
          do j = jc(n-1)+1,jc(n)
            al(j) = al(j)*da(n)*da(ir(j))
          end do ! j
        end do ! n
      endif

      end
