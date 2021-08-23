c$Id: uplagm.f,v 1.1 2006/11/20 20:32:37 rlt Exp $
      subroutine uplagm(du,ulagr,lagre,ie,ix)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Accumlate Lagrange multiplier unknowns

c      Inputs:
c         du(*)      - Increment to solution
c         lagre(*)   - Lagrange multiplier equation numbers
c         ie(nie,*)  - Element group control data
c         ix(nen1,*) - Element connection data

c      Outputs:
c         ulagr(*)   - Lagrange multiplier values
c-----[--+---------+---------+---------+---------+---------+---------+-]
      implicit   none

      include   'cdata.h'
      include   'cdat1.h'
      include   'sdata.h'

      integer    i, ma,mm, n,nn
      integer    lagre(*),ie(nie,*),ix(nen1,*)
      real*8     ulagr(*),du(*)

      save

      nn = 0
      do n = 1,numel
        ma = ix(nen1,n)
        if(ie(nie-8,ma).gt.0) then
          mm = lagre(n) - 1
          do i = 1,ie(nie-8,ma)
            ulagr(nn+i) = ulagr(nn+i) + du(mm+i)
          end do ! i
        endif
        nn = nn + ndl
      end do ! n

      end
