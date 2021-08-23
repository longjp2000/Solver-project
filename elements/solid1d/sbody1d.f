c$Id: sbody1d.f,v 1.1 2006/11/20 20:32:56 rlt Exp $
      subroutine sbody1d(d,xl, r,ndm,ndf ,isw)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Body force computation

c      Inputs:

c      Outputs:
c-----[--.----+----.----+----.-----------------------------------------]

      implicit   none

      include   'eldata.h'
      include   'elbody.h'

      integer    ndm, ndf, isw, l,lint,i
      real*8     d(*), xl(ndm,*), r(ndf,*), sg(2,4), shp(2,4), bf(3)
      real*8     rr, xsj

c     Body force computations

      if(isw.eq.15) then
        bf(1) = bodyf(1)
      else
        call sbodyf(d, bf)
      endif

c     Set quadrature order

      if(nel.eq.2) then
        lint = 2
      else
        lint = 3
      endif
      if(nint(d(182)).eq.1) then
        call int1dn(lint,sg)
      else
        call int1d(lint,sg)
      endif

c     Compute body loadings

      do l = 1,lint
        call shp1d(sg(1,l),xl,shp,ndm,nel,xsj)
        xsj = xsj*sg(2,l)

        if(nint(d(16)).eq.3) then
          rr = 0.0d0
          do i = 1,nel
            rr = rr + shp(2,i)*xl(1,i)
          end do ! i
          xsj = xsj*rr
        endif

c       Compute residual

        do i = 1,nel
          r(1,i)  = r(1,i) + shp(2,i)*bf(1)*xsj
        end do ! j
      end do ! l

      end
