c$Id: sbody2d.f,v 1.1 2006/11/20 20:32:56 rlt Exp $
      subroutine sbody2d(d,xl,ix, r,ndm,ndf ,isw)

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

      logical    quad
      integer    ndm, ndf, isw, l,lint,i, ix(*)
      real*8     d(*), xl(ndm,*), r(ndf,*), body(3), sg(3,25), shp(3,9)
      real*8     el(4,7), rr, xsj

c     Body force computations

      if(isw.eq.15) then
        body(1) = bodyf(1)
        body(2) = bodyf(2)
      else
        call sbodyf(d, body)
      endif

c     Set quadrature order

      if(nel.eq.3) then
        if(d(182).gt.0.0d0) then
          call tint2dn(nel,lint,el)
        else
          l = 1
        call tint2d (l,lint,el)
        endif
        quad = .false.
      elseif(nel.eq.6 .or. nel.eq.7) then
        if(d(182).gt.0.0d0) then
          call tint2dn(nel,lint,el)
        else
          l = 7
        call tint2d (l,lint,el)
        endif
        quad = .false.
      else
        quad = .true.
        if(nint(d(182)).gt.0) then
          call int2dn(nel,lint,sg)
        else
          l = nint(d(5))
          call int2d (l,lint,sg)
        endif
      endif

c     Compute body loadings

      do l = 1,lint
        if(quad) then
          call shp2d(sg(1,l),xl,shp,xsj,ndm,nel,ix,.false.)
          xsj = xsj*sg(3,l)
        else
          call trishp(el(1,l),xl,ndm,nel-4,xsj,shp)
          xsj = xsj*el(4,l)
        endif

        if(nint(d(16)).eq.3) then
          rr = 0.0d0
          do i = 1,nel
            rr = rr + shp(3,i)*xl(1,i)
          end do ! i
          xsj = xsj*rr
        endif

c       Compute residual

        do i = 1,nel
          r(1,i)  = r(1,i) + shp(3,i)*body(1)*xsj
          r(2,i)  = r(2,i) + shp(3,i)*body(2)*xsj
        end do ! j
      end do ! l

      end
