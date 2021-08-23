c$Id: sbody3d.f,v 1.1 2006/11/20 20:32:56 rlt Exp $
      subroutine sbody3d(d,xl, r,ndm,ndf ,isw)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Body force computations

c      Inputs:
c        d(*)      - Material set parameters
c        xl(ndm,*) - Element nodal coordinates
c        ndm       - Mesh dimension
c        ndf       - Dofs/node
c        isw       - Option flag

c      Outputs:
c        r(ndf,*)  - Body force nodal vector
c-----[--.----+----.----+----.-----------------------------------------]

      implicit   none

      include   'eldata.h'
      include   'elbody.h'

      integer    ndm, ndf, isw, l,lint,i, ord
      real*8     d(*), xl(ndm,*), r(ndf,*)
      real*8     body(3), sg(4,125), sv(5,16), shp(4,27), xsj

c     Body force computations

      if(isw.eq.15) then
        body(1) = bodyf(1)
        body(2) = bodyf(2)
        body(3) = bodyf(3)
      else
        call sbodyf(d, body)
      endif

c     Set quadrature order

      if(nel.eq.4) then
        ord = 1
        if(nint(d(182)).gt.0) then
          call tint3dn(nel,lint,sv)
        else
          l =  2
          call tint3d (l,lint,sv)
        endif
      elseif(nel.eq.10) then
        ord =  2
        l   =  3
        call tint3d(l,lint,sv)
      elseif(nel.eq.11) then
        ord =  2
        if(nint(d(182)).gt.0) then
          call tint3dn(nel,lint,sv)
        else
          l =  4
          call tint3d(l,lint,sv)
        endif
      else
        ord = 0
        if(nint(d(182)).gt.0) then
          call int3dn(nel,lint,sg)
        else
          l = nint(d(5))
          call int3d(l,lint,sg)
        endif
      endif

c     Compute body loadings

      do l = 1,lint
        if(ord.gt.0) then
          call tetshp(sv(1,l),xl,ndm,ord,xsj,shp)
          xsj = xsj*sv(5,l)
        else
          call shp3d(sg(1,l),xsj,shp,xl,ndm,nel)
          xsj = xsj*sg(4,l)
        endif

c       Compute residual

        do i = 1,nel
          r(1,i)  = r(1,i) + shp(4,i)*body(1)*xsj
          r(2,i)  = r(2,i) + shp(4,i)*body(2)*xsj
          r(3,i)  = r(3,i) + shp(4,i)*body(3)*xsj
        end do ! i
      end do ! l

      end
