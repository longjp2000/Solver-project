c$Id: pside1.f,v 1.1 2006/11/20 20:32:37 rlt Exp $
      subroutine pside1(nr,xs,tr,side,is,ns,ndm,shp,rt, x, styp)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose: Construct one dimensional side interpolation for coords.

c     Inputs:
c       nr        - Number of increments on side
c       xs(3,*)   - Nodal values of interpolation function
c       side      - side number (check sign)
c       is(*)     - List of side nodes
c       ns        - Order of Lagrange polynomial for side
c       ndm       - Spatial dimension of mesh
c       shp(*)    - Shape functions for nodal values
c       rt(3,*)   - Temporary storage for polar coordinates
c       styp      - Type of edge: 0 = Lagrange interpolation
c                                 1 = Radius/angle interpolation
c                                 2 = Segmental lines
c                                 3 = Eliptical interpolation

c     Outputs:
c       x(ndm,ip) - Coordinates of points
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'cdata.h'
      include  'pointer.h'
      include  'comblk.h'

      integer   k,n, nr,ns,side,ndm, styp, is(*)
      real*8    xs(3,*), shp(*), rt(3,*), x(ndm,*), xx(3), tr(3,4)

      save

c     Lagrange interpolation

      if(styp.eq.0) then

        call interp1(nr,xs,side,is,ns,ndm,shp, x)

c     R/theta interpolation

      elseif(styp.eq.1) then

        call arcint1(nr,xs,side,is,ns,ndm,shp, x,rt)

c     Straight segment interpolation

      elseif(styp.eq.2) then

        call segint1(xs,side,is,ns,ndm,shp, x)

c     Eliptical interpolation

      elseif(styp.eq.3) then

        call elpint1(nr,xs,side,is,ns,ndm,shp, x,rt)

c     Bernstein interpolation

      elseif(styp.eq.4) then

        call bernst1(nr,xs,side,is,ns,ndm, x)

      endif

c     Transform to current frame

      do n = 1,nr+1
        do k = 1,ndm
          xx(k) = x(k,n)
        end do ! k
        do k = 1,ndm
          x(k,n) = tr(k,4)+tr(k,1)*xx(1)+tr(k,2)*xx(2)+tr(k,3)*xx(3)
        end do ! k
      end do ! n

      end
