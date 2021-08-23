c$Id: shp3dc.f,v 1.1 2006/11/20 20:32:37 rlt Exp $
      subroutine shp3dc(ss, shp)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose:

c      Inputs:
c         ss(*)    - Gauss points

c      Outputs:
c         shp(4,*) - Cubic shape functions
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      integer    i, ir(64), is(64), it(64)
      real*8     nr(4),dr(4), ns(4),ds(4), nt(4),dt(4), ss(3), shp(4,*)

      save

      data       ir /1,2,2,1, 1,2,2,1, 1,2,2,1, 1,2,2,1,
     &               3,4,2,2,4,3,1,1,3,4,4,3,
     &               3,4,2,2,4,3,1,1,3,4,4,3,
     &               3,4,2,2,4,3,1,1,3,4,4,3,
     &               3,4,2,2,4,3,1,1,3,4,4,3/

      data       is /1,1,2,2, 1,1,2,2, 1,1,2,2, 1,1,2,2,
     &               1,1,3,4,2,2,4,3,3,3,4,4,
     &               1,1,3,4,2,2,4,3,3,3,4,4,
     &               1,1,3,4,2,2,4,3,3,3,4,4,
     &               1,1,3,4,2,2,4,3,3,3,4,4/

      data       it /1,1,1,1, 2,2,2,2, 3,3,3,3, 4,4,4,4,
     &               1,1,1,1,1,1,1,1,1,1,1,1,
     &               2,2,2,2,2,2,2,2,2,2,2,2,
     &               3,3,3,3,3,3,3,3,3,3,3,3,
     &               4,4,4,4,4,4,4,4,4,4,4,4/

c     Set 1-d shape functions for each local direction

      nr(3) =  9.0d0*(1.d0 - ss(1)*ss(1))*(1.d0 - 3.d0*ss(1))/16.d0
      nr(4) =  9.0d0*(1.d0 - ss(1)*ss(1))*(1.d0 + 3.d0*ss(1))/16.d0
      nr(1) =  0.5d0 * (1.d0 - ss(1)) - (2.d0*nr(3) + nr(4))/3.d0
      nr(2) =  0.5d0 * (1.d0 + ss(1)) - (2.d0*nr(4) + nr(3))/3.d0

      ns(3) =  9.0d0*(1.d0 - ss(2)*ss(2))*(1.d0 - 3.d0*ss(2))/16.d0
      ns(4) =  9.0d0*(1.d0 - ss(2)*ss(2))*(1.d0 + 3.d0*ss(2))/16.d0
      ns(1) =  0.5d0 * (1.d0 - ss(2)) - (2.d0*ns(3) + ns(4))/3.d0
      ns(2) =  0.5d0 * (1.d0 + ss(2)) - (2.d0*ns(4) + ns(3))/3.d0

      nt(3) =  9.0d0*(1.d0 - ss(3)*ss(3))*(1.d0 - 3.d0*ss(3))/16.d0
      nt(4) =  9.0d0*(1.d0 - ss(3)*ss(3))*(1.d0 + 3.d0*ss(3))/16.d0
      nt(1) =  0.5d0 * (1.d0 - ss(3)) - (2.d0*nt(3) + nt(4))/3.d0
      nt(2) =  0.5d0 * (1.d0 + ss(3)) - (2.d0*nt(4) + nt(3))/3.d0

      dr(3) =  9.0d0*( 9.d0*ss(1)*ss(1) - 2.d0*ss(1) - 3.d0)/16.d0
      dr(4) =  9.0d0*(-9.d0*ss(1)*ss(1) - 2.d0*ss(1) + 3.d0)/16.d0
      dr(1) = -0.5d0 - (2.d0*dr(3) + dr(4))/3.d0
      dr(2) =  0.5d0 - (2.d0*dr(4) + dr(3))/3.d0

      ds(3) =  9.0d0*( 9.d0*ss(2)*ss(2) - 2.d0*ss(2) - 3.d0)/16.d0
      ds(4) =  9.0d0*(-9.d0*ss(2)*ss(2) - 2.d0*ss(2) + 3.d0)/16.d0
      ds(1) = -0.5d0 - (2.d0*ds(3) + ds(4))/3.d0
      ds(2) =  0.5d0 - (2.d0*ds(4) + ds(3))/3.d0

      dt(3) =  9.0d0*( 9.d0*ss(3)*ss(3) - 2.d0*ss(3) - 3.d0)/16.d0
      dt(4) =  9.0d0*(-9.d0*ss(3)*ss(3) - 2.d0*ss(3) + 3.d0)/16.d0
      dt(1) = -0.5d0 - (2.d0*dt(3) + dt(4))/3.d0
      dt(2) =  0.5d0 - (2.d0*dt(4) + dt(3))/3.d0

c     Set local 3-d shape functions

      do i = 1,64
        shp(1,i) = dr(ir(i)) * ns(is(i)) * nt(it(i))
        shp(2,i) = nr(ir(i)) * ds(is(i)) * nt(it(i))
        shp(3,i) = nr(ir(i)) * ns(is(i)) * dt(it(i))
        shp(4,i) = nr(ir(i)) * ns(is(i)) * nt(it(i))
      end do ! i

      end
