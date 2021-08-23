c$Id: pltaxs.f,v 1.1 2006/11/20 20:33:12 rlt Exp $
      subroutine pltaxs(xi,ndm,ct)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Draw vectors on screen for direction of coord. axes

c      Inputs:
c         xi(*)     - Location on origin for axes
c         ndm       - Spatial dimension of mesh
c         ct        - Size factor for plot

c      Outputs:
c         none      - Plot outputs to screen/file
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      integer   ndm,m,n
      real*8    ct,fac1,fac2,fac3, dd(3,3),xx(3,5),xi(ndm)

      save

c     Perspective projecton of axes

      do m = 1,3
        do n = 1,3
          dd(n,m) = 0.0d0
        end do ! n
      end do ! m
      do m = 1,ndm
        dd(m,m) = ct
      end do ! m

c     Compute plot location for axes

      do m = 1,ndm
        call pzero(xx,15)
        do n = 1,ndm
          fac1 = dd(n,m)
          xx(n,1) = xi(n)
          xx(n,2) = xx(n,1) + fac1
          xx(n,5) = xx(n,2)
        end do ! n
        fac1 = dd(1,m)*0.1d0
        fac2 = dd(2,m)*0.1d0
        fac3 = dd(3,m)*0.1d0
        xx(1,3) = xx(1,2) - 3.d0*fac1 -  fac2 - fac3
        xx(2,3) = xx(2,2) - 3.d0*fac2 +  fac1 + fac3
        xx(3,3) = xx(3,2) - 3.d0*fac3 +  fac1 + fac2
        xx(1,4) = xx(1,2) - 3.d0*fac1 +  fac2 + fac3
        xx(2,4) = xx(2,2) - 3.d0*fac2 -  fac1 - fac3
        xx(3,4) = xx(3,2) - 3.d0*fac3 -  fac1 - fac2

c       Plot vector

        call plotl(xx(1,1),xx(2,1),xx(3,1),3)
        do n = 2,5
          call plotl(xx(1,n),xx(2,n),xx(3,n),2)
        end do ! n
        call plotl(xx(1,2),xx(2,2),xx(3,2),3)

c       Add label

        call plabl(m)

      end do ! M

      end
