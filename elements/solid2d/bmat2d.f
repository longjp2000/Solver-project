c$Id: bmat2d.f,v 1.1 2006/11/20 20:32:56 rlt Exp $
      subroutine bmat2d(c,r,shp,g,bbar)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Sets B-bar matrix for 2-d problems

c      Inputs:
c         c         - Constant for plane = 0; for axisymm = 1
c         r         - Radius for axisymmetrix (= 1 for plane)
c         shp(3)    - Shape function and derivatives
c         g(2)      - b-bar integrals

c      Outputs:
c         bbar(4,2) - B-bar matrix for a node
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      real*8    c,r,bb1,bb2,sh3, shp(3),g(2),bbar(4,2)

c     Mixed modification to form B-bar

      sh3 = c*shp(3)/r
      bb1 = (g(1) - shp(1) - sh3)*0.3333333333333333d0
      bb2 = (g(2) - shp(2)      )*0.3333333333333333d0

c     B-bar matrix for plane and axisymmetric problems

      bbar(1,1) = bb1 + shp(1)
      bbar(2,1) = bb1
      bbar(3,1) = bb1 + sh3
      bbar(4,1) = shp(2)
      bbar(1,2) = bb2
      bbar(2,2) = bb2 + shp(2)
      bbar(3,2) = bb2
      bbar(4,2) = shp(1)

      end
