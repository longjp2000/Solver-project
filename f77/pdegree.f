c$Id: pdegree.f,v 1.1 2006/11/20 20:34:00 rlt Exp $
      subroutine pdegree(angle, sind,cosd)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Compute sin and cos in terms of degree angle.
c               F77 version

c      Input:
c         angle  - Angle in degrees

c      Outputs:
c         sind   - Sine of angle
c         cosd   - Cosine of angle
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      real*8     angle, sind,cosd

      sind = sin(atan(1.d0)*angle/45.d0)
      cosd = cos(atan(1.d0)*angle/45.d0)

      end
