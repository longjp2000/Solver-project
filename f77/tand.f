c$Id: tand.f,v 1.1 2006/11/20 20:34:00 rlt Exp $
      function tand(x)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Compute tangent for degrees 'x'
c-----[--.----+----.----+----.-----------------------------------------]

      implicit none

      real*8 tand, x

      tand = tan(atan(1.d0)*x/45.d0)

      end
