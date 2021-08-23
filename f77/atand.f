c$Id: atand.f,v 1.1 2006/11/20 20:34:00 rlt Exp $
      function atand(x)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Compute arctangent in degrees
c-----[--.----+----.----+----.-----------------------------------------]

      implicit none

      real*8 atand, x

      atand = 45.d0/atan(1.d0)*atan(x)

      end
