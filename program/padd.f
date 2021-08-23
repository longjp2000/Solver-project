c$Id: padd.f,v 1.1 2006/11/20 20:32:37 rlt Exp $
      function padd(val)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose: Accumulate a parameter

c     Input:
c        val   - Value to accumulate

c     Output:
c        padd  - Accumulated or zeroed value
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      real*8    padd, val, xval

      save

      data      xval /0.0d0/

c     Look at parameter

      if(val.eq.0.0d0) then
        xval = 0.0d0
      else
        xval = xval + val
      endif

      padd = xval

      end
