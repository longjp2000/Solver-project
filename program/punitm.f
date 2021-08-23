c$Id: punitm.f,v 1.1 2006/11/20 20:32:37 rlt Exp $
      subroutine punitm(aa,nn,unit)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose: Scale 'aa' array by 'unit'

c     Inputs:

c     Outputs:
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      integer    n,nn
      real*8     aa(nn), unit

      save

      do n = 1,nn
        aa(n) = aa(n)*unit
      end do ! n

      end
