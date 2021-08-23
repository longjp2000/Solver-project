c$Id: piacel.f,v 1.1 2006/11/20 20:32:37 rlt Exp $
      subroutine piacel(ml,dr,a,neq)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Compute starting acceleration for transient problems
c               with diagonal (lumped) mass type arrays

c      Inputs:
c         ml(*)    - Diagonal mass type array
c         dr(*)    - Residual
c         neq      - Number of active equations

c      Outputs:
c         a(*)     - Initial acceleration
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      integer   n,neq
      real*8    ml(*),dr(*),a(*)

      save

c     Compute starting acceleration

      do n = 1,neq
        if(ml(n).ne.0.0d0) then
          a(n) = dr(n)/ml(n)
        else
          a(n) = 0.0d0
        endif
      end do ! n

      end
