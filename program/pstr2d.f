c$Id: pstr2d.f,v 1.1 2006/11/20 20:32:37 rlt Exp $
      subroutine pstr2d(sig,pp)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c     1. Use 'pi' from 'pconstant.h'                        14/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Compute principal stresses for 2-d problems.

c      Input:
c         sig(1) - Stresses in order: sig-xx, sig-yy, sig-zz, sig-xy
c      Output:
c         pp(1)  - Principal stresses in order: sig-1, sig-2, angle
c                  (degrees): sig-xx to sig-1
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'pconstant.h'

      real*8    sig(4), pp(3), xi1, xi2, rho

      save

      xi1     = (sig(1) + sig(2))*0.5d0
      xi2     = (sig(1) - sig(2))*0.5d0
      rho     = sqrt(xi2*xi2 + sig(4)*sig(4))
      pp(1)   = xi1 + rho
      pp(2)   = xi1 - rho
      if(xi2.ne.0.0d0) then
        pp(3) = 90.0d0*atan2(sig(4),xi2)/pi
      else
        pp(3) = 45.0d0
      endif

      end
