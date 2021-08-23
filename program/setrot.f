c$Id: setrot.f,v 1.1 2006/11/20 20:32:37 rlt Exp $
      subroutine setrot ( x , mo, xlm , thk , numnp )

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Obtain directors from the original mesh x(1-3,*) and
c               second mesh xlm(1-3,3,*), normalize and set trans-
c               formation using the exponential map. The norm is
c               multiplied by 2 and stored in thk(*).

c      Inputs:
c         x(3,*)  - Nodal coordinates for mesh
c         mo(*)   - Rotational update type
c         numnp   - Number of nodes in mesh

c      Outputs:
c         xlm(*)  - Normalized nodal directors
c         thk(*)  - Nodal thicknesses
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      integer   i,n, numnp, mo(*)
      real*8    xn, x1, x2, x3, fac, x(3,*),xlm(9,6,*),thk(*)

      save

      do n = 1 , numnp

        if(mo(n).eq.-1 .or. mo(n).eq.-5) then

c         Compute Fiber:

          x1         = xlm(7,1,n) - x(1,n)
          x2         = xlm(8,1,n) - x(2,n)
          x3         = xlm(9,1,n) - x(3,n)

c         Normalize:

          xn         = sqrt ( x1**2 + x2**2 + x3**2 )
          if (xn.gt.1.d-6) then
            xlm(7,1,n) = x1 / xn
            xlm(8,1,n) = x2 / xn
            xlm(9,1,n) = x3 / xn
            thk(n)     = xn * 2.d0
          else
            xlm(7,1,n) = 0.d0
            xlm(8,1,n) = 0.d0
            xlm(9,1,n) = 1.d0
            thk(n)       = 1.d0
          endif

c         Assemble [Lambda]:

          if (xlm(9,1,n).gt.0.d0) then
            fac          =  1.d0 / ( 1.d0 + xlm(9,1,n) )
            xlm(1,1,n) =  xlm(9,1,n) + fac*xlm(8,1,n)*xlm(8,1,n)
            xlm(4,1,n) =             - fac*xlm(8,1,n)*xlm(7,1,n)
            xlm(2,1,n) =             - fac*xlm(7,1,n)*xlm(8,1,n)
            xlm(5,1,n) =  xlm(9,1,n) + fac*xlm(7,1,n)*xlm(7,1,n)
            xlm(3,1,n) = -xlm(7,1,n)
            xlm(6,1,n) = -xlm(8,1,n)
          else
            fac          =  1.d0 / ( 1.d0 - xlm(9,1,n) )
            xlm(1,1,n) = -xlm(9,1,n) + fac*xlm(8,1,n)*xlm(8,1,n)
            xlm(4,1,n) =               fac*xlm(8,1,n)*xlm(7,1,n)
            xlm(2,1,n) =             - fac*xlm(7,1,n)*xlm(8,1,n)
            xlm(5,1,n) =  xlm(9,1,n) - fac*xlm(7,1,n)*xlm(7,1,n)
            xlm(3,1,n) =  xlm(7,1,n)
            xlm(6,1,n) = -xlm(8,1,n)
          endif
          do i = 1,9
            xlm(i,2,n) = xlm(i,1,n)
            xlm(i,3,n) = xlm(i,1,n)
            xlm(i,6,n) = xlm(i,1,n)
          end do ! i
        endif
      end do ! n

      end
