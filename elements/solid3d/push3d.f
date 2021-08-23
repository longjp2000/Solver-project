c$Id: push3d.f,v 1.1 2006/11/20 20:32:56 rlt Exp $
      subroutine push3d(a,f, a1)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose: Push forward (or pull back) second rank tensor

c     Inputs:
c       a(6)   - Second rank tensor stored as
c                     | a(1)  a(4)  a(6) |
c                a  = | a(4)  a(2)  a(5) |
c                     | a(6)  a(5)  a(3) |
c       f(3,3) - Deformation gradient

c     Outputs:
c       a1(6)  - Second rank tensor stored same as 'a'
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      integer   i
      real*8    a(6),f(3,3),t(3,3),  a1(6)

c     Compute push forward (or pull back)

      do i = 1,3

c     Compute f*a by rows -> store in t

        t(i,1) = f(i,1)*a(1) + f(i,2)*a(4) + f(i,3)*a(6)
        t(i,2) = f(i,1)*a(4) + f(i,2)*a(2) + f(i,3)*a(5)
        t(i,3) = f(i,1)*a(6) + f(i,2)*a(5) + f(i,3)*a(3)
      end do

c     Compute t*f-T by columns -> store in a1

      a1(1) = t(1,1)*f(1,1) + t(1,2)*f(1,2) + t(1,3)*f(1,3)
      a1(2) = t(2,1)*f(2,1) + t(2,2)*f(2,2) + t(2,3)*f(2,3)
      a1(3) = t(3,1)*f(3,1) + t(3,2)*f(3,2) + t(3,3)*f(3,3)
      a1(4) = t(2,1)*f(1,1) + t(2,2)*f(1,2) + t(2,3)*f(1,3)
      a1(5) = t(3,1)*f(2,1) + t(3,2)*f(2,2) + t(3,3)*f(2,3)
      a1(6) = t(3,1)*f(1,1) + t(3,2)*f(1,2) + t(3,3)*f(1,3)

      end
