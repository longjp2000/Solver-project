c$Id: apre3f.f,v 1.1 2006/11/20 20:32:56 rlt Exp $
      subroutine apre3f(d,detf, sig,aa,u, jsw)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose: Finite deformation volumetric elastic model
c                     W(J,be) = U(J)
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      integer   i,j, jsw
      real*8    detf,detfi,u,up,upp,ha,hp,hpp,press,d(*),sig(10),aa(6,6)

      save

c     Divide deviatoric stress and tangent by jacobian

      detfi = 1.d0/detf

      do i = 1,6
        sig(i) = sig(i)*detfi
        do j = 1,6
          aa(i,j) = aa(i,j)*detfi
        end do ! j
      end do ! i

c     Compute pressure and volumetric moduli

      call fengy3(d,detf, u,up,upp,ha,hp,hpp, jsw)

c     Pressure and tangent (not mixed pressure)

      press =  up
      upp   =  upp * detf

      do i = 1,3
        sig(i) = sig(i) + press
      end do ! i

c     Add volumetric correction to aa

      aa(1,1) = aa(1,1) - press + upp
      aa(1,2) = aa(1,2) + press + upp
      aa(1,3) = aa(1,3) + press + upp

      aa(2,1) = aa(2,1) + press + upp
      aa(2,2) = aa(2,2) - press + upp
      aa(2,3) = aa(2,3) + press + upp

      aa(3,1) = aa(3,1) + press + upp
      aa(3,2) = aa(3,2) + press + upp
      aa(3,3) = aa(3,3) - press + upp

      aa(4,4) = aa(4,4) - press
      aa(5,5) = aa(5,5) - press
      aa(6,6) = aa(6,6) - press

      end
