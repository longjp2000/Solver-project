c$Id: mnrv3f.f,v 1.1 2006/11/20 20:32:56 rlt Exp $
      subroutine mnrv3f(d,detf,bb, sig,aa,xlamd,ha,estore)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose: Finite deformation elasticity Mooney-Rivlin model

c     Inputs:
c         d(100)     Material parameters
c         d(21)      Bulk  modulus
c         d(22)      Shear modulus
c         detf       Jacobian determinant at t_n+1
c         bb(6)      Left Cauchy-Green tensor

c     Outputs:
c         sig(6)     CAUCHY stress tensor
c         aa(6,6)    CAUCHY (spatial) elastic moduli
c         estore     Stored energy density
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      integer   i,j,k,l, ii,jj, jsw, in(3,3)
      real*8    detf, press,u,up,upp, mu, mu2, ha,hp,hpp, estore
      real*8    xlamd,logj, nu,nu2,nu4, trb, trb2
      real*8    d(*),sig(6),aa(6,6),bb(6),b2(6)

      data      in / 1, 4, 6, 4, 2, 5, 6, 5, 3 /

c     Compute pressure and its derivative

      jsw = nint(d(170))
      call fengy3(d,detf,u,up,upp, ha,hp,hpp, jsw)
      press =  up  + xlamd * hp
      upp   = (upp + xlamd * hpp) * detf

c     Set CAUCHY stresses and elastic tangent moduli

      nu  =  d(22)*d(23)/detf
      mu  =  d(22)/detf - nu
      nu2 =  nu + nu
      mu2 =  mu + mu

c     First invariant term

      do i = 1,3
        sig(i  )    = mu*bb(i) - mu + press
        sig(i+3)    = mu*bb(i+3)
      end do ! i

c     Second invariant term

      nu2   = nu  + nu
      nu4   = nu2 + nu2
      b2(1) = bb(1)*bb(1) + bb(4)*bb(4) + bb(6)*bb(6)
      b2(2) = bb(4)*bb(4) + bb(2)*bb(2) + bb(5)*bb(5)
      b2(3) = bb(6)*bb(6) + bb(5)*bb(5) + bb(3)*bb(3)
      b2(4) = bb(1)*bb(4) + bb(4)*bb(2) + bb(6)*bb(5)
      b2(5) = bb(4)*bb(6) + bb(2)*bb(5) + bb(5)*bb(3)
      b2(6) = bb(6)*bb(1) + bb(5)*bb(4) + bb(3)*bb(6)
      trb   = bb(1) + bb(2) + bb(3)
      do i = 1,3
        sig(i  ) = sig(i  ) + nu*(trb*bb(i)   - b2(i  )) - nu2
        sig(i+3) = sig(i+3) + nu*(trb*bb(i+3) - b2(i+3))
      end do ! i

c     Tangent terms from Kronnecker delta

      do i = 1,3
        aa(i  ,i  ) = mu2 + nu4 - press + upp
        aa(i+3,i+3) = mu  + nu2 - press
      end do ! i

c     Tangent terms from second invariant

      do l = 1,3
        do k = l,3
          jj = in(k,l)
          do j = 1,3
            do i = j,3
              ii = in(i,j)
              aa(ii,jj) = aa(ii,jj) + nu2*(bb(ii)*bb(jj)
     &                  - 0.5d0*(bb(in(i,k))*bb(in(j,l))
     &                         + bb(in(i,l))*bb(in(j,k))))
            end do ! i
          end do ! j
        end do ! k
      end do ! l

c     Add volumetric correction to aa

      upp     = press   + upp
      aa(1,2) = aa(1,2) + upp
      aa(2,1) = aa(1,2)
      aa(1,3) = aa(1,3) + upp
      aa(3,1) = aa(1,3)
      aa(2,3) = aa(2,3) + upp
      aa(3,2) = aa(2,3)

c     Compute stored energy density

      logj   =  log(abs(detf))
      trb2   =  b2(1)*b2(1) + b2(2)*b2(2) + b2(3)*b2(3)
     &       + (b2(4)*b2(4) + b2(5)*b2(5) + b2(6)*b2(6))*2.d0
      estore =  u + 0.5d0*d(22)*((1.d0 - d(23))*(trb-3.0d0-2.0d0*logj)
     &            + d(23)*(0.5d0*(trb*trb - trb2)-3.0d0-4.0d0*logj))

      end
