c$Id: wder3f.f,v 1.1 2006/11/20 20:32:56 rlt Exp $
      subroutine wder3f(d,bpr, tautil,atilp)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c     1. Use constants from 'pconstant.h'                   14/11/2006
c-----[--.----+----.----+----.-----------------------------------------]

c     Stored Energy Function in Principal Stretches
c          W       = w(xlam_1) + w(xlam_2) + w(xlam_3)

c          xlam _i = (Je)^(-1/3)*lambda_i (deviatoric stretches)

c     Ogden strain energy functions in terms of Seth strains
c          w_i  = c-alpha/(m-alpha)^2*(xlam_i**m-alpha - 1.) -> wengy

c     Inputs:
c          d(22) = first  ci
c          d(23) = first  ni
c          d(24) = second ci
c          d(25) = second ni
c          d(26) = third  ci
c          d(27) = third  ni
c          d(28) = no. terms in series expansion of W(xlam) = 3, max

c          bpr(3) principal values of left Cauchy-Green tensor

c     Outputs:
c          tautil(3)   Principal deviatoric Kirchoff stresses
c          atilp(6,6)  Deviatoric tangent matrix in principal basis
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'cdamag.h'
      include  'pconstant.h'

      integer   nn, i, j, k, n
      real*8    xjthrd, c1, c2, c3, c4, tol
      real*8    d(*),bpr(3),xlamt(3),xlam2(3),tautil(3)
      real*8    atilp(6,6),dwtil(4,3),dtdl(3,3)

      save

      data      tol    /  1.0d-06 /

c     Initialize arrays

      call pzero(dwtil,12)
      call pzero(dtdl,9)
      call pzero(atilp,36)

c     Compute J**(-1/3)

      
      xjthrd = (bpr(1)*bpr(2)*bpr(3))**(-one6)

c     Compute modified (Flory) stretches

      do i = 1,3
        xlamt(i) = xjthrd*sqrt(bpr(i))
        xlam2(i) = xlamt(i)*xlamt(i)
      end do ! i

c     Compute strain energy function derivatives

      wengy = 0.0d0
      nn    = d(28)

      do n = 1,nn
        c1   = d(20+2*n)
        c2   = d(20+2*n+1)
        do i = 1,3
          wengy      = wengy      + (c1/c2)*(xlamt(i)**c2 - 1.d0)
          dwtil(1,i) = dwtil(1,i) + c1*xlamt(i)**(c2-1.d0)
          dwtil(2,i) = dwtil(2,i) + (c2-1.d0)*c1*xlamt(i)**(c2-2.d0)
          dwtil(3,i) = dwtil(3,i)
     &                  + (c2-2.d0)*(c2-1.d0)*c1*xlamt(i)**(c2-3.d0)
          dwtil(4,i) = dwtil(4,i)
     &        + (c2-3.d0)*(c2-2.d0)*(c2-1.d0)*c1*xlamt(i)**(c2-4.d0)
        end do ! i
      end do ! n


c     Compute deviatoric principal Kirchoff stresses and
c     their derivatives wrt xlamt_j

      do i=1,3
        tautil(i) = xlamt(i)*dwtil(1,i) - ( xlamt(1)*dwtil(1,1)
     &                                    + xlamt(2)*dwtil(1,2)
     &                                    + xlamt(3)*dwtil(1,3) )*one3
        dtdl(i,i) = dwtil(2,i)*xlamt(i) + dwtil(1,i)
        do j = 1,3
          dtdl(i,j) = dtdl(i,j)
     &              - ( dwtil(2,j)*xlamt(j) + dwtil(1,j) )*one3
        end do ! j
      end do ! i

c     Compute deviatoric tangent matrix in principal basis

      do i = 1,3
        k = 1+mod(i,3)

        if( abs(xlamt(i)-xlamt(k)).gt.tol )then
          atilp(i+3,i+3) =
     &         ( xlam2(i)*tautil(k) - xlam2(k)*tautil(i) )/
     &         ( xlam2(k) - xlam2(i) )
        else
          c3=1.d0/(xlamt(i)+xlamt(k))
          c4=xlamt(k)-xlamt(i)
          atilp(i+3,i+3) = c3*xlamt(i)*xlamt(k) * ( - dwtil(1,i)
     &                   + dwtil(2,i)*xlamt(i)
     &                   + dwtil(3,i)*xlamt(i)*c4*0.5d0
     &                   + dwtil(4,i)*xlamt(i)*c4**2*one6 )
     &                 + ( dwtil(1,1)*xlamt(1)
     &                   + dwtil(1,2)*xlamt(2)
     &                   + dwtil(1,3)*xlamt(3) )*one3
        endif

        atilp(i,i) = -2.d0*tautil(i)
        do j = 1,3
          atilp(i,j) = atilp(i,j) + xlamt(j)*dtdl(i,j)
     &               - ( xlamt(1)*dtdl(i,1)
     &                 + xlamt(2)*dtdl(i,2)
     &                 + xlamt(3)*dtdl(i,3) )*one3
        end do ! j
      end do ! i

      end
