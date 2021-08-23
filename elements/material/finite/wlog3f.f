c$Id: wlog3f.f,v 1.1 2006/11/20 20:32:56 rlt Exp $
      subroutine wlog3f( d,bpr, epsd,taup,aap,w )

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c     1. Use constants from 'pconstant.h'                   14/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose: Compute stress and tangent in principal stretch form.

c              Stored Energy Function in Principlal Stretches
c                   W       = w(lam_1) + w(lam_2) + w(lam_3)
c                   taup_i  = [w(lam_i)]' * lambda_i
c                   gamm_ij = lambda_j * [taup_i],j

c              Logarithmic strain energy function
c                   w    = mu * [log(lambda)]**2

c     Inputs:
c          d(*)     Material parameters
c          bpr(3)   Principal stretch - squared

c     Outputs:
c          epsd(3)  Principal deviatoric strains
c          taup(3)  Principal deviatoric Kirchhoff stresses
c          aap(6,6) Moduli in principal basis
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'pconstant.h'

      integer   i, j
      real*8    jthird, twomu, fourmu, w
      real*8    bpr(3),d(*),lamt(3), epsd(3), taup(3), aap(6,6)

      real*8    tol
      data      tol / 1.0d-08 /

c     Compute Je**(-1/3) and volume preserving stretches

      jthird = (bpr(1)*bpr(2)*bpr(3))**(-one6)
      do i = 1,3
        lamt(i) = jthird*sqrt(bpr(i))
      end do ! i

c     Compute predictor strain energy function derivatives

      twomu  = 2.d0   * d(22)
      w      = 0.d0
      do i = 1,3
        epsd(i) = log(lamt(i))
        taup(i) = twomu * epsd(i)
        w       = w +  d(22)*epsd(i)**2
      end do ! i

c     Load material part of moduli in principal stretches

      do i = 1,6
        do j = 1,6
          aap(j,i) = 0.0d0
        end do ! j
      end do ! i

      twomu  = twomu * one3
      fourmu = twomu + twomu

      do i = 1,3
        aap(i,i) = fourmu - 2.d0*taup(i)
      end do ! i

      aap(1,2) = - twomu
      aap(2,1) = - twomu

      aap(2,3) = - twomu
      aap(3,2) = - twomu

      aap(3,1) = - twomu
      aap(1,3) = - twomu

c     Compute deviatoric tangent matrix in principal basis

      do i = 1,3
        j = mod(i,3) + 1

        if (abs(bpr(i) - bpr(j)) .gt. tol) then
          aap(i+3,i+3) = ( bpr(i)*taup(j) - bpr(j)*taup(i) )
     &                 / ( bpr(j) - bpr(i) )
        else
          aap(i+3,i+3) = aap(i,i) - aap(j,i)
        endif

      end do ! i

      end
