c$Id: modmnrv.f,v 1.1 2006/11/20 20:32:56 rlt Exp $
      subroutine modmnrv(d,detf,b1, sig,aa,xlamd,ha,engy)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c     1. Use constants from 'pconstant.h'                   14/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Compressible Mooney-Rivlin model
c                (J_2/3 regularization)
c                _ __      _                  _             _
c              W(J,bb) = U(J) + 0.5*mu*(1-c)*(I_1 - 3) + c*(I_2 - 3)
c                  __
c                  bb  = J^(-2/3) * b1
c                  _     __                    __
c                  I_1 = bb:1           (trace bb)
c                  _          _   _     __ __
c                  I_2 = 0.5*(I_1*I_1 - bb:bb)

c      Inputs:

c      Outputs:
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'pconstant.h'
      include  'elcount.h'

      integer   i,j,k,l, ii,jj, jsw, in(3,3)
      real*8    detf, detfi, j23, trbb3,trbb2, bdi, mub1,mub2,mub3
      real*8    nub1,nub2,nub3,nub4, i1b,i2b
      real*8    u, up, upp, ha, hp, hpp, press, xlamd, engy
      real*8    d(*), b1(6), sig(6), aa(6,6)
      real*8    bb(6),bd(6),b2(6)

      save

      data      in / 1, 4, 6, 4, 2, 5, 6, 5, 3 /

c     Compute deviatoric be

      detfi = 1.d0/detf
      j23   = detfi**two3
      do i = 1,6
        bb(i) = b1(i) * j23
      end do ! i

c     Compute bb:bb

      b2(1) = bb(1)*bb(1) + bb(4)*bb(4) + bb(6)*bb(6)
      b2(2) = bb(4)*bb(4) + bb(2)*bb(2) + bb(5)*bb(5)
      b2(3) = bb(6)*bb(6) + bb(5)*bb(5) + bb(3)*bb(3)
      b2(4) = bb(1)*bb(4) + bb(4)*bb(2) + bb(6)*bb(5)
      b2(5) = bb(4)*bb(6) + bb(2)*bb(5) + bb(5)*bb(3)
      b2(6) = bb(6)*bb(1) + bb(5)*bb(4) + bb(3)*bb(6)
      trbb2 =  b2(1)**2 + b2(2)**2 + b2(3)**2
     &      + (b2(4)**2 + b2(5)**2 + b2(6)**2)*2.d0

      i1b   = bb(1) + bb(2) + bb(3)
      i2b   =  0.5d0*(i1b*i1b - trbb2)

      trbb3 = i1b   * one3
      do i = 1,3
        bd(i  ) = bb(i  ) - trbb3
        bd(i+3) = bb(i+3)
      end do ! i

c     Compute deviatoric Kirchhoff stress tensor.

      nub1 = d(23)*d(22)
      mub1 = d(22) - nub1
      do i = 1,6
        sig(i) = mub1 * bd(i) + nub1*(i1b*bd(i) - b2(i))
      end do ! i
      do i = 1,3
        sig(i) = sig(i) - two3*nub1*i2b
      end do ! i

c     Compute tangent tensor
c                                  __             __     _
c     Rank one update: -2/3 mu * ( bd x g +  g x  bd ) / J

      mub3 = two3 * mub1
      do i = 1,6
        bdi = mub3 * bd(i)
        do j = 1,3
          aa(i,j) =  aa(i,j) - bdi
          aa(j,i) =  aa(j,i) - bdi
        end do ! i
      end do ! i
c                       __                     _
c     Deviatoric term 2 mu [ I - 1/3 g x g ] / J

      mub1 = mub1 * trbb3
      mub2 = mub1 + mub1
      mub3 = mub2 * one3

      do i = 1,3
        aa(i  ,i  ) = aa(i  ,i  ) + mub2
        aa(i+3,i+3) = aa(i+3,i+3) + mub1
        do j = 1,3
          aa(i ,j ) = aa(i ,j )   - mub3
        end do ! i
      end do ! i

c     Tangent terms from second invariant

      nub4 = nub1 + nub1
      nub3 = two3 * nub4
      do i = 1,6
        bdi = nub3 * (i1b*bb(i) - b2(i))
        do j = 1,3
          aa(i,j) =  aa(i,j) - bdi
          aa(j,i) =  aa(j,i) - bdi
        end do ! i
      end do ! i

      nub1 = two3 * i2b * nub4
      nub2 = nub1 + nub1
      nub3 = two3 * nub1
      do i = 1,3
        aa(i  ,i  ) = aa(i  ,i  ) + nub2
        aa(i+3,i+3) = aa(i+3,i+3) + nub1
        do j = 1,3
          aa(i ,j ) = aa(i ,j )   + nub3
        end do ! i
      end do ! i

      do l = 1,3
        do k = l,3
          jj = in(k,l)
          do j = 1,3
            do i = j,3
              ii = in(i,j)
              aa(ii,jj) = aa(ii,jj) + nub4*(bb(ii)*bb(jj)
     &                  - 0.5d0*(bb(in(i,k))*bb(in(j,l))
     &                         + bb(in(i,l))*bb(in(j,k))))
            end do ! i
          end do ! j
        end do ! k
      end do ! l

c     Compute spatial deviatoric stresses and material moduli

      do i = 1,6
        sig(i) = sig(i) * detfi
        do j = 1,6
          aa(i,j) = aa(i,j) * detfi
        end do ! i
      end do ! i

c     Compute pressure and volumetric moduli

      jsw = nint(d(170))
      call fengy3(d,detf, u,up,upp,ha,hp,hpp,jsw)

c     Pressure and tangent (not mixed pressure)

      press =  up  + xlamd * hp
      upp   = (upp + xlamd * hpp) * detf

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

c     Add pressure terms

      do i = 1,3
        sig(i) = sig(i) + press
      end do ! i

c     Compute stored energy density

      engy  =  u + d(22)*(1.d0 - d(23))*(i1b - 3.d0)*0.5d0
     &                   + d(22)*d(23) *(i2b - 3.d0)*0.5d0

      end
