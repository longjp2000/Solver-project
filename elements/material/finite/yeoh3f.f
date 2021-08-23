c$Id: yeoh3f.f,v 1.1 2006/11/20 20:32:56 rlt Exp $
      subroutine yeoh3f(d,detf,bb, sig,dd, xlamd,ha, estore)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c     1. Use constants from 'pconstant.h'                   14/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose: Compressible Yeoh model
c                (J_2/3 regularization)
c                   _ __      _                    __
c                 W(J,be) = U(J) + 0.5*mu*[J^(-2/3)be:1 - 3)
c                                        _              _
c                                + k1*(I_1 - 2)^2 + k2*(I_1 - 3)^3
c                 _                 __
c                 I_1     = J^(-2/3)be:1 - 3)

c     Input:
c          d(*)    -  Program material parameters (ndd)
c          detf    -  Determinant of deforamtion gradient
c          bb(6)   -  Left Cauchy-Green deformation tensor
c          xlamd   -  Augmented "penalty" value

c     Output:
c          sig(*)  -  Stresses at point.
c                     N.B. 1-d models use only sig(1)
c          dd(6,*) -  Current material tangent moduli
c                     N.B. 1-d models use only dd(1,1) and dd(2,1)
c          ha      -  Augmented function
c          estore  -  Stored energy
c-----[--.----+----.----+----.-----------------------------------------]
      implicit none

      include  'elcount.h'
      include  'pconstant.h'

      integer   i,j, jsw
      real*8    bb(6),be(6), detf,d(*), sig(*),dd(6,*)

      real*8    detfi,j23,trbe3,bei, mub1,mub2,mub3,mub4, i1bar, estore
      real*8    u,up,upp, ha,hp,hpp, k1,k2, press, xlamd

      save

c     Compute deviatoric be

      detfi = 1.d0/detf
      j23   = detfi**two3
      do i = 1,6
        be(i) = bb(i) * j23
      end do
c     call mprint(be,1,6,1,'be')

      i1bar = be(1) + be(2) + be(3)
      trbe3 = i1bar * one3
      i1bar = i1bar - 3.0d0
      be(1) = be(1) - trbe3
      be(2) = be(2) - trbe3
      be(3) = be(3) - trbe3

c     Compute deviatoric Kirchhoff stress tensor.

      mub1 = d(22) ! Yeoh mu
      k1   = d(23)  ! Yeoh k1
      k2   = d(24)  ! Yeoh k2
      mub4 = mub1*4.0d0*(k1 + 3.0d0*k2*i1bar)
      mub1 = mub1*(1.d0 + i1bar*(2.0d0*k1 + 3.0d0*k2*i1bar))
      do i = 1,6
        sig(i) = mub1 * be(i)
      end do

c     Compute tangent tensor

      mub3 = two3 * mub1
      do i = 1,6
        bei = mub4 * be(i)
        do j = 1,6
          dd(i,j) = bei*be(j)
        end do ! j
      end do
      do i = 1,6
        bei = mub3 * be(i)
        do j = 1,3
          dd(i,j) =  dd(i,j) - bei
          dd(j,i) =  dd(j,i) - bei
        end do
      end do

      mub1 = mub1 * trbe3
      mub2 = mub1 + mub1
      mub3 = mub2 * one3

      do i = 1,3
        dd(i  ,i  ) = dd(i  ,i  ) + mub2
        dd(i+3,i+3) = dd(i+3,i+3) + mub1
        do j = 1,3
          dd(i ,j ) = dd(i ,j )   - mub3
        end do
      end do

c     Compute spatial deviatoric stresses and material moduli

      do i = 1,6
        sig(i) = sig(i) * detfi
        do j = 1,6
          dd(i,j) = dd(i,j) * detfi
        end do
      end do

c     Compute pressure and volumetric moduli

      jsw = nint(d(170))
      jsw = max(1,min(3,jsw))
      call fengy3(d,detf, u,up,upp,ha,hp,hpp,jsw)

c     Pressure and tangent (not mixed pressure)

      press =  up  + xlamd * hp
      upp   = (upp + xlamd * hpp) * detf

c     Add volumetric correction to dd

      dd(1,1) = dd(1,1) - press + upp
      dd(1,2) = dd(1,2) + press + upp
      dd(1,3) = dd(1,3) + press + upp

      dd(2,1) = dd(2,1) + press + upp
      dd(2,2) = dd(2,2) - press + upp
      dd(2,3) = dd(2,3) + press + upp

      dd(3,1) = dd(3,1) + press + upp
      dd(3,2) = dd(3,2) + press + upp
      dd(3,3) = dd(3,3) - press + upp

      dd(4,4) = dd(4,4) - press
      dd(5,5) = dd(5,5) - press
      dd(6,6) = dd(6,6) - press

c     Add pressure terms

      do i = 1,3
        sig(i) = sig(i) + press
      end do

c     Energy

      estore = 0.0d0  ! Need to add
      end
