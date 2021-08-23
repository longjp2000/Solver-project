c$Id: modlsd.f,v 1.1 2006/11/20 20:32:56 rlt Exp $
      subroutine modlsd(d,ta,eps,h1,h2,nh,istrt, dd,sig,isw)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose: Small Deformation Constitutive Equation Driver

c     Input parameters
c          d(*)      -  up to ndd-nud-1 material parameters
c          eps(9,3)  -  current strains and fluxes at point
c          h(nh)     -  history terms at point
c          nh        -  number of history terms
c          istrt     -  Start state: 0 = elastic; 1 = last solution
c          im        -  material type
c     Ouput parameters
c          dd(6,6,5) -  current material tangent moduli
c                       Coupled problems:    | dd_1   dd_2 |
c                                            | dd_3   dd_4 |
c                       Rayleigh damping: dd_5
c          sig(6)    -  stresses at point.
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'cdat1.h'
      include  'complx.h'
      include  'elcount.h'
      include  'pmod2d.h'
      include  'sdata.h'

      logical   state, conv
      integer   ii,it,i,j,nh,istrt,isw,ict, ntm, umat,uprm
      real*8    ta,e33
      real*8    d(*),eps(9,*),h1(*),h2(*),dd(6,6,5),sig(*),theta(3)

      save

c     Extract analysis type: 1=plane stress; 2=plane strain; 3=axisy;
c                            8=axisy+torsion

      stype = d(16)
      if(stype.eq.1) then
        e33 = eps(3,1)
        ict = 0
      endif

c     Check for user model

      uprm  = ndd-nud
      umat  = int(d(uprm)) - 100

c     Set number stress/strain components

      if(ndm.eq.3 .or. stype.eq.8) then
        ntm = 6
      elseif(ndm.eq.2) then
        ntm = 4
      else
        ntm = 1
      end if

c     Plane stress return point

100   continue

c     Zero stress and dd arrays

      do i = 1,6
        sig(i) = 0.0d0
        do j = 1,6
          dd(j,i,1) = 0.0d0
          dd(j,i,2) = 0.0d0
          dd(j,i,3) = 0.0d0
          dd(j,i,4) = 0.0d0
          dd(j,i,5) = 0.0d0
        end do ! j
      end do ! i

c     Set constant initial stresses

      if(nint(d(160)).eq.1) then
        do i = 1,6
          sig(i) = d(160+i)
        end do !
      end if

c     Program material models

      if(umat.lt.0) then

c       Set model type

        plasfl = nint(d(40)).eq.1 .or. nint(d(40)).eq.3
        viscfl = nint(d(40)).eq.2

c       Move h1 to h2

        do i = 1,nh
          h2(i) = h1(i)
        end do ! i

c       P l a s t i c i t y

        if(plasfl) then

          if(nint(d(40)).eq.1) then

c           Plane stress plasticity

            if (stype.eq.1 .or.stype.eq.4) then

              call epps2d(d,eps,h2,h2(4),h2(7),istrt, sig,dd,dd(1,1,5))
              state = h2(8).eq.0.0d0
              it    = 2

c           Plane strain, axisymmetric or 3D plasticity

            else

              call mises(d,eps,h2(3),h2,ntm,istrt, sig,dd,dd(1,1,5))
              state = h2(2).eq.0.0d0
              it    = 2
              if(.not.state .and. isw.eq.8) sig(10) = h2(1)

            endif

c         G e n e r a l i z e d    P l a s t i c i t y

          elseif(nint(d(40)).eq.3) then

            call gplas3d(d,eps,h2(4),h2,ntm,istrt, sig,dd,dd(1,1,5))
            state = h2(2).eq.0.0d0
            it    = 3

          endif

c       E l a s t i c i t y

        else

c         Piezoelectric case

          if(d(150).gt.0.0d0) then

            call pzstrs(d,eps(7,1),eps,sig,dd)

c         Thermoelastic case

          else

            call estrsd(d,ta,eps,sig,dd,dd(1,1,5))
            nomats(1,1) = nomats(1,1) + 1

          endif

        end if

c       V i s c o e l a s t i c i t y

        if(viscfl) then

c         Complex modulus form

          if(cplxfl) then
            call cvisco(d,eps,sig,dd,dd(1,1,3))

c         Time moduli integration

          else
            call viscoe(d,ta,eps,h2(1),h2(ntm+1),ntm,sig,dd,dd(1,1,5))
          endif

        end if

c     U s e r    M o d e l    I n t e r f a c e

      else

c       Compute trace to pass to user module

        theta(1) = eps(1,1) + eps(2,1) + eps(3,1)
        theta(2) = eps(1,2) + eps(2,2) + eps(3,2)
        theta(3) = eps(1,3) + eps(2,3) + eps(3,3)

        ii = 1 ! Permits later addition of quadrature point.
        call umodel(umat,eps,theta,ta,d(1),d(uprm+1),h1(1),h2(1),nh,
     &              ii,istrt,sig,dd,isw)

      end if

c     Plane stress modification

      if(stype.eq.1) then
        ict       = ict + 1
        call fpstrs(sig,dd,e33,ict, conv,.false.)
        eps(3,1)  = e33
        if(.not.conv .and. ict.lt.6) go to 100
      endif

c     User model

      if(umat.gt.0) then ! Currently does nothing

c     Plastic update

      elseif(plasfl) then
        if(state) then
          nomats(2,it) = nomats(2,it) + 1
        else
          nomats(1,it) = nomats(1,it) + 1
        endif

c     Viscoelastic update

      elseif(viscfl) then

        nomats(1,4) = nomats(1,4) + 1

c     Piezoelectric case

      elseif(d(150).gt.0.0d0) then
        nomats(1,5) = nomats(1,5) + 1

c     Elastic update

      else
        nomats(1,1) = nomats(1,1) + 1
      endif

      end
