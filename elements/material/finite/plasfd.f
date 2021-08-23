c$Id: plasfd.f,v 1.1 2006/11/20 20:32:56 rlt Exp $
      subroutine plasfd(d, detf, fn1, fn, epp, be, epl, bn,
     &                  ntm, istrt, dd, sig, isw, state)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c     1. Use constants from 'pconstant.h'                   14/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose: Finite Deformation Isotropic I1-J2-J3 Plasticity Models
c              in Principal Logarithmic Stretches
c              Incremental  Lagrangian  Formulation

c     Inputs:

c       d(*)      -  Material parameters
c       detf      -  Jacobian determinant (F)        at t_n+1
c       fn1(3,3)  -  Deformation gradient            at t_n+1
c       fn(3,3)   -  Deformation gradient            at t_n
c                 History Variables
c       epp       -  Cumulative plastic strain       at t_n
c       be(*)     -  Left Cauchy-Green tensor        at t_n
c       epl(*)    -  Plastic Strain for Hardening    at t_n

c       istrt     -  Start state: 0 = elastic; 1 = last solution

c     Outputs:

c                 History Variables
c       epp       -  Cumulative plastic strain       at t_n+1
c       be(*)     -  Left  Cauchy-Green tensor       at t_n+1
c       epl(*)    -  Plastic Strain for Hardening    at t_n+1
c                 History Initialization Variables
c       bn(*)     -  Left Cauchy-Green tensor        at t_n   (isw=14 only)

c       sig(*)    -  Cauchy stress tensor
c       dd(6,6)   -  Cauchy (spatial) elastic moduli
c-----[--.----+----.----+----.-----------------------------------------]
c     MATERIAL CONSTANTS:
c       d(21)    -  Elastic Bulk  modulus
c       d(22)    -  Elastic Shear modulus

c       d(44)    -  Isotropic Hardening Modulus - H_iso
c       d(45)    -  Kinematic Hardening Modulus - H_kin

c       d(46)    -  1 = VON MISES      [ sqrt(2*J2)                           ]
c       d(41)    -  Yield stress (ep = 0)(Y0)
c       d(42)    -  Yield stress (ep =00)(Yi)
c       d(43)    -  Delta value  (ep = 0)(delta)

c       d(46)    -  2 = DRUCKER PRAGER [ sqrt(2*J2) + (1/3) * aa * I1         ]
c       d(41)    -  Yield stress         (Y0)
c       d(42)    -  aa parameter         (aa)

c       d(46)    -  3 = PRAGER LODE    [ sqrt(2*J2) + sqrt(27/2) * bb * J3/J2 ]
c       d(41)    -  Yield stress         (Y0)
c       d(42)    -  Lode angle parameter (bb)

c     VARIABLES:
c       vol       - volumetric strain   -   vol = ( eps : 1 ) / 3
c       be(*)     - Left Cauchy-Green tensor
c       nn_t(3,3) - principal directions (by columns) tensor form
c       ll2(3)    - squares of principal stretches =  lambda^2

c       tau(3)    - principal values   total    Kirchhoff stress
c       tt(3)     - principal values deviatoric Kirchhoff stress
c       pp        - pressure         volumetric Kirchhoff stress

c       dtde(3,3) - Kirchhoff stress derivative
c                   dtde(a,b)   = [ d tau_a / d (lambda_b^2) ] * lambda_b^2
c                               = [ d tau_a / d ( eps_b ) ] / 2.0d0

c       dd_l(6,6)   - spatial elastic moduli in principal basis

c     FUNCTIONS:
c       yield     - yield function value
c       yfunct    - yield function
c                    flg = .false.  only yield function
c                    flg = .true.   yield function and derivatives
c-----[--.----+----.----+----.-----------------------------------------]
      implicit none

      include 'counts.h'
      include 'iofile.h'
      include 'pconstant.h'
      include 'rdata.h'

      logical  conv, state
      integer  ntm,istrt, i, j, a, b, rot, it, itmax, isw, p1(6), p2(6)
      real*8   d(*), f(3,3), be(*),bn(*),epl(*),dd(6,6), sig(*)
      real*8   fn1(3,3),fn(3,3),fi(3,3),detfr,detf, detr
      real*8   pp, kk , gg, TwoG, K3inv, G2inv, Ypr, dyld
      real*8   ll2_tr(3), ll2(3), tau(3), tt(3)
      real*8   eps_tr(3), th_tr, vol_tr, ee_tr(3), alp_n(3)
      real*8   eps_e(3),         vol_e,  ee_e(3) , alp(3), ta(3)
      real*8   Y0, I1, J2, J3, f1, f2, f3, f11, f22, f33, f12, f13, f23
      real*8   gam, Hk, Hkr, aatr, epse
      real*8   err, dsol(7),res(7), tres(7,7),fss(3,3)
      real*8   dtde(3,3), ff(6,6), dd_l(6,6), epp, eppn
      real*8   nn_t(3,3), nn(3), aa(3), bb(3)
      real*8   tolc, One2, d_el, xx,f112,f123, yield,yfunct

      data     p1    /1,2,3,1,2,3/
      data     p2    /1,2,3,2,3,1/
      data     One2  / 0.5d0 /
      data     dd_l  / 36*0.0d0 /

c     INITIALIZE HISTORY VARIABLES: ISW = 14

      if (isw.eq.14) then        ! Initialize non-zero quantities
        do a = 1,3
          be(a) = 1.0d0
          bn(a) = 1.0d0
        end do ! a
        return
      end if

c     READING MATERIAL PARAMETERS  &  CONSTANT SETUP

      kk    = d(21)
      gg    = d(22)
      Y0    = d(41)

      Hk    = d(45)
      if(Hk.gt.0.0d0) then
        Hkr = 1.d0/Hk
      else
        Hkr = 0.0d0
      endif

      TwoG  = 2.0d0 * gg

      itmax = 50                      ! Iteration maximum #
      tolc  = 1.d-09                  ! Tolerance for convergence

c     Check state for iterations

      if(niter.eq.0) then         ! First iteration in step
        if(istrt.eq.0) then       ! Elastic state requested
          state = .false.
          dyld  =  0.0d0
        else                      ! Last state requested
          state = .true.
          dyld  =  1.0d-08*d(41)
        endif
      else                        ! Not first iteration in step
        state = .true.
        dyld  =  0.0d0
      endif

c     KINEMATIC COMPUTATIONS  (Incremental Lagrangian formulation)

c     Inverse of fn(3,3) * J_n

      fi(1,1) = fn(2,2)*fn(3,3) - fn(2,3)*fn(3,2)
      fi(2,1) = fn(2,3)*fn(3,1) - fn(2,1)*fn(3,3)
      fi(3,1) = fn(2,1)*fn(3,2) - fn(2,2)*fn(3,1)

      fi(1,2) = fn(3,2)*fn(1,3) - fn(3,3)*fn(1,2)
      fi(2,2) = fn(3,3)*fn(1,1) - fn(3,1)*fn(1,3)
      fi(3,2) = fn(3,1)*fn(1,2) - fn(3,2)*fn(1,1)

      fi(1,3) = fn(1,2)*fn(2,3) - fn(1,3)*fn(2,2)
      fi(2,3) = fn(1,3)*fn(2,1) - fn(1,1)*fn(2,3)
      fi(3,3) = fn(1,1)*fn(2,2) - fn(1,2)*fn(2,1)

      detfr   = 1.d0/(fn(1,1)*fi(1,1)+fn(1,2)*fi(2,1)+fn(1,3)*fi(3,1))

      do j = 1,3
        do i = 1,3
          f(i,j) = (fn1(i,1)*fi(1,j)
     &            + fn1(i,2)*fi(2,j)
     &            + fn1(i,3)*fi(3,j))*detfr
        end do ! i
      end do ! j

c     Compute Elastic left Cauchy-Green tensor:    be = F * Cp * F^T

      call pushr2(f,be,be,1.0d0)      ! Update configuration

c     Compute principal stretches and directions

      nn_t(1,1) = be(1)
      nn_t(2,2) = be(2)
      nn_t(3,3) = be(3)

      nn_t(1,2) = be(4)
      nn_t(2,1) = be(4)

      nn_t(2,3) = be(5)
      nn_t(3,2) = be(5)

      nn_t(1,3) = be(6)
      nn_t(3,1) = be(6)

      call eig3(nn_t,ll2_tr,rot)

c     COMPUTE TRIAL KIRCHHOFF STRESS  ( pressure and deviator )

      do a = 1, 3
        ee_tr(a) = log( sqrt(ll2_tr(a)) )  ! log ( lambda(a)^TR )
      end do ! a
      th_tr = ee_tr(1) + ee_tr(2) + ee_tr(3)
      do a = 1,3
        ee_tr(a) = ee_tr(a) - one3*th_tr
      end do ! a

      pp   = kk * th_tr                    ! Pressure: K*th_tr
      eppn = epp

      do a = 1, 3
        tt(a)    = TwoG   * ee_tr(a)
        tau(a)   = tt(a)  + pp
        alp(a)   = Hk * epl(a)
        alp_n(a) = alp(a)
      end do ! a

c     Deviatoric: ta =  tt - alp_dev

      aatr = (alp(1) + alp(2) + alp(3))*one3
      do a = 1,3
        ta(a)   = tt(a) - alp(a) + aatr
      end do ! a

c     CHECK ELASTIC / PLASTIC STEP

c     Compute stress invariant and yield function

      I1 =  3.0d0 * (pp - aatr)
      J2 = ( ta(1)*ta(1) + ta(2)*ta(2) + ta(3)*ta(3) ) * One2
      J3 = ( ta(1)*ta(1)*ta(1) + ta(2)*ta(2)*ta(2) +
     &                           ta(3)*ta(3)*ta(3)   ) * one3

      yield = yfunct(d,epp,I1,J2,J3, f1,f2,f3,
     &               f11,f22,f33,f12,f13,f23,Ypr,.false.)

c     Check yield

      if ( (yield+dyld .gt. 0.0d0) .and. state ) then

c     PLASTIC STEP  -->  RETURN MAP

        conv   = .false.
        it     =  1
        gam    =  0.0d0

        K3inv  = one3 /  kk
        G2inv  = One2 /  gg

        vol_tr = th_tr * one3
        vol_e  = K3inv * pp

        do a = 1, 3
          eps_tr(a) = ee_tr(a) + vol_tr
          ee_e(a)   = G2inv * tt(a)
          eps_e(a)  = ee_e(a)  + vol_e
        end do ! a
        epse = 1.d0/(abs(eps_tr(1)) + abs(eps_tr(2)) + abs(eps_tr(3)))

        d_el = ( K3inv - G2inv ) * one3

        do while ((.not.conv).and.(it.le.itmax))

          yield = yfunct(d,epp,I1,J2,J3, f1,f2,f3,
     &                   f11,f22,f33,f12,f13,f23,Ypr,.true.)

          xx  = f1 - f3 * two3 * J2

          do a = 1, 3
            nn(a)    = xx + ( f2 + f3 * ta(a) ) * ta(a)
            res(a)   = eps_e(a) - eps_tr(a) + gam * nn(a)
            res(a+3) = (alp_n(a) - alp(a))*Hkr + gam * nn(a)
c           res(a+3) = 0.0d0
          end do ! a

          res(7) = yield

          err  = (abs(res(1)) + abs(res(2)) + abs(res(3)))*epse
     &         +  abs(res(7))/Y0
          conv = err .lt. tolc

c         Construct local tangent matrix

          do a = 1, 3
            aa(a) = f23   * ta(a) + f13
            bb(a) = ta(a) * ta(a) - two3 * J2
          end do ! a

          f112 = f11 - one3 * f2
          f123 = f12 - two3 * f3

          do a = 1, 3
            do b = a, 3
              fss(a,b) = ( f112 + f22 * ta(a) * ta(b)
     &                 +   f33 * bb(a) * bb(b)
     &                 +  f123 * ( ta(b) + ta(a) )
     &                 + aa(a) * bb(b) + bb(a) * aa(b) ) * gam
              fss(b,a) = fss(a,b)
            end do ! b
            fss(a,a) = fss(a,a) + gam * ( f2 + 2.0d0 * f3 * ta(a) )
          end do ! a

          do a = 1, 3
            do b = 1, 3
              tres(a,b) = fss(a,b) + d_el
            end do ! b
            tres(a,a) = tres(a,a) + G2inv
          end do ! a

c         Kinematic and Isotropic Hardening

          if(Hk.gt.0.0d0) then

            do a = 1,3
              do b = 1,3
                tres(a,b+3)   = -fss(a,b)
                tres(a+3,b)   = -fss(a,b)
                tres(a+3,b+3) =  fss(a,b)
              end do ! b
              tres(a+3,a+3) = tres(a+3,a+3) + 1.d0/Hk
              tres(a  ,7) =  nn(a)
              tres(a+3,7) = -nn(a)
              tres(7,a  ) =  nn(a)
              tres(7,a+3) = -nn(a)
            end do ! a

            tres(7,7) = -Ypr

            call invert(tres,7,7)

            do i = 1, 7
              dsol(i) = - tres(i,1)*res(1)
     &                  - tres(i,2)*res(2)
     &                  - tres(i,3)*res(3)
     &                  - tres(i,4)*res(4)
     &                  - tres(i,5)*res(5)
     &                  - tres(i,6)*res(6)
     &                  - tres(i,7)*res(7)
            end do ! i

c         Isotropic hardening only

          else
            do a = 1,3
              tres(a,4) = nn(a)
              tres(4,a) = nn(a)
            end do ! a

            tres(4,4) = -Ypr

            call invert(tres,4,7)

            do i = 1, 4
              dsol(i) = - tres(i,1)*res(1)
     &                  - tres(i,2)*res(2)
     &                  - tres(i,3)*res(3)
     &                  - tres(i,4)*res(7)
            end do ! i
            dsol(7) = dsol(4)
            dsol(4) = 0.0d0
            dsol(5) = 0.0d0
            dsol(6) = 0.0d0
          endif

c         Update Kirchhoff stress and plastic flow

          tau(1) = tau(1) + dsol(1)
          tau(2) = tau(2) + dsol(2)
          tau(3) = tau(3) + dsol(3)
          gam    = gam    + dsol(7)

c         Update accumulated plastic strain

          epp    = eppn + gam

c         Update Back Stress

          alp(1) = alp(1) + dsol(4)
          alp(2) = alp(2) + dsol(5)
          alp(3) = alp(3) + dsol(6)

c         Update vol.-dev. Kirchhoff stress and stress invariants

          pp     = ( tau(1) + tau(2) + tau(3) ) * one3
          tt(1)  = tau(1) - pp
          tt(2)  = tau(2) - pp
          tt(3)  = tau(3) - pp

          aatr = (alp(1) + alp(2) + alp(3))*one3
          do a = 1,3
            ta(a)   = tt(a) - alp(a) + aatr
          end do ! a

          I1 =  3.0d0 * (pp - aatr)
          J2 = ( ta(1)*ta(1)       + ta(2)*ta(2)       +
     &                               ta(3)*ta(3)         ) * One2
          J3 = ( ta(1)*ta(1)*ta(1) + ta(2)*ta(2)*ta(2) +
     &                               ta(3)*ta(3)*ta(3)   ) * one3

c         Update vol.-dev. logarithmic strain

          vol_e  = K3inv * pp

          do a = 1, 3
            ee_e(a)  = G2inv * tt(a)
            eps_e(a) = ee_e(a) + vol_e
          end do ! a

          it = it + 1

        end do ! while

c       Warning: check convergence

        if(.not.conv .and. niter.gt.0) then
          write(  *,*) ' *WARNING* No convergence in PLASFD',err,tolc
          write(iow,*) ' *WARNING* No convergence in PLASFD',err,tolc
          call plstop()
        endif

c       Update elastic left Cauchy-Green tensor and plastic acc. strain

        do a = 1, 3
          ll2 (a) = exp( 2.0d0 * eps_e(a) )
        end do ! a
        do i = 1, 6
          be(i) = ll2(1) * nn_t(p1(i),1)*nn_t(p2(i),1)
     &          + ll2(2) * nn_t(p1(i),2)*nn_t(p2(i),2)
     &          + ll2(3) * nn_t(p1(i),3)*nn_t(p2(i),3)
        end do ! i

c       Update plastic strains

        epl(1) = epl(1) + gam * nn(1)
        epl(2) = epl(2) + gam * nn(2)
        epl(3) = epl(3) + gam * nn(3)

c       Compute elasto-plastic tangent

        do b = 1, 3
          do a = 1, 3
            dtde(a,b) = 0.5d0 * tres(a,b)
          end do ! a
        end do ! b

      else

c     ELASTIC STEP  ( only tangent computation )

        state  = .false.                   ! Indicate elastic on return

        dtde(1,1) = One2 * kk  + two3 * gg
        dtde(1,2) = One2 * kk  - one3 * gg
        dtde(1,3) = dtde(1,2)
        dtde(2,1) = dtde(1,2)
        dtde(2,2) = dtde(1,1)
        dtde(2,3) = dtde(1,2)
        dtde(3,1) = dtde(1,2)
        dtde(3,2) = dtde(1,2)
        dtde(3,3) = dtde(1,1)

      end if

c     COMPUTE CAUCHY STRESS

      detr = 1.d0 / detf
      do i = 1, ntm
        sig(i) = (tau(1) * nn_t(p1(i),1)*nn_t(p2(i),1)
     &         +  tau(2) * nn_t(p1(i),2)*nn_t(p2(i),2)
     &         +  tau(3) * nn_t(p1(i),3)*nn_t(p2(i),3)) * detr
      end do ! i

c     TANGENT TRANSFORMATION

c     Material tangent (computation in the principal basis)

      do a = 1, 3

C       Upper 3x3 block of dd_l()

        do b = 1, 3
          dd_l(b,a) = 2.0d0 * dtde(b,a)
        end do ! b
        dd_l(a,a) = dd_l(a,a) - 2.0d0 * tau(a)

C       Lower 3x3 block of dd_l() [ diagonal block ]

        b = mod(a,3) + 1
        if (abs(ll2_tr(a)-ll2_tr(b)).gt.tol) then
          dd_l(a+3,a+3) = ( ll2_tr(a)*tau(b) - ll2_tr(b)*tau(a) ) /
     &                    ( ll2_tr(b) - ll2_tr(a))
        else
          dd_l(a+3,a+3) =  dtde(a,a) - dtde(b,a) - tau(a)
        endif

      end do ! a

c     Transform matrix to standard basis

      call tranr4(nn_t,nn_t,ff)
      call pushr4(ff,ff,dd_l,dd,detf)

      end

      function yfunct(d,epp,I1,J2,J3,f1,f2,f3,
     &                f11,f22,f33,f12,f13,f23,Ypr,flg)

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c     1. Use constants from 'pconstant.h'                   14/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
      implicit none

      include 'iofile.h'
      include 'pconstant.h'

      real*8  yfunct, yield

      real*8  d(*),epp
      real*8  I1, J2, J3
      real*8  f1, f2, f3, f11, f22, f33, f12, f13, f23, Ypr
      real*8  J2_12, J2_1, J2_2
      real*8  YY, aa, bb, ss, c1, c2
      logical flg

      integer iyield

      f1  = 0.0d0
      f2  = 0.0d0
      f3  = 0.0d0
      f11 = 0.0d0
      f22 = 0.0d0
      f33 = 0.0d0
      f12 = 0.0d0
      f13 = 0.0d0
      f23 = 0.0d0

      iyield = d(46)
      YY     = d(41) + d(44)*epp
      Ypr    = d(44)

c     VON MISES        YIELD FUNCTION

      if (iyield.eq.1) then

        aa    = (d(41) - d(42))*exp(-d(43)*epp)
        YY    = d(42) + aa + d(44)*epp
        Ypr   = Ypr - d(43)*aa
        ss    = sqrt( 2.0d0 * J2 )
        yield = ss - YY

c       Compute yield function derivatives

        if (flg) then

          if (ss.ne.0.0d0) then
            f2  =   1.0d0 / ss
            f22 = - f2**3
          end if

        end if

c     DRUCKER PRAGER   YIELD FUNCTION

      elseif (iyield.eq.2) then

        aa    = d(42)
        ss    = sqrt( 2.0d0 * J2 )
        yield = ss + one3 * aa * I1 - YY

c       Compute yield function derivatives

        if (flg) then

          f1  = one3 * aa

          if (ss.ne.0.0d0) then
            f2  =   1.0d0 / ss
            f22 = - f2**3
          end if

        end if

c     PRAGER-LODE      YIELD FUNCTION

      elseif (iyield.eq.3) then

        if (J2.gt.0.0d0) then

c         Read in material parameters

          bb = d(42)

          c1 = 1.0d0 / sqrt(2.0d0)
          c2 = sqrt(13.5d0) * bb

c         Compute yield function

c         yield = J2 * ( 1.0d0 + bb * J3 * J2 **(-1.5d0)) - YY
          yield = ( sqrt(2.0d0*J2) + c2 * J3 / J2 ) - YY

c         Compute yield function derivatives

          if (flg) then

            J2_1    =   1.d0/J2
            J2_2    =   J2_1**2
            J2_12   =   sqrt(J2_1)

            f2      =   c1 * J2_12            - c2 * J3 * J2_2
            f3      =                           c2      * J2_1
            f22     = - c1 * J2_12**3 * 0.5d0 + c2 * J3 * J2_1**3 * 2.d0
            f23     =                         - c2      * J2_2
          end if

        else
          yield = -YY
        end if

c     NO YIELD FUNCTION SPECIFIED

      else

        yield = 0.0d0
        write (iow,9000)
        if (ior.lt.0) then
          write (*,9000)
        endif

      end if

      yfunct = yield

c     Format

9000  format(' *** NO YIELD FUNCTION SPECIFIED ***')

      end
