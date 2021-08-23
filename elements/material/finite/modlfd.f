c$Id: modlfd.f,v 1.1 2006/11/20 20:32:56 rlt Exp $
      subroutine modlfd(d,f,finv,df,detf,ta,hn,hn1,nh,istrt, dd,sig,bb,
     &                  xlamd,ha, bbar, isw)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose: Driver for finite deformation constitutive models

c     Inputs:
c       d(*)         - Material parameters array
c       f(3,3)       - Deformation gradient
c       finv(3,3)    - Inverse deformation gradient
c       df(3,3)      - Incremental deformation gradient
c       detf         - Determinant of deformation gradient
c       ta           - Temperature at point
c       hn(*)        - History parameters at t_n
c       nh           - Number history parameters/stress point
c       istrt        - Start state: 0 = elastic; 1 = last solution
c       xlamd        - Augmentation "penalty" value
c       bbar         - Flag (true for B-bar, false for others)
c       isw          - Solution option from elements

c     Outputs
c       hn1(*)       - History parameters at t_n+1
c       dd(*,*,5)    - Expanded material moduli for mixed computation
c                      N.B. Computed from spatial tangent, aa(6,6,5);
c                      Otherwise copy of aa.
c       sig(10)      - Cauchy stress values at t_n+1
c       bb(6)        - Left Cauchy-Green deformation tensor
c       ha           - Augmentation function
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'cdat1.h'
      include  'elcount.h'
      include  'elengy.h'
      include  'sdata.h'

      logical   bbar, plasfl,viscfl, state, conv
      integer   nh,istrt,isw,jsw, ii,i,j, imat, ntm, uprm,umat, stype
      integer   ict
      real*8    ta,xlamd,ha, detf22,f33
      real*8    d(*),  bb(6),  detf(*),f(3,*),finv(3,3),df(3,3)
      real*8    hn(*), hn1(*), dd(*), sig(*), aa(6,6,5)

      save

c     Set number of active stress components

      stype = nint(d(16))

      if(stype.eq.1) then
        detf22  = f(1,1)*f(2,2) - f(1,2)*f(2,1)
        f33     = f(3,3)
        ict     = 0
      endif
      if(ndm.eq.3 .or. stype.eq.8) then
        ntm = 6
      elseif(ndm.eq.2) then
        ntm = 4
      elseif(ndm.eq.1) then
        ntm = 1
      endif

c     Set material model type and user material pointers

      uprm  = ndd-nud
      umat  = int(d(uprm)) - 100

c     Compute Left Cauchy-Green deformation tensor

      bb(1) = f(1,1)*f(1,1) + f(1,2)*f(1,2) + f(1,3)*f(1,3)
      bb(2) = f(2,1)*f(2,1) + f(2,2)*f(2,2) + f(2,3)*f(2,3)
      bb(3) = f(3,1)*f(3,1) + f(3,2)*f(3,2) + f(3,3)*f(3,3)
      bb(4) = f(1,1)*f(2,1) + f(1,2)*f(2,2) + f(1,3)*f(2,3)
      bb(5) = f(2,1)*f(3,1) + f(2,2)*f(3,2) + f(2,3)*f(3,3)
      bb(6) = f(1,1)*f(3,1) + f(1,2)*f(3,2) + f(1,3)*f(3,3)

c     Plane stress return point

100   continue

c     Zero stress and moduli

      do i = 1,6
        sig(i) = 0.0d0
        do j = 1,6
          aa(j,i,1) = 0.0d0
          aa(j,i,2) = 0.0d0
          aa(j,i,3) = 0.0d0
          aa(j,i,4) = 0.0d0
          aa(j,i,5) = 0.0d0
        end do ! j
      end do ! i

c     Set constant initial stress state

      if(nint(d(160)).eq.1) then
        do i = 1,6
          sig(i) = d(160+i)
        end do ! i
      endif

c     Program material models

      if(umat.lt.0) then

c       Set model type

        plasfl = nint(d(40)).eq.1 .or. nint(d(40)).eq.3
        viscfl = nint(d(40)).eq.2

c       Move history variables (in hn-array) to hn1-array

        do i = 1,nh
          hn1(i) = hn(i)
        end do ! i

        imat = d(20)

c       COMPUTE STRESS AND TANGENTS

c       Standard neo-Hookean elastic model

        if    (imat.eq.1) then
          call stnh3f(d,detf,bb, sig,aa,xlamd,ha,estore)

c       Modified neo-Hookean model

        elseif(imat.eq.2) then
          call neoh3f(d,f,finv,detf,bb, hn1,ntm, sig,aa,xlamd,ha,estore)

c       Ogden model

        elseif(imat.eq.3) then
          jsw = nint(d(170))
          call nalp3f(d,f,finv,detf,bb, hn1,ntm, sig,aa,estore, 1,jsw)

c       Finite stretch plasticity model

        elseif(imat.eq.4) then
          if(plasfl) then
            call plasfd(d,detf(1),f,f(1,4), hn1(1),hn1(2),hn1(8),hn(2),
     &                  ntm,istrt, aa,sig,isw,state)
          else
            jsw = nint(d(170))
            call nalp3f(d,f,finv,detf,bb, hn1,ntm, sig,aa,estore, 2,jsw)
          endif

c       Saint-Venant-Kirchhoff model (energy conserving capability)

        elseif(imat.eq.5 .or. imat.eq.6) then
          call stvk(d,detf,f,df,sig,aa, estore)

c       Fung Pseudo-exponential model

        elseif(imat.eq.7) then
          call pfung(d,f,detf,sig,aa, estore)

c       Mooney-Rivlin stress and tangents

        elseif(imat.eq.9) then
          call mnrv3f(d,detf,bb, sig,aa,xlamd,ha,estore)

c       Modified Mooney-Rivlin stress and tangents

        elseif(imat.eq.10) then
          call modmnrv(d,detf,bb, sig,aa,xlamd,ha,estore)

c       Arruda-Boyce model

        elseif(imat.eq.11) then

          call arruda(d,detf,bb, sig,aa,xlamd,ha,estore)

c       Yeoh model

        elseif(imat.eq.12) then

          call yeoh3f(d,detf,bb, sig,aa,xlamd,ha,estore)

        endif

c       Plane stress modification

        if(stype.eq.1) then
          ict       = ict + 1
          call fpstrs(sig,aa,f33,ict, conv,.true.)
          f(3,3)    = f33
          bb(3)     = f33*f33
          finv(3,3) = 1.d0/f33
          detf(1)   = detf22*f33
          if(.not.conv .and. ict.lt.6) go to 100
        endif

c       Plastic update

        if(plasfl) then
          if(state) then
            nomats(2,2) = nomats(2,2) + 1
          else
            nomats(1,2) = nomats(1,2) + 1
          endif

c       Viscoelastic update

        elseif(viscfl) then

          nomats(1,4) = nomats(1,4) + 1

c       Elastic update

        else
          nomats(1,1) = nomats(1,1) + 1
        endif

c     U s e r    M o d e l    I n t e r f a c e

      elseif(umat.le.100) then

        ii = 1 ! Permits later passing of point number
        call umodel(umat,f,detf,ta,d(1),d(uprm+1),hn(1),hn1(1),nh,
     &              ii,istrt, sig,aa, isw)

c     U s e r    F i n i t e    M o d e l    I n t e r f a c e

      else

        umat = umat - 100
        call umodelf(umat,f,finv,df,detf,bb,ta,d(1),d(uprm+1),
     &               hn(1),hn1(1),nh,ntm,istrt, sig,aa, xlamd,ha, isw)

      endif

c     Project aa to D-matrix for B-bar

      if(bbar) then

        call dmatmx ( aa, dd )

      else

        call pmove  ( aa, dd, 36 )

      end if

      end

      subroutine fpstrs(sig,aa,f33,ict, conv,finite)

      implicit   none

      logical    conv, finite
      integer    ict
      real*8     sig(*), aa(6,6), df33,f33, dsig, a3inv, tol

      save

      data       tol /1.d-10/

c     Compute stress

      dsig  = sig(3)

c     Return if either stress or modulus is zero

      if(aa(3,3).eq.0.0d0) then
        conv = .true.
        return
      endif

c     Perform iteration

      a3inv =  1.d0/aa(3,3)
      if(finite) then
        df33  = -dsig*a3inv*f33
      else
        df33  = -dsig*a3inv
      endif
      f33   =  f33 + df33

      if(abs(df33).le.tol*abs(f33)) then
        conv = .true.
      else
        conv = .false.
      endif

c     Reduce tangent array when converged

      if(conv .or. ict.eq.6) then
        aa(1,1) = aa(1,1) - aa(1,3)*a3inv*aa(3,1)
        aa(2,1) = aa(2,1) - aa(2,3)*a3inv*aa(3,1)
        aa(4,1) = aa(4,1) - aa(4,3)*a3inv*aa(3,1)
        aa(1,2) = aa(1,2) - aa(1,3)*a3inv*aa(3,2)
        aa(2,2) = aa(2,2) - aa(2,3)*a3inv*aa(3,2)
        aa(4,2) = aa(4,2) - aa(4,3)*a3inv*aa(3,2)
        aa(1,4) = aa(1,4) - aa(1,3)*a3inv*aa(3,4)
        aa(2,4) = aa(2,4) - aa(2,3)*a3inv*aa(3,4)
        aa(4,4) = aa(4,4) - aa(4,3)*a3inv*aa(3,4)

        aa(1,3) = 0.0d0
        aa(2,3) = 0.0d0
        aa(4,3) = 0.0d0

        aa(3,1) = 0.0d0
        aa(3,2) = 0.0d0
        aa(3,4) = 0.0d0

        aa(3,3) = 0.0d0
      endif

      end
