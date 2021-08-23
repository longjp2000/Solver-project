c$Id: inmate.f,v 1.1 2006/11/20 20:32:56 rlt Exp $
      subroutine inmate(d,tdof,ndv,eltype)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c     1. Use 'pi' from 'pconstant.h'                        14/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Input material parameters for material models

c      Inputs:
c         ietype    - Element type indicator
c                     1: Mechanical solid
c                     2: Truss
c                     3: Frame
c                     4: Plate
c                     5: Shell/Membrane
c                     6: Thermal solid
c                     7: Three-dimensional
c                     8: Unspecified

c      Outputs:
c         d(*)      - Material set parameters
c         tdof      - Dof to specify temperature location in
c                     coupled thermo-mechanical problems
c         ndv       - Number of element history variables
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'cdata.h'
      include  'cdat1.h'
      include  'elplot.h'
      include  'hdata.h'
      include  'iofile.h'
      include  'pconstant.h'
      include  'pglob1.h'
      include  'pmod2d.h'
      include  'refnd.h'
      include  'refng.h'
      include  'rigid1.h'
      include  'sdata.h'
      include  'strnum.h'
      include  'umac1.h'

      character text(2)*15,wd(9)*13,qd(2)*9
      logical   pcomp,erflag,errck,pinput,tinput,inputs,incomp,damgfl
      logical   cflag,eflag,fflag,iflag,sflag,tflag,uflag,efl,flg,oned
      logical   uparfl
      integer   eltype,ietype,imat,ndv,nvis,tdof,i,ii,j,jj,nn,umat,uprm
      integer   n1,n3,nlay,ntub,nlob,nseg,nrct,nqud,nhd,ntm,naug,nsiz
      real*8    bulk,e1,e2,e3,nu12,nu23,nu31,g12,g23,g31,rad,sigc,sigt
      real*8    ev(10),dd(6,6),alp(3),d(*)

      save

      data      wd/'Plane Stress' ,'Plane Strain' ,'Axisymmetric',
     &             'Plate & Shell','Truss & Frame','Thermal',
     &             '3-Dimensional','Axisy+Torsion','Unspecified' /

      data      qd/'G a u s s' ,'N o d a l'/

c     PART 1: Set default values

      tdof   = gtdof
      etype  = 1
      dtype  = gdtype
      ietype = abs(eltype)

      if(dtype.gt.0 .or. ietype.ne.1) then
        fflag  = .false.
        sflag  = .true.
      else
        fflag  = .true.
        sflag  = .false.
      endif

      cflag  = .false.
      eflag  = .false.
      erflag = .false.
      hflag  = .false.
      iflag  = .false.
      incomp = .false.
      oned   = .false.
      plasfl = .false.
      tflag  = .false.
      uflag  = .false.
      damgfl = .false.
      viscfl = .false.
      uparfl = .false.

      if(nen.gt.11.and. ndm.eq.2) then
        ii = 4
      elseif(nen.gt.4 .and. ndm.eq.2 .or.
     &   nen.gt.8 .and. ndm.eq.3 ) then
        ii = 3
      else
        ii = 2
      endif
      jj   = ii - 1
      imat = 1
      lref = 0
      naug = 0
      nvis = 0
      nlay = 0
      nqud = 0
      nrct = 0
      nseg = 0
      ntub = 0
      nlob = 2

c     Default rigid body number (denotes flexible)

      nrmatl= 0

c     Default Mass: Consistent

      d( 7) = 1.d0

c     Default angle in degrees

      rad   = pi/180.0d0
      d(31) = d(31)/rad

c     Solid type

      if(ietype.eq.1 .or. ietype.eq.6) then
        stype  = g2type
        d(170) = 1.d0  ! Default volume model

c     Truss/frame type

      elseif(ietype.eq.2 .or.ietype.eq.3) then
        if(ietype.eq.2) oned   = .true.
        stype = 5
        lref  = gref
        sref  = 0
        if(lref.eq.1) then
          do i = 1,ndm
            refx(i) = grefx(i)
            tref(i) = 0.0d0
          end do ! i
        else
          do i = 1,ndm
            refx(i) = gtref(i)
            tref(i) = 0.0d0
          end do ! i
        end if

c     Plate/shell type

      elseif(ietype.eq.4) then
        stype = 4
        ii    = 3
      elseif(ietype.eq.5) then
        stype = 4
        ii    = 2
      elseif(ietype.eq.7) then
        stype = 7
      endif

c     Angular velocity and Rayleigh damping set

      if(ietype.ne.6) then
        d(65) = gomega(1)
        d(77) = gray(1)
        d(78) = gray(2)
      endif

c     Quadrature type

      if(gquadn.gt.0.0d0) then
        d(182) = 1.d0                         ! Nodal Quadrature 
      endif

      d(14)  = 1.d0
      alp(1) = 0.d0

c     PART 2: Poll for legal input records: Stop on blank record

      inputs = .true.

      do while(inputs)

c       Input record

        errck = .true.
        do while(errck)
          errck = tinput(text,2,ev,10)
        end do ! while

c       Plate/Shell Thickness

        if    (ietype.ne.2.and.ietype.ne.3
     &                    .and. pcomp(text(1),'thic',4)) then
          d(14) = ev(1)
          if(ev(2).eq.0.0d0) then
            d(37) = 5.d0/6.d0
          else
            d(37) = ev(2)
          endif
          d(102) = ev(3)

c       Mass density

        elseif(pcomp(text(1),'dens',4)) then
          d(4)  = ev(1)
          if(ev(2).eq.0.0d0) then
            d(8) = 1.d0
          else
            d(8) = max(0.0d0,ev(2))
          endif

c       Damping factor

        elseif(pcomp(text(1),'damp',4)) then
          if(pcomp(text(2),'rayl',4)) then
            d(77) = ev(1)
            d(78) = ev(2)
          else
            d(70) = ev(1)
          endif

c       Mass matrix type

        elseif(pcomp(text(1),'mass',4)) then
          if(pcomp(text(2),'lump',4) .or.
     &       pcomp(text(2),'diag',4)) then
            d(7) = 0.0d0
          elseif(pcomp(text(2),'cons',4)) then
            d(7) = 1.0d0
          elseif(pcomp(text(2),'off',3)) then
            d(7) = -1.0d0
          else
            d(7) = ev(1)
          endif

c       Method number

        elseif(pcomp(text(1),'meth',4)) then

          d(80) = ev(1)
          d(81) = ev(2)

c       Rigid body specifier

        elseif(pcomp(text(1),'rigi',4)) then
          nrmatl = nint(ev(1))
          nrbody = max(nrbody,nrmatl)

c       Hierarchical flag

        elseif(pcomp(text(1),'hier',4)) then

          if(pcomp(text(2),'off',3)) then
            d(120) = 0.0d0   ! Set formulation to normal
          else
            d(120) = 1.0d0   ! Set formulation to hierarchical
          endif

c       Normal surface loading

        elseif(((ietype.eq.4.or.ietype.eq.5) .or.
     &          (ietype.eq.3.and.ndm.eq.2)) .and.
     &           pcomp(text(1),'load',4)  ) then

          d(10)  = ev(1)
          if(pcomp(text(2),'axia',4) .or. pcomp(text(2),'tang',4)) then
            d(83) = ev(2)
          elseif(pcomp(text(2),'foll',4)) then
            d(68) = 1.0d0
          else
            d(68) = 0.0d0
          endif

c       Solid material body forces

c       elseif(ietype.ne.4.and.ietype.ne.6
c    &                    .and. pcomp(text(1),'body',4)) then
        elseif(pcomp(text(1),'body',4)) then

          if(pcomp(text(2),'heat',4)) then
            d(66) = ev(1)
          else
            d(11) = ev(1)
            d(12) = ev(2)
            d(13) = ev(3)
            if(pcomp(text(2),'loca',4)) then
              d(69) = 2.0d0
            elseif(pcomp(text(2),'gfol',4)) then
              d(69) = 3.0d0
            elseif(pcomp(text(2),'lfol',4)) then
              d(69) = 4.0d0
            else
              d(69) = 1.0d0
            endif
          endif

c       Orthotropic material principal direction angle

        elseif((ietype.eq.1.or.ietype.eq.4.or.ietype.eq.5)
     &                    .and. pcomp(text(1),'angl',4)) then

          if(pcomp(text(2),'pola',4)) then
            d(85) = 1.0d0
            d(86) = ev(1)
            d(87) = ev(2)
            d(88) = ev(3)
          else
            d(31) = ev(1)
          endif

c       Transient indicator for element use

        elseif(pcomp(text(1),'tran',4)) then

          if(pcomp(text(2),'back',4)) then
            d(89) = 1
          elseif(pcomp(text(2),'newm',4)) then
            d(89) = 2
          else
            d(89) = 3
          endif

c       Frame/Truss: Cross section properties

        elseif((ietype.eq.2.or.ietype.eq.3)
     &                    .and. pcomp(text(1),'cros',4)) then

          cflag = .true.

          do i = 1,7
            d(31+i) = ev(i)
          end do ! i

          if(ev(5).eq.0.0d0) then
            d(36) = d(33) + d(34)
          endif
          if(ev(6).eq.0.0d0) then
            d(37) = 5.d0/6.d0
          endif
          if(ev(7).eq.0.0d0) then
            d(38) = 5.d0/6.d0
          endif

c       Frame/Truss: Layer section properties

        elseif((ietype.eq.2 .or. ietype.eq.3)
     &         .and. pcomp(text(1),'laye',4)) then

          if(ndm.eq.2 .or. ietype.eq.2) then
            oned            = .true.
            nlay            =  nlay + 1
            if(nlay.gt.13) then
              write(ilg,4010)
              write(iow,4010)
              call plstop()
            endif
            d(101+2*nlay) =  ev(1)
            d(102+2*nlay) =  ev(2)
            nlob          =  min(6,max(nlob,nint(ev(3))))

c           Set area/inertia to prevent error message

            d(32)           = -999.d0
            d(33)           = -999.d0
          else
            write(ilg,4016)
            write(iow,4016)
            if(ior.lt.0) write(*,4016)
            erflag = .true.
          endif

c       Truss/Frame: Shaped section properties

        elseif((ietype.eq.2 .or. ietype.eq.3) .and.
     &               pcomp(text(1),'sect',4)) then

          oned    = .true.
          if    (pcomp(text(2),'laye',4)) then
            if(ndm.eq.2 .or. ietype.eq.2) then
              nlay            =  nlay + 1
              if(nlay.gt.13) then
                write(ilg,4010)
                write(iow,4010)
                call plstop()
              endif
              d(101+2*nlay) =  ev(1)
              d(102+2*nlay) =  ev(2)
              nlob          =  min(6,max(nlob,nint(ev(3))))
            else
              write(ilg,4016)
              write(iow,4016)
              if(ior.lt.0) write(*,4016)
              erflag = .true.
            endif

          elseif(pcomp(text(2),'tube',4)) then

            ntub   =  max(4,nint(ev(3)))
            nlob   =  min(6,max(nlob,nint(ev(4))))
            d(100) =  1
            d(101) =  ntub         ! number of sector
            d(102) =  nlob         ! quadrature/sector
            d(103) =  ev(1)        ! radius
            d(104) =  ev(2)        ! thickness

          elseif(pcomp(text(2),'rect',4)) then

            nrct            =  nrct + 1
            if(nrct.gt.5) then
              write(ilg,4011)
              write(iow,4011)
              call plstop()
            endif
            n1              =  max(3,min(5,nint(ev(5))))
            n3              =  max(3,min(5,nint(ev(6))))
            nqud            =  nqud + n3*n1
            d(100)          =  2
            d(101+5*nrct-4) =  ev(1)
            d(101+5*nrct-3) =  ev(2)
            d(101+5*nrct-2) =  ev(3)
            d(101+5*nrct-1) =  ev(4)

            d(101+5*nrct  ) =  10*n1 + n3

          elseif(pcomp(text(2),'wide',4)) then
            d(100)  = 3
            nqud    = 12
            do i = 1,6
              d(100+i) = ev(i)
            end do ! i
          elseif(pcomp(text(2),'chan',4)) then
            d(100)  = 4
            nqud    = 12
            do i = 1,6
              d(100+i) = ev(i)
            end do ! i
          elseif(pcomp(text(2),'angl',4)) then
            d(100)  = 5
            nqud    = 8
            do i = 1,6
              d(100+i) = ev(i)
            end do ! i
          elseif(pcomp(text(2),'circ',4)) then
            d(100)  = 6
            i       = nint(ev(2))
            call int2dc(i,nqud,tt)
            do i = 1,6
              d(100+i) = ev(i)
            end do ! i
          else
            write(iow,*) ' *ERROR* - no section ',text(2)
            write(ilg,*) ' *ERROR* - no section ',text(2)
            call plstop()
          endif

c         Set area/inertia to prevent error message

          d(32)           = -999.d0
          d(33)           = -999.d0

c       Piezo-electric model inputs

        elseif(pcomp(text(1),'piez',4)) then

          if(pcomp(text(2),'perm',4)) then
            d(151) = ev(1)
            d(152) = ev(2)
            d(153) = ev(3)
          elseif(pcomp(text(2),'line',4) .or.
     &           pcomp(text(2),'lin+',4)) then
            d(150) = max(d(150),1.d0)
            d(154) = ev(1)
            d(155) = ev(2)
            d(156) = ev(3)
          elseif(pcomp(text(2),'lin-',4)) then
            d(150) = 4.d0
            d(157) = ev(1)
            d(158) = ev(2)
            d(159) = ev(3)
          endif

c       Frame: Reference node/vector set

        elseif(pcomp(text(1),'refe',4)) then

          if    (pcomp(text(2),'node',4)) then
            lref = 1
            do i = 1,ndm
              refx(i) = ev(i)
            end do ! i
          elseif(pcomp(text(2),'vect',4)) then
            lref = 2
            do i = 1,ndm
              refx(i) = ev(i)
            end do ! i
          elseif(pcomp(text(2),'pola',4)) then
            lref = 3
            do i = 1,ndm
              refx(i) = ev(i)
            end do ! i
          elseif(pcomp(text(2),'axia',4)) then
            lref = 4
            do i = 1,ndm
              refx(i) = ev(i)
            end do ! i
          elseif(pcomp(text(2),'shea',4)) then
            sref = 1
            do i = 1,ndm
              tref(i) = ev(i)
            end do ! i
          else
            lref = 0
            write(ilg,4008)
            write(iow,4008)
            if(ior.lt.0) then
              write(*,4008)
            end if
            erflag = .true.
          endif

c       Frame: No Shear Option

        elseif((ietype.eq.3 .or. ietype.eq.4)
     &          .and. pcomp(text(1),'shea',4)) then

          if(pcomp(text(2),'off',3)) then
            d(79) = 1.0d0
          else
            d(79) = 0.0d0
          endif

c       Truss/Frame: Nonlinear flag

        elseif((ietype.eq.2.or.ietype.eq.3)
     &                    .and. pcomp(text(1),'nonl',4)) then

          if(ev(1).eq.0.0d0) then
            d(39) = 1.0d0
          else
            d(39) = 0.0d0
          endif

c       Kinematics type: Small deformation

        elseif(pcomp(text(1),'smal',4)) then

          dtype = 1
          if(ietype.eq.1) then
            sflag = .true.
            fflag = .false.
          endif

c       Kinematics type: Finite deformation

        elseif(pcomp(text(1),'fini',4)) then

          dtype = -1
          if(ietype.eq.1) then
            fflag = .true.
            sflag = .false.
            if(pcomp(text(2),'volu',4)) then
              d(170) = max(1,min(3,nint(ev(1))))
            endif
          endif

c       Element type: Displacement

        elseif(pcomp(text(1),'disp',4)) then

          etype = 1

c       Element type: Mixed (B-Bar)

        elseif(pcomp(text(1),'mixe',4)) then

          etype = 2

          if(pcomp(text(2),'enha',4)) then
            etype = 5
          endif

c       Element type: Enhanced Strain

        elseif(pcomp(text(1),'enha',4)) then

          etype = 3

          if(pcomp(text(2),'mixe',4)) then
            etype = 5
          endif

c       Element type: Energy-Momentum Conserving

        elseif(pcomp(text(1),'cons',4)) then

          etype = 4

c       Element type: Co-rotational Formulation

        elseif(pcomp(text(1),'coro',4)) then

          etype = 6

c       Initial data

        elseif(pcomp(text(1),'init',4)) then

c         Constant initial stress/strain state

          if( pcomp(text(2),'stre',4)) then
            d(160) = 1
            do i = 1,6
              d(160+i) = ev(i)
            end do ! i
          elseif( pcomp(text(2),'stra',4)) then
            d(160) = 3
            do i = 1,6
              d(160+i) = ev(i)
            end do ! i
          elseif( pcomp(text(2),'augm',4)) then
            d(160) = 2
          endif

c       Quadrature data

        elseif((ietype.eq.1.or.ietype.ge.5)
     &              .and. pcomp(text(1),'quad',4)) then

c         Set nodal quadrature flag

          if(pcomp(text(2),'noda',4) .or. pcomp(text(2),'node',4)) then
            d(182) = 1.0d0
          elseif(pcomp(text(2),'loba',4)) then
            d(182) = 0.0d0
            jj     = -1
          elseif(pcomp(text(2),'gaus',4)) then
            d(182) = 0.0d0
          endif

c         Limit quadrature for built-in elements

          if(ev(1).ne.0.0d0 .or. ev(2).ne.0.0d0) then
            if(eltype.gt.0) then
              ii = max(-1,min(5,nint(ev(1))))
              jj = max( 1,min(5,nint(ev(2))))
            else
              ii = nint(ev(1))
              jj = nint(ev(2))
            endif
          endif

c       Temperature dof

        elseif(ietype.ne.6 .and.
     &     (pcomp(text(1),'temp',4) .or. pcomp(text(1),'volt',4))) then

          tdof = ev(1)

c       Input solution type

        elseif((ietype.eq.1 .or. ietype.eq.6)  .and.
     &                     pcomp(text(1),'plan',4)) then
          if( pcomp(text(2),'stre',4)) then
            stype = 1
          else
            stype = 2
          endif

        elseif((ietype.eq.1 .or. ietype.eq.6)  .and.
     &                     pcomp(text(1),'axis',4)) then
          if( pcomp(text(2),'tors',4)) then
            stype = 8
          else
            stype = 3
          endif

c       Penalty parameter

        elseif(pcomp(text(1),'pena',4)) then

          d(60)  = ev(1)
          d(121) = ev(2)

c       Incompressible flag

        elseif(pcomp(text(1),'inco',4)) then

          if(pcomp(text(2),'off',3)) then
            incomp = .false.
          else
            incomp = .true.
          endif

c       Thermal properties

        elseif(ietype.ne.6 .and. pcomp(text(1),'ther',4)) then

          tflag = .true.

c         Orthotropic input

          if(pcomp(text(2),'orth',4)) then

            alp(1) = ev(1)
            alp(2) = ev(2)
            alp(3) = ev(3)
            d(9)   = ev(4)

c         Transversely isotropic input

          elseif(pcomp(text(2),'tran',4)) then

            alp(1) = ev(1)
            alp(2) = ev(2)
            alp(3) = ev(1)
            d(9)   = ev(3)

c         Isotropic inputs

          else

            alp(1) = ev(1)
            alp(2) = ev(1)
            alp(3) = ev(1)
            d(9)   = ev(2)

          endif

c       Input error estimator value

        elseif(pcomp(text(1),'adap',4)) then

          if(pcomp(text(2),'erro',4)) then
            if(ev(1).eq.0.0d0) then
              d(50) = 0.05
            else
              d(50) = ev(1)
            endif
          endif

c       Elastic properties

        elseif(ietype.ne.6 .and. pcomp(text(1),'elas',4)) then

          eflag  = .true.

c         Transverse isotropy inputs

          if(pcomp(text(2),'tran',4)) then

            iflag = .false.
            imat  =  5

            e1   = ev(1)
            e2   = ev(2)
            e3   = ev(1)
            nu12 = ev(3)
            nu23 = ev(3)
            nu31 = ev(4)
            g12  = ev(5)
            g23  = ev(5)
            g31  = ev(1)/(2.d0 + nu31 + nu31)

c         Orthotropic inputs

          elseif(pcomp(text(2),'orth',4)) then

            iflag = .false.
            imat  =  5

            e1   = ev(1)
            e2   = ev(2)
            e3   = ev(3)
            nu12 = ev(4)
            nu23 = ev(5)
            nu31 = ev(6)
            g12  = ev(7)
            g23  = ev(8)
            g31  = ev(9)

c         Finite Elastic Models

c         Regular compressible Neo-Hookean

          elseif(ietype.eq.1 .and. pcomp(text(2),'neoh',4)) then

            imat  =  1
            dtype = -1
            fflag = .true.
            sflag = .false.

            e1    = ev(1)
            nu12  = ev(2)

c         Modified compressible Neo-Hookean

          elseif(ietype.eq.1 .and. pcomp(text(2),'mneo',4)) then

            imat  =  2
            dtype = -1
            fflag = .true.
            sflag = .false.

            e1    = ev(1)
            nu12  = ev(2)

c         Ogden compressible model

          elseif(ietype.eq.1 .and. pcomp(text(2),'ogde',4)) then

            imat  =  3
            dtype = -1
            fflag = .true.
            sflag = .false.

            bulk  = ev(1)
            g12   = ev(2)
            nu12  = ev(3)
            g23   = ev(4)
            nu23  = ev(5)
            g31   = ev(6)
            nu31  = ev(7)
            e1    = max(bulk,abs(g12),abs(g23),abs(g31))

c         Regular compressible Mooney-Rivlin

          elseif(ietype.eq.1 .and. pcomp(text(2),'moon',4) .or.
     &                             pcomp(text(2),'mmoo',4)) then

            if(pcomp(text(2),'moon',4)) then
              imat  =  9
            else
              imat  =  10
            endif
            dtype = -1
            fflag = .true.
            sflag = .false.

            e1    = ev(1)
            nu12  = ev(2)
            nu23  = ev(3)

c         Fung type exponential model

          elseif(pcomp(text(2),'fung',4)) then

            imat   =  7
            dtype  = -1
            fflag  = .true.
            sflag  = .false.

            d(30) = ev(1)
            do i = 1,9
              d(20+i) = ev(i+1)
            end do ! i
            e1     = d(30)

c         Full inputs

          elseif(pcomp(text(2),'modu',4) .or.
     &           pcomp(text(2),'comp',4)) then

            imat  =  8

c           Input by rows (store transpose)

            nsiz   = nint(ev(1))
            if(nsiz.eq.0 .or. nsiz.gt.6) then
              write(ilg,4006) text(2),nsiz
              write(iow,4006) text(2),nsiz
              call plstop()
            endif
            d(200) = nsiz
            do i =  1,nsiz
              errck = pinput(dd(1,i),nsiz)
            end do ! i

c           Convert to moduli

            if(pcomp(text(2),'comp',4)) then
              do i = 1,nsiz
                do j = 1,i
                  dd(i,j) = dd(j,i)
                end do ! j
              end do ! i
              call invert(dd,nsiz,6)
            endif

c           Save triangular part

            n1 = 200
            do i = 1,nsiz
              do j = 1,i
                n1 = n1 + 1
                d(n1) = dd(j,i)
              end do ! j
            end do ! i
            nsiz = n1
            e1   = d(201)

c         Compressible Arruda-Boyce model

          elseif(ietype.eq.1 .and. pcomp(text(2),'arru',4)) then

            imat  = 11
            fflag = .true.
            sflag = .false.
            dtype = -1

            e1    = ev(1)
            nu12  = ev(2)
            nu23  = ev(3)

c         Compressible Yeoh model

          elseif(ietype.eq.1 .and. pcomp(text(2),'yeoh',4)) then

            imat  = 12
            fflag = .true.
            sflag = .false.
            dtype = -1

            e1    = ev(1)
            nu12  = ev(2)
            nu23  = ev(3)
            nu31  = ev(4)

c         Isotropic inputs

          else

            imat  = 4
            if(pcomp(text(2),'stve',4).or.pcomp(text(2),'stvk',4)) then
              imat  = 5
              dtype = -1
              fflag = .true.
              sflag = .false.
            elseif(pcomp(text(2),'cons',4)) then
              imat  = 6
              dtype = -1
              etype =  4
              fflag = .true.
              sflag = .false.
            endif
            iflag = .true.

            e1    = ev(1)
            e2    = e1
            e3    = e1
            nu12  = ev(2)

            nu23  = nu12
            nu31  = nu12
            g12   = 0.5d0*e1/(1.d0 + nu12)
            g23   = g12
            g31   = g12
            bulk  = e1/(1.d0 - 2.d0*min(0.49999999d0,nu12))/3.d0

          endif

c       Tension/Compression only properties

        elseif(ietype.eq.2 .and. pcomp(text(1),'tens',4)) then
          d(167) = 1.d0
        elseif(ietype.eq.2 .and. pcomp(text(1),'comp',4)) then
          d(167) = -1.d0

c       Activation indicators

        elseif(pcomp(text(1),'acti',4)) then

          if(pcomp(text(1),'ther',4)) then
            if(ev(1).eq.0.0d0) then
              d(168) =  ev(1)
            else
              d(168) = -1.0d0
            endif
          else ! mechanical
            if(ev(1).eq.0.0d0) then
              d(169) =  ev(1)
            else
              d(169) = -1.0d0
            endif
          endif

c       Viscoelastic-damgage properties

        elseif(ietype.ne.6 .and. pcomp(text(1),'dama',4)) then

          damgfl = .true.
          d(58)  = ev(1)
          d(59)  = ev(2)

c       Viscoelastic properties

        elseif(ietype.ne.6 .and. pcomp(text(1),'visc',4)) then

          viscfl = .true.
          plasfl = .false.
          d(40)  = 2.d0  ! Constitution is now Viscoelastic
          nvis   = nvis + 1
          if(nvis.le.3) then
            d(50+2*nvis-1) = ev(1)
            d(50+2*nvis  ) = ev(2)
            d(57)          = nvis
          else
            write(ilg,4015)
            write(iow,4015)
            call plstop()
          endif

c       Plasticity properties

        elseif(ietype.ne.6 .and. pcomp(text(1),'plas',4)) then

          plasfl = .true.
          viscfl = .false.
          d(40)  = 1.d0  ! Constitution is now Plastic

c         Mises yield/loading function

          if(pcomp(text(2),'mise',4)) then

            d(41) = ev(1)
            d(42) = ev(2)
            d(43) = ev(3)
            if(d(41).eq.0.0d0) d(41) = 1.d+20
            if(d(42).eq.0.0d0) d(42) = d(41)
            if(d(43).eq.0.0d0) d(43) = 1.d0
            d(46) = 1.d0  ! Mises flag

c         Drucker-Prager yield/loading function

          elseif(pcomp(text(2),'druc',4)) then

            d(41) = ev(1)
            d(42) = ev(2)
            if(d(42).eq.0.0d0) d(42) = d(41)
            d(46) = 2.d0  ! Drucker-Prager flag

c         Drucker-Prager yield/loading function

          elseif(pcomp(text(2),'lode',4)) then

            d(41) = ev(1)
            d(42) = ev(2)
            if(d(42).eq.0.0d0) d(42) = d(41)
            d(46) = 3.d0  ! Drucker-Prager flag

c         Hardening parameters

          elseif(pcomp(text(2),'hard',4)) then

            d(44)  = ev(1)
            d(45)  = ev(2)

c         Linear segment hardening parameters

          elseif(pcomp(text(2),'segm',4)) then

            nseg          = nseg + 1
            if(nseg.gt.6) then
              write(ilg,4012)
              write(iow,4012)
              call plstop()
            endif
            d(130)        = nseg
            d(128+3*nseg) = ev(1)
            d(129+3*nseg) = ev(2)
            d(130+3*nseg) = ev(3)

c         Viscoplastic rate parameter

          elseif(pcomp(text(2),'rate',4) .or.
     &           pcomp(text(2),'time',4)) then

            d(180) = ev(1)
            d(181) = max(1.d0,ev(2))

c        Generalized modal parameters

          elseif(pcomp(text(2),'gene',4)) then

            d(40)  = 3.d0  ! Constitution is now Generalized Plastic
            d(41)  = ev(1)
            d(42)  = ev(2)
            d(43)  = ev(3)

          endif

c       Angular velocity: rad/sec

        elseif(ietype.ne.6 .and. pcomp(text(1),'omeg',4)) then

        if(pcomp(text(2),'cycl',4)) then
          d(65) = ev(1)*2.d0*pi
        else
          d(65) = ev(1)
        endif

c       Fourier Heat Conduction properties

        elseif(pcomp(text(1),'four',4)) then

          hflag = .true.
          d(67) = 1.0d0 ! Heat constitution added

c         Surface confection (2-d only)

          if(pcomp(text(2),'conv',4)) then

            d(127) = ev(1)    ! Surface convection      (h)
            d(128) = ev(2)    ! Free stream temperature (T_inf)

c         Reference absolute temperature

          elseif(pcomp(text(2),'temp',4)) then

            d(129) = ev(1)    ! Reference absolute temperature (T_ref)

c         Orthotropic inputs

          elseif(pcomp(text(2),'orth',4)) then

            d(61) = ev(1)     ! k_1 - conductivity
            d(62) = ev(2)     ! k_2 - conductivity
            d(63) = ev(3)     ! k_3 - conductivity
            d(64) = ev(4)     ! c   - specific heat

c         Isotropic inputs

          else

            d(61) = ev(1)
            d(62) = ev(1)
            d(63) = ev(1)
            d(64) = ev(2)

          endif

c         Thermal dof for mechanical solid and truss problems

          if(ietype.eq.1 .and. tdof.eq.0 .or. ietype.eq.2) then

            tdof = ndm + 1

          endif

c       Ground motion acceleration factors/proportional load numbers

        elseif(ietype.ne.6 .and. pcomp(text(1),'grou',4)) then

          do i = 1,ndm
            d(70+i) = ev(2*i-1)
            d(73+i) = ev(2*i)
          end do ! i

c       Element (Lagrange multiplier) unknowns

        elseif(pcomp(text(1),'elem',4).or.pcomp(text(1),'lagr',4)) then

          d(230) = ev(1)

c       Consititutive solution start state: Default = elastic (0)

        elseif(pcomp(text(1),'star',4)) then

          d(84)  = ev(1)

c       User Material Parameters

        elseif(pcomp(text(1),'upar',4)) then

          uparfl = .true.
          do i = 1,10
            d(230+i) = ev(i)
          end do ! i

c       User Material Model interface

        elseif(pcomp(text(1),'ucon',4).or.pcomp(text(1),'fcon',4)) then

c         Default user constitutive equation number

          umat    = 0
          uprm    = ndd - nud
          n1      = 0
          n3      = 0

          uct     = 'read'
          call uconst(text(2),ev, d(1), d(uprm+1),n1,n3, umat)

          if(pcomp(text(1),'ucon',4)) then
            d(uprm) = umat + 100
          else
            d(uprm) = umat + 200
            dtype   = -1
          endif

c         Activate user program models

          uflag  = .true.
          if(e1.eq.0.0d0) then
            if(ietype.eq.3) then
              write(ilg,6000)
              write(iow,6000)
              call plstop()
            endif
            e1 = 1.0d0
          endif

c         Increase number of history terms/quadrature point

          nh1    = nh1 + n1
          nh3    = nh3 + n3

c       Check end of data

        elseif(pcomp(text(1),'    ',4)) then

c         Transfer to sets and checks

          inputs = .false.

        endif

      end do ! while

c     Number of stress/strain history terms/pt

      if(ndm.eq.3 .or. stype.eq.8) then
        ntm = 6
      elseif(ndm.eq.2) then
        ntm = 4
      else
        ntm = 1
      endif

c     PART 3: Set final parameters and output material state

c     Set moduli

      if(ietype.ne.6) then


c       Small deformation options

        if(sflag) then

          if(imat.ne.8) then
            d(1)    = e1
            d(2)    = nu12
            d(3)    = alp(1)

            dd(1,1) =  1.d0/e1
            dd(2,2) =  1.d0/e2
            dd(3,3) =  1.d0/e3

            dd(1,2) = -dd(1,1)*nu12
            dd(1,3) = -dd(3,3)*nu31
            dd(2,3) = -dd(2,2)*nu23

c           1-Dimensional Models

            if(stype.eq.5) then
              dd(2,2) =  1.d0
              dd(3,3) =  1.d0
              dd(1,2) =  0.d0
              dd(1,3) =  0.d0
              dd(2,3) =  0.d0

c           Plane Stress Models

            elseif(stype .eq.4) then

              d(90)   =  dd(1,3)
              d(91)   =  dd(2,3)

              dd(3,3) =  1.d0
              dd(2,3) =  0.d0
              dd(1,3) =  0.d0

            endif

            dd(2,1) = dd(1,2)
            dd(3,1) = dd(1,3)
            dd(3,2) = dd(2,3)

c           Mechanical modulus properties

            if(.not.incomp) then
              call invert(dd,3,6)
              if(min(dd(1,1),dd(2,2),dd(3,3)).lt.0) then
                write(iow,4005) ((dd(i,j),j=1,3),i=1,3)
                write(ilg,4005) ((dd(i,j),j=1,3),i=1,3)
                if(ior.lt.0) then
                  write(*,4005) ((dd(i,j),j=1,3),i=1,3)
                endif
              endif
            elseif(iflag) then
              dd(1,1) =  4.d0*g12/3.d0
              dd(2,2) =  dd(1,1)
              dd(3,3) =  dd(1,1)
              dd(1,2) = -dd(1,1)*0.5d0
              dd(2,3) =  dd(1,2)
              dd(1,3) =  dd(1,2)
            endif

c           Set for plane stress

            if(stype.eq.5) then
              dd(3,3) = 0.0d0
            endif

c           Save moduli

            d(21) = dd(1,1)
            d(22) = dd(2,2)
            d(23) = dd(3,3)
            d(24) = dd(1,2)
            d(25) = dd(2,3)
            d(26) = dd(1,3)
            d(27) = g12
            d(28) = g23
            d(29) = g31

c           Thermal properties

            d(47) = dd(1,1)*alp(1) + dd(1,2)*alp(2) + dd(1,3)*alp(3)
            d(48) = dd(2,1)*alp(1) + dd(2,2)*alp(2) + dd(2,3)*alp(3)
            d(49) = dd(3,1)*alp(1) + dd(3,2)*alp(2) + dd(3,3)*alp(3)

c           Set for plane stress problems

            if(stype .eq.4) then
              d(23) = 0.0d0
              d(49) = 0.0d0
              d(92) = alp(3)
            endif
          endif ! imat.ne.8

c         Output parameters for element

          if(plasfl) then
            jj   = ii
            d(6) = jj
          endif

c         Output elastic properties

          if(eflag) then
            if(stype.eq.5 .or. iflag) then
              write(iow,2000) wd(stype),e1,nu12
              if(ior.lt.0) then
                write(*,2000) wd(stype),e1,nu12
              endif
            elseif(imat.eq.8) then
              write(iow,2063) (d(i),i=201,nsiz)
              if(ior.lt.0) then
                write(*,2063) (d(i),i=201,nsiz)
              endif
            else
              write(iow,2001) wd(stype),e1,e2,e3,nu12,nu23,nu31,
     &                        g12,g23,g31,d(31)
              if(ior.lt.0) then
                write(*,2001) wd(stype),e1,e2,e3,nu12,nu23,nu31,
     &                        g12,g23,g31,d(31)
              endif
              if(nint(d(85)).gt.0) then
                write(iow,2065) (d(i),i=86,85+ndm)
                if(ior.lt.0) then
                  write(*,2065) (d(i),i=86,85+ndm)
                endif
              endif
            endif
            if(ietype.eq.1 .or.ietype.eq.5) then
              if(d(182).gt.0.0d0) then
                i = 2
              else
                i = 1
              endif
              write(iow,2019) qd(i),ii,jj
              if(ior.lt.0) then
                write(*,2019) qd(i),ii,jj
              endif
            endif
          elseif(.not.uflag) then
            write(ilg,4000)
            write(iow,4000)
            if(ior.lt.0) write(*,4000)
            erflag = .true.
          endif

c         Output thermal expansions

          if(tflag) then
            write(iow,2002) alp,d(9),tdof
            if(ior.lt.0) then
              write(*,2002) alp,d(9),tdof
            endif
            if(d(129).ne.0.0d0) then
              write(iow,2076) d(129)
              if(ior.lt.0) then
                write(*,2076) d(129)
              endif
            endif
          endif

c         Output fourier heat conduction properties

          if(hflag) then
            write(iow,2020) wd(stype),(d(i),i=61,64),d(66)
            if(d(127).gt.0.0d0) then
              write(iow,2073) d(127),d(128)
            endif
            if(ior.lt.0) then
              write(*,2020) wd(stype),(d(i),i=61,64),d(66)
              if(d(127).gt.0.0d0) then
                write(*,2073) d(127),d(128)
              endif
            endif
          endif

c       Finite deformation options

        elseif(fflag)then

c         Output Regular NeoHookean

          if(imat.eq.1) then

            bulk = e1/(1.d0 - nu12*2.d0)/3.d0
            g12  = e1/(1.d0 + nu12)/2.d0
            if(ior.lt.0) then
              write(*,2010) ' ',wd(stype),e1,nu12,bulk,g12
            endif
            write(iow,2010) ' ',wd(stype),e1,nu12,bulk,g12
            if(ior.lt.0) then
              write(*,2062) nint(d(170))
            endif
            write(iow,2062) nint(d(170))

c           Compute Lame' parameters

            d(1)  = e1
            d(2)  = nu12
            d(3)  = alp(1)

            d(21) = bulk - 2.d0/3.d0*g12
            d(22) = g12

c         Output Modified NeoHookean

          elseif(imat.eq.2) then

            bulk = e1/(1.d0 - nu12*2.d0)/3.d0
            g12  = e1/(1.d0 + nu12)/2.d0
            if(ior.lt.0) then
              write(*,2010) ' Modified ',wd(stype),e1,nu12,bulk,g12
            endif
            write(iow,2010) ' Modified ',wd(stype),e1,nu12,bulk,g12
            if(ior.lt.0) then
              write(*,2062) nint(d(170))
            endif
            write(iow,2062) nint(d(170))

c           Set material parameters

            d(1)  = e1
            d(2)  = nu12
            d(3)  = alp(1)

            d(21) = bulk
            d(22) = g12

c         Output Ogden

          elseif(imat.eq.3) then

            nn = 1
            if(nu23.ne.0.0d0) nn = 2
            if(nu31.ne.0.0d0) nn = 3
            if(nn.eq.1) then
              if(ior.lt.0) then
                write(*,2024) bulk,g12,nu12
              endif
              write(iow,2024) bulk,g12,nu12
            elseif(nn.eq.2) then
              if(ior.lt.0) then
                write(*,2024) bulk,g12,nu12,g23,nu23
              endif
              write(iow,2024) bulk,g12,nu12,g23,nu23
            elseif(nn.eq.3) then
              if(ior.lt.0) then
                write(*,2024) bulk,g12,nu12,g23,nu23,g31,nu31
              endif
              write(iow,2024) bulk,g12,nu12,g23,nu23,g31,nu31
            endif
            if(ior.lt.0) then
              write(*,2062) nint(d(170))
            endif
            write(iow,2062) nint(d(170))

c           Set material parameters

            d(21) = bulk
            d(22) = g12
            d(23) = nu12
            d(24) = g23
            d(25) = nu23
            d(26) = g31
            d(27) = nu31
            d(28) = nn

c         Logarithmic stretch model: Like linear elastic

          elseif(imat.eq.4) then

            if(ior.lt.0) then
              write(*,2027) e1,nu12,bulk,g12
            endif
            write(iow,2027) e1,nu12,bulk,g12

c         Set material parameters

            d(1)  = e1
            d(2)  = nu12
            d(3)  = alp(1)

            d(21) = bulk
            d(22) = g12

c         St. Venant-Kirchhoff

          elseif(imat.eq.5 .or. imat.eq.6) then

            if(iflag) then
              write(iow,2000) wd(stype),e1,nu12
              if(ior.lt.0) then
                write(*,2000) wd(stype),e1,nu12
              endif
            else
              write(iow,2001) wd(stype),e1,e2,e3,nu12,nu23,nu31,
     &                        g12,g23,g31,d(31)
              if(ior.lt.0) then
                write(*,2001) wd(stype),e1,e2,e3,nu12,nu23,nu31,
     &                        g12,g23,g31,d(31)
              endif
              if(d(85).gt.0) then
                write(iow,2065) (d(i),i=1,86,85+ndm)
                if(ior.lt.0) then
                  write(*,2065) (d(i),i=1,86,85+ndm)
                endif
              endif
            endif

c           Set material parameters

            d(1)    = e1
            d(2)    = nu12
            d(3)    = alp(1)

            dd(1,1) =  1.d0/e1
            dd(2,2) =  1.d0/e2

            dd(1,2) = -dd(1,1)*nu12
            dd(2,1) =  dd(1,2)

            if(stype.ne.1) then
              dd(3,3) =  1.d0/e3

              dd(1,3) = -dd(3,3)*nu31
              dd(2,3) = -dd(2,2)*nu23

              dd(3,1) =  dd(1,3)
              dd(3,2) =  dd(2,3)
            else
              dd(3,3) =  1.d0

              dd(1,3) =  0.d0
              dd(2,3) =  0.d0

              dd(3,1) =  0.d0
              dd(3,2) =  0.d0
            endif

c           Mechanical modulus properties

c           if(.not.incomp .or. stype.eq.1) then
            if(.not.incomp .or. stype.eq.4) then
              call invert(dd,3,6)
            endif

c           Save moduli

            d(21) = dd(1,1)
            d(22) = dd(2,2)
            d(23) = dd(3,3)
            d(24) = dd(1,2)
            d(25) = dd(2,3)
            d(26) = dd(3,1)
            d(27) = g12
            d(28) = g23
            d(29) = g31

c           Thermal properties

            d(47) = dd(1,1)*alp(1) + dd(1,2)*alp(2) + dd(1,3)*alp(3)
            d(48) = dd(2,1)*alp(1) + dd(2,2)*alp(2) + dd(2,3)*alp(3)
            if(stype.ne.1) then
              d(49) = dd(3,1)*alp(1) + dd(3,2)*alp(2) + dd(3,3)*alp(3)
            else
              d(49) = 0.d0
            endif

c         Fung model

          elseif(imat.eq.7) then
            if(ietype.eq.5) then
              write(iow,2054) wd(stype),(d(20+i),i=0,4)
              if(ior.lt.0) then
                write(*,2054) wd(stype),(d(20+i),i=0,4)
              endif
            else
              write(iow,2055) wd(stype),(d(20+i),i=0,9)
              if(ior.lt.0) then
                write(*,2055) wd(stype),(d(20+i),i=0,9)
              endif
            endif

c         Output Regular Mooney-Rivlin

          elseif(imat.eq.9) then

            bulk = e1/(1.d0 - nu12*2.d0)/3.d0
            g12  = e1/(1.d0 + nu12)/2.d0
            if(ior.lt.0) then
              write(*,2074) ' ',wd(stype),e1,nu12,bulk,g12,nu23
            endif
            write(iow,2074) ' ',wd(stype),e1,nu12,bulk,g12,nu23
            if(ior.lt.0) then
              write(*,2062) nint(d(170))
            endif
            write(iow,2062) nint(d(170))

c           Compute Lame' parameters

            d(1)  = e1
            d(2)  = nu12
            d(3)  = alp(1)

            d(21) = bulk - 2.d0/3.d0*g12
            d(22) = g12
            d(23) = nu23

c         Output Modified Mooney-Rivlin

          elseif(imat.eq.10) then

            d(21) = e1/(1.d0 - nu12*2.d0)/3.d0
            d(22) = e1/(1.d0 + nu12)/2.d0
            d(23) = nu23
            if(ior.lt.0) then
              write(*,2074) ' Modified ',wd(stype),e1,nu12,
     &                      (d(i),i=21,23)
            endif
            write(iow,2074) ' Modified ',wd(stype),e1,nu12,
     &                      (d(i),i=21,23)
            if(ior.lt.0) then
              write(*,2062) nint(d(170))
            endif
            write(iow,2062) nint(d(170))

            d(1)  = e1
            d(2)  = nu12
            d(3)  = alp(1)

c         Output compressible Arruda-Boyce

          elseif(imat.eq.11) then

            bulk = e1/(1.d0 - nu12*2.d0)/3.d0
            g12  = e1/(1.d0 + nu12)/2.d0
            if(ior.lt.0) then
              write(*,2077) ' Modified ',wd(stype),e1,nu12,nu23,bulk,g12
            endif
            write(iow,2077) ' Modified ',wd(stype),e1,nu12,nu23,bulk,g12
            if(ior.lt.0) then
              write(*,2062) nint(d(170))
            endif
            write(iow,2062) nint(d(170))

c           Set material parameters

            d(1)  = e1
            d(2)  = nu12
            d(3)  = alp(1)

            d(21) = bulk
            d(22) = g12/(1+ 0.6d0*nu23 + 99.d0*nu23*nu23/175.d0)
            d(23) = nu23

c         Output compressible Yeoh

          elseif(imat.eq.12) then

            bulk = e1/(1.d0 - nu12*2.d0)/3.d0
            g12  = e1/(1.d0 + nu12)/2.d0
            if(ior.lt.0) then
              write(*,2078) ' Modified ',wd(stype),e1,nu12,nu23,bulk,g12
            endif
            write(iow,2078) ' Modified ',wd(stype),e1,nu12,nu23,bulk,g12
            if(ior.lt.0) then
              write(*,2062) nint(d(170))
            endif
            write(iow,2062) nint(d(170))

c           Set material parameters

            d(1)  = e1
            d(2)  = nu12
            d(3)  = alp(1)

            d(21) = bulk
            d(22) = g12
            d(23) = nu23
            d(24) = nu31

          endif
        endif

c       Output tension/compression only states

        if(d(167).gt.0.0d0) then
          if(ior.lt.0) then
            write(*,2052)
          endif
          write(iow,2052)
        elseif(d(167).lt.0) then
          if(ior.lt.0) then
            write(*,2053)
          endif
          write(iow,2053)
        endif

c       Output activation indicators

        if(d(168).ne.0.0d0) then
          if(ior.lt.0) then
            write(*,2070)
          endif
          write(iow,2070)
        endif
        if(d(169).ne.0) then
          if(ior.lt.0) then
            write(*,2071)
          endif
          write(iow,2071)
        endif

c       Output Transient type

        if(d(89).gt.0.0d0) then
          if(ior.lt.0) then
            write(*,2072) nint(d(89))
          endif
          write(iow,2072) nint(d(89))
        endif

c       Output shell/plate thickness

        if(stype.ne.5 .and. stype.ne.7) then
          write(iow,2018) d(14)
          if(nint(d(102)).gt.1) write(iow,2051) nint(d(102))
        endif
        if(stype.eq.4) then
          write(iow,2021) d(10)
          if(dtype.lt.0 .and. d(68).gt.0.0d0) then
            write(iow,2033)
          endif
        endif
        if(ior.lt.0) then
          if(stype.ne.5 .and. stype.ne.7) then
            write(*,2018) d(14)
            if(nint(d(102)).gt.1) write(*,2051) nint(d(102))
          endif
          if(stype.eq.4) then
            write(*,2021) d(10)
            if(dtype.lt.0 .and. d(68).gt.0.0d0) then
              write(iow,2033)
            endif
          endif
        endif

c       Output density and body loading

        if(ior.lt.0) then
          write(*,2029) d(4),d(11),d(12),d(13)
        endif
        write(iow,2029) d(4),d(11),d(12),d(13)

c       Output constant initial stresses

        if(nint(d(160)).eq.1) then
          if(ietype.eq.2 .or. ietype.eq.3) then
            j = 1
          else
            j = 6
          endif
          write(iow,2044) (d(160+i),i=1,j)
          if(ior.lt.0) then
            write(*,2044) (d(160+i),i=1,j)
          endif
        elseif(nint(d(160)).eq.3) then
          if(ietype.eq.2 .or. ietype.eq.3) then
            j = 1
          else
            j = 6
          endif
          write(iow,2045) (d(160+i),i=1,j)
          if(ior.lt.0) then
            write(*,2045) (d(160+i),i=1,j)
          endif
        elseif(nint(d(160)).eq.2) then
          if(ietype.eq.2 .or. ietype.eq.3) then
            naug = 1
            write(iow,2014)
            if(ior.lt.0) then
              write(*,2014)
            endif
          endif
        endif

c       Output frame distributed loading type

        if(d(69).eq.1.0d0) then
          write(iow,2034)
          if(ior.lt.0) write(*,2034)
        elseif(d(69).eq.2.0d0) then
          write(iow,2035)
          if(ior.lt.0) write(*,2035)

        elseif(d(69).eq.3.0d0) then
          write(iow,2036)
          if(ior.lt.0) write(*,2036)
        elseif(d(69).eq.4.0d0) then
          write(iow,2037)
          if(ior.lt.0) write(*,2037)
        endif

c       Output angular velocity

        if(d(65).ne.0.0d0) then
          if(ior.lt.0) then
            write(*,2030) d(65)
          endif
          write(iow,2030) d(65)
        endif

c       Output ground acceleration factors

        flg = .false.
        efl = .false.
        do i = 1,ndm
          if(nint(d(73+i)).gt.0 ) flg = .true.
          if(nint(d(73+i)).gt.10) then
            efl = .true.
          endif
        end do ! i
        if(flg) then
          write(iow,2023) (d(70+i),nint(d(73+i)),i=1,ndm)
          if(ior.lt.0) then
            write(*,2023) (d(70+i),nint(d(73+i)),i=1,ndm)
          endif
        endif

        if(efl) then
          write(ilg,4007)
          write(iow,4007)
          if(ior.lt.0) write(*,4007)
          erflag = .true.
        endif

c       Output section properties

        if(cflag) then
          if(ietype.eq.2) then
            j = 32
          else
            j = 38
          endif
          write(iow,2015) (d(i),i=32,j)
          if(ior.lt.0) then
            write(*,2015) (d(i),i=32,j)
          endif

          if(ndm.eq.2 .and. ietype.eq.3) then
            if(d(33).lt.0.0d0) then
              write(ilg,4009) d(33)
              write(iow,4009) d(33)
              if(ior.lt.0) then
                write(*,4009) d(33)
              endif
              erflag = .true.
            endif
          elseif(ndm.eq.3 .and. ietype.eq.3) then
            ev(1) = d(33)*d(34) - d(35)*d(35)
            ev(2) = 1.d-8*max(abs(d(33)),abs(d(34)))
            if((ev(1).lt.ev(2) .and. d(35).ne.0.0d0)
     &                          .or. ev(1).le.0.0d0) then
              write(ilg,4009) ev(1)
              write(iow,4009) ev(1)
              if(ior.lt.0) then
                write(*,4009) ev(1)
              endif
              erflag = .true.
            endif
          endif
        endif

c       Output layer data

        if(nlay.gt.0) then
          d(101) = nlay
          d(102) = nlob
          nqud   = nlay*nlob
          if(ior.lt.0) then
            write(*,2038) (i,d(101+2*i),d(101+2*i+1),i=1,nlay)
            write(*,2039) nlob
          endif
          write(iow,2038) (i,d(101+2*i),d(101+2*i+1),i=1,nlay)
          write(iow,2039) nlob
          if(d(37).eq.0.0d0) d(37) = 5.d0/6.d0

c         Compute location of centroid and shift input data

c         call bmcent(nlay,d(102))

c         Output shifted layer data

c         if(ior.lt.0) then
c           write(*,2038) (i,d(101+2*i),d(101+2*i+1),i=1,nlay)
c           write(*,2039) nlob
c         endif
c         write(iow,2038) (i,d(101+2*i),d(101+2*i+1),i=1,nlay)
c         write(iow,2039) nlob
        endif

c       Set default shear factor for 3-d beams

        if(nint(d(100)) .gt. 0) then
          if(d(37).eq.0.0d0) d(37) = 5.d0/6.d0
          if(d(38).eq.0.0d0) d(38) = 5.d0/6.d0
        endif

c       Output tube data

        if(ntub.gt.0) then
          if(ior.lt.0) then
            write(*,2040) d(103),d(104),ntub,nlob
          endif
          write(iow,2040) d(103),d(104),ntub,nlob
        endif

c       Output rectangle data

        if(nrct.gt.0) then
          d(101) = nrct
          if(ior.lt.0) then
            write(*,2042) ((d(96+5*j+i),i=1,4),nint(d(101+5*j))/10,
     &                      mod(nint(d(101+5*j)),10),j=1,nrct)
          endif
          write(iow,2042) ((d(96+5*j+i),i=1,4),nint(d(101+5*j))/10,
     &                      mod(nint(d(101+5*j)),10),j=1,nrct)
        endif

c       Output shaped section data

        if(nint(d(100)).eq.3) then   ! Wide flange
          write(iow,2048) (d(100+i),i=1,6)
          if(ior.lt.0) then
            write(*,2048) (d(100+i),i=1,6)
          endif
        elseif(nint(d(100)).eq.4) then   ! Channel
          write(iow,2049) (d(100+i),i=1,6)
          if(ior.lt.0) then
            write(*,2049) (d(100+i),i=1,6)
          endif
        elseif(nint(d(100)).eq.5) then   ! Angle
          write(iow,2050) (d(100+i),i=1,4)
          if(ior.lt.0) then
            write(*,2050) (d(100+i),i=1,4)
          endif
        elseif(nint(d(100)).eq.6) then   ! Solid circle
          d(102) = max(1.d0,min(4.d0,d(102)))  ! quadrature: 1 le nqd le 4
          write(iow,2059) d(101),nint(d(102))
          if(ior.lt.0) then
            write(*,2059) d(101),nint(d(102))
          endif
        endif

c       Output plasticity parameters

        if(plasfl) then

c         Small deformation plasticity (also for one-d problems)

          if(sflag .or. oned) then

            if(nint(d(40)).eq.1) then
              write(iow,2003) (d(i),i=41,45)
              if(ior.lt.0) then
                write(*,2003) (d(i),i=41,45)
              endif
              if(d(180).gt.0.0d0) then
                write(iow,2075) d(180),nint(d(181))
                if(ior.lt.0) then
                  write(*,2075) d(180),nint(d(181))
                endif
              endif
              nn  = 1
              nhd = ntm + 1  ! epp, ep(i),i=1,ntm
            elseif(nint(d(40)).eq.3) then
              write(iow,2047) (d(i),i=41,45)
              if(ior.lt.0) then
                write(*,2047) (d(i),i=41,45)
              endif
              nn  = 2
              nhd = ntm + 1  ! epp, ep(i),i=1,ntm
            endif

c         Finite deformation plasticity

          elseif(fflag) then

            if(imat.eq.4) then
              nhd = 2 + 6 + 3 ! epp, dete, be(i),i=1,6, e_pl(i),i=1,3
              nn  = 1
              if(nint(d(46)).eq.1) then
                sigt  = d(41)
                sigc  = d(42)
                d(41) = sqrt(2.d0/3.d0)*sigt
                d(42) = sqrt(2.d0/3.d0)*sigc
                write(iow,2056) sigt,sigc,d(43), d(44),d(45)
              elseif(nint(d(46)).eq.2) then
                sigt  = d(41)
                sigc  = d(42)
                d(41) = sqrt( 8.d0/3.d0)*(sigc * sigt)/(sigt + sigc)
                d(42) = sqrt(18.d0/3.d0)*(sigc - sigt)/(sigt + sigc)
                write(iow,2057) sigt,sigc,d(41),d(42),d(44),d(45)
              elseif(nint(d(46)).eq.3) then
                sigt  = d(41)
                sigc  = d(42)
                d(41) = sqrt( 8.d0/3.d0)*(sigc * sigt)/(sigt + sigc)
                d(42) =                  (sigc - sigt)/(sigt + sigc)
                write(iow,2058) sigt,sigc,d(41),d(42),d(44),d(45)
                if(abs(d(42)).gt. 0.125d0) then
                  write(ilg,4013)
                  write(iow,4013)
                  if(ior.lt.0) then
                    write(*,4013)
                  endif
                endif

              endif
            else
              write(ilg,4014)
              write(iow,4014)
              call plstop()
            endif

          endif ! fflag/sflag

c         One dimensional model history storage

          if(oned) then
            if(nseg.eq.0) then
              nh1 = nh1 + nn + 3   ! ep(1); epp; state each layer
            else
              write(iow,2041) (d(i),i=131,130+3*nseg)
              if(ior.lt.0) then
                write(*,2041) (d(i),i=131,130+3*nseg)
              endif
              nh1 = nh1 + 4   ! ep(1); alp(1); epp; state each layer
            endif

c         Multi dimensional model history storage

          else
            if(ndm.eq.3 .or. stype.eq.8) then
              nh1 = nh1 + nhd*nn + 1
            else
              if(stype.eq.1 .and. sflag) then
                nh1 = nh1 + nn + 7   ! ep(3), beta(3), ep; state (plane stress)
              else
                if(fflag) then
                  nh1 = nh1 + nhd*nn + 1
                else
                  nh1 = nh1 + 5*nn + 1 ! ep(4); ep; state        (plane strain)
                endif
              endif
            endif
          endif

        endif ! plasfl

c       Output visco-elastic properties

        if(viscfl) then
          write(iow,2004) (d(i),i=51,50+2*nvis)
          if(ior.lt.0) then
            write(*,2004) (d(i),i=51,50+2*nvis)
          endif
          if(fflag) then
            i = 1
          else
            i = 0
          endif
          if(oned) then
            nh1 = nh1 + i + 1 + nvis   ! xi; eps_n; hh(*) each point
          else
            if(ndm.eq.3 .or. stype.eq.8) then
              nh1 = nh1 + i + 6 + 6*nvis  ! xi; eps_n(6); hh(6) each point
            else
              nh1 = nh1 + i + 4 + 4*nvis  ! xi; eps_n(4); hh(4) each point
            endif
          endif
        endif ! viscfl

c       Output damage properties

        if(damgfl) then
          write(iow,2028) d(58),d(59)
          if(ior.lt.0) then
            write(*,2028) d(58),d(59)
          endif
        endif ! damgfl

c       Constitutive start (.ne.0 for nonclassical elastic)

        if(d(84).gt.0) then

          write(iow,2066) d(84)
          if(ior.lt.0) then
            write(*,2066) d(84)
          endif

        endif

c       Kinematics type

        if(ietype.ne.6) then
          if(dtype.gt.0) then
            write(iow,2005)
            if(ior.lt.0) then
              write(*,2005)
            endif
          else
            write(iow,2006)
            if(ior.lt.0) then
              write(*,2006)
            endif
          endif
        endif

c       Hierarchical formulation

        if(d(120).gt.0) then

          write(iow,2081)
          if(ior.lt.0) then
            write(*,2081)
          endif

        endif

c       Output Piezo-electric properties

        if(d(150).gt.0.0d0) then
          if(tdof.eq.0) tdof = ndm + 1
          write(iow,2043) tdof,(d(i),i=151,155 + nint(d(150)))
          if(ior.lt.0) then
            write(*,2043) tdof,(d(i),i=151,155 + nint(d(150)))
          endif
        endif

c       Element type

        if(ietype.eq.1) then
          if(etype.eq.1) then
            write(iow,2012) 'Displacement'
            if(ior.lt.0) then
              write(*,2012) 'Displacement'
            endif
          elseif(etype.eq.2) then
            write(iow,2012) 'Mixed'
            if(ior.lt.0) then
              write(*,2012) 'Mixed'
            endif
          elseif(etype.eq.3) then
            write(iow,2012) 'Enhanced'
            if(ior.lt.0) then
              write(*,2012) 'Enhanced'
            endif
          elseif(etype.eq.4) then
            write(iow,2012) 'Energy Conserving'
            if(ior.lt.0) then
              write(*,2012) 'Energy Conserving'
            endif
          elseif(etype.eq.5) then
            write(iow,2012) 'Mixed-Enhanced'
            if(ior.lt.0) then
              write(*,2012) 'Mixed-Enhanced'
            endif
          elseif(etype.eq.6) then
            write(iow,2012) 'Co-Rotational'
            if(ior.lt.0) then
              write(*,2012) 'Co-Rotational'
            endif
          endif
          if(incomp) then
            write(iow,2067)
            if(ior.lt.0) then
              write(*,2067)
            endif
          endif
        elseif(ietype.eq.3) then
          if(dtype.gt.0) then
            if(nint(d(79)).eq.0.0d0) then
              write(iow,2012) 'Linear Displacement'
              if(ior.lt.0) then
                write(*,2012) 'Linear Displacement'
              endif
            else
              write(iow,2012) 'Cubic Displacement'
              if(ior.lt.0) then
                write(*,2012) 'Cubic Displacement'
              endif
            endif
          else
            if(nint(d(79)).eq.0.0d0) then
              write(iow,2012) 'Linear Displacement'
              if(ior.lt.0) then
                write(*,2012) 'Linear Displacement'
              endif
              if(ndm.eq.3) then
                if(etype.eq.2 .or. etype.eq.3) then
                  write(iow,2013) 'Incremental Rotation'
                  if(ior.lt.0) then
                    write(*,2013) 'Incremental Rotation'
                  endif
                elseif(etype.eq.4) then
                  write(iow,2013) 'Energy Conserving'
                  if(ior.lt.0) then
                    write(*,2013) 'Energy Conserving'
                  endif
                elseif(etype.eq.6) then
                  write(iow,2013) 'Co-Rotational'
                  if(ior.lt.0) then
                    write(*,2013) 'Co-Rotational'
                  endif
                else
                  write(iow,2013) 'Simo-Vu Quoc Displacement'
                  if(ior.lt.0) then
                    write(*,2013) 'Simo-Vu Quoc Displacement'
                  endif
                endif
              endif
            else
              write(iow,2012) 'Cubic Displacement, 2nd Order'
              if(ior.lt.0) then
                write(*,2012) 'Cubic Displacement, 2nd Order'
              endif
            endif
          endif
        endif

c     Thermal element only

      elseif(ietype.eq.6) then

        write(iow,2020) wd(stype),(d(i),i=61,64),d(66),d(4)
        if(d(127).gt.0.0d0) then
          write(iow,2073) d(127),d(128)
        endif
        if(d(182).gt.0.0d0) then
          j = 2
        else
          j = 1
        endif
        write(iow,2019) qd(j),ii,jj
        if(ior.lt.0) then
          write(*,2020) wd(stype),(d(i),i=61,64),d(66),d(4)
          if(d(127).gt.0.0d0) then
            write(*,2073) d(127),d(128)
          endif
          write(*,2019) qd(j),ii,jj
        endif

      endif

c     Shear deformation on/off

      if(ietype.eq.3 .or. ietype.eq.4) then
        if(d(79).gt.0.0d0) then
          write(iow,2013) 'No Shear Deformation'
          if(ior.lt.0) then
            write(*,2013) 'No Shear Deformation'
          endif
        else
          write(iow,2013) 'Includes Shear Deformation',d(37)
          if(ior.lt.0) then
            write(*,2013) 'Includes Shear Deformation',d(37)
          endif
        endif
      endif

c     Mass type

      if(d(4).gt.0.0d0) then
        if(d(7).eq.0.0d0) then
          write(iow,2007)
          if(ior.lt.0) then
            write(*,2007)
          endif
        elseif(d(7).eq.1.0d0) then
          write(iow,2008)
          if(ior.lt.0) then
            write(*,2008)
          endif
        else
          write(iow,2009) d(7)
          if(ior.lt.0) then
            write(*,2009) d(7)
          endif
        endif
        if(ietype.eq.3 .or. ietype.eq.5.and.dtype.lt.0) then
          write(iow,2032) d(8)
          if(ior.lt.0) then
            write(*,2032) d(8)
          endif
        endif
      endif

c     Damping factor

      if(d(70).gt.0.0d0) then
        write(iow,2046) d(70)
        if(ior.lt.0) then
          write(*,2046) d(70)
        endif
      endif

c     Rayleigh Damping factors

      if(max(abs(d(77)),abs(d(78))).gt.0.0d0) then
        write(iow,2060) d(77),d(78)
        if(ior.lt.0) then
          write(*,2060) d(77),d(78)
        endif
      endif

c     Number of element variables

      if(nint(d(230)).ne.0) then
        write(iow,2069) nint(d(230))
        if(ior.lt.0) then
          write(iow,2069) nint(d(230))
        endif
      endif

c     Method type factors

      if(max(abs(d(80)),abs(d(81))).gt.0.0d0) then
        write(iow,2061) d(80),d(81)
        if(ior.lt.0) then
          write(*,2061) d(80),d(81)
        endif
      endif

c     Error value

      if(d(50).ne.0.0d0) then
        write(iow,2011) d(50)
        if(ior.lt.0) then
          write(*,2011) d(50)
        endif
      endif

      if(d(39).ne.0.0d0) then

        write(iow,2016)
        if(ior.lt.0) then
          write(*,2016)
        endif
      endif

c     Output penalty value

      if(d(60).ne.0.0d0) then
        write(iow,2022) d(60)
        if(ior.lt.0) then
          write(*,2022) d(60)
        endif
      endif

c     Set beam/shell kappa

      if((ietype.eq.4 .or. ietype.eq.5)      .and.
     &   (nint(d(79)).eq.0) .and. dtype.lt.0) then
        write(iow,2017) d(37)
        if(ior.lt.0) then
          write(*,2017) d(37)
        endif
      endif

c     Output reference coordinate/vector

      if    (sref.eq.1) then

        d(93) = sref
        do i = 1,2
          d(93+i) = tref(i)
        end do ! i

        write(iow,2031) (i,tref(i),i=1,2)
        if(ior.lt.0) then
          write(*,2031) (i,tref(i),i=1,2)
        endif

      endif


c     Output reference coordinate/vector

      if(lref.gt.0) then
        d(96) = lref
        do i = 1,3
          d(96+i) = refx(i)
        end do ! i
      endif
      if    (lref.eq.1) then

        write(iow,2025) (i,refx(i),i=1,ndm)
        if(ior.lt.0) then
          write(*,2025) (i,refx(i),i=1,ndm)
        endif

      elseif(lref.eq.2) then

        write(iow,2026) (i,refx(i),i=1,ndm)
        if(ior.lt.0) then
          write(*,2026) (i,refx(i),i=1,ndm)
        endif

      elseif(lref.eq.3) then

        write(iow,2079)
        if(ior.lt.0) then
          write(*,2079)
        endif

      elseif(lref.eq.4) then

        write(iow,2064)
        if(ior.lt.0) then
          write(*,2064)
        endif

      endif

c     Output rigid body indicator

      if(nrmatl.gt.0) then
        write(iow,2080) nrmatl
        if(ior.lt.0) then
          write(*,2080) nrmatl
        endif
      endif

c     Output user parameters

      if(uparfl) then

        j = 0
        do i = 1,10
          if(d(230+i).ne.0.0d0) then
            j = i
          endif
        end do ! i

        write(iow,2068) (i,d(230+i),i=1,j)
        if(ior.lt.0) then
          write(*,2068) (i,d(230+i),i=1,j)
        endif

      endif

c     Save quadrature

      d( 5) = ii
      d( 6) = jj

c     Save types

      istp  = stype
      d(15) = nh1 + naug
      d(16) = stype
      d(17) = etype
      d(18) = dtype
      d(19) = tdof
      d(20) = imat
      d(31) = d(31)*rad

c     Solid (1) or shell/membrane (5)

      if(ietype.eq.1 .or. ietype.eq.5) then
        i     = nh1
        if(ndm.eq.2) then
          nh1   = nh1*ii*ii
          if(etype.eq.5) nh1 = nh1 + 2 ! mixed-enhance unknown
        elseif(ndm.eq.3) then
          nh1   = nh1*ii*ii*ii
          if(etype.eq.5) nh1 = nh1 + 9 ! mixed-enhance unknown
        endif
        if(etype.eq.3 .and. sflag) nh1 = nh1 + 5     ! linear element
        if(etype.eq.3 .and. fflag) nh1 = nh1 + i     ! extra quadrature
        if(stype.eq.1)             nh1 = nh1 + ii*ii ! store F_33/q_pt

c     Truss (2) or Frame (3)

      elseif(ietype.eq.2 .or. ietype.eq.3) then
        d(82) = max(1,nqud)
        if(nint(d(79)).eq.1) then
          nh1 = nh1*3
        endif
        if(oned.and.nlay.gt.1) then
          nh1 = nh1*((nlay-1)*(nlob-1) + 1)  ! Total data/cross section
        elseif(oned.and.ntub.gt.1) then
          nh1 = nh1*ntub*(nlob-1)            ! Total data/cross section
        elseif(oned.and.nqud.gt.1) then
          nh1 = nh1*nqud                     ! Total data/cross section
        endif

c     Plate (4)

      elseif(ietype.eq.4) then
        nh1   = nh1*ii

c     Three dimensional (7)

      elseif(ietype.eq.7) then
        nh1   = nh1*ii**3
      endif

c     Set history for saving element variables

      d(149) = nh1 ! Total number of variables on frame cross section
      nh3    = nh3 + ndv

c     Set number of element variables

      nlm  = nint(d(230))

c     Check for warnings

      if(d(4).eq.0.0d0) then
        write(iow,3000)
        if(ior.lt.0) then
          write(*,3000)
        endif
      endif

c     Check for errors

      if     (ietype.eq.1) then
        if(e1.eq.0.0d0) then
          write(ilg,4001)
          write(iow,4001)
          if(ior.lt.0) then
            write(*,4001)
          endif
        endif
      elseif (ietype.eq.2) then
        if(e1.eq.0.0d0 .or. d(32).eq.0.0d0) then
          write(ilg,4002) 'Truss',e1,d(32)
          write(iow,4002) 'Truss',e1,d(32)
          if(ior.lt.0) then
            write(*,4002) 'Truss',e1,d(32)
          endif
        endif
      elseif (ietype.eq.3) then
        if(e1.eq.0.0d0 .or. d(32).eq.0.0d0 .or. d(33).eq.0.0d0) then
          write(ilg,4002) 'Frame',e1,d(32),d(33)
          write(iow,4002) 'Frame',e1,d(32),d(33)
          if(ior.lt.0) then
            write(*,4002) 'Frame',e1,d(32),d(33)
          endif
        endif
      elseif (ietype.eq.4) then
        if(e1.eq.0.0d0 .or. d(37).eq.0.0d0 .or. d(14).eq.0.0d0) then
          write(ilg,4003) 'Plate',e1,d(14),d(37)
          write(iow,4003) 'Plate',e1,d(14),d(37)
          if(ior.lt.0) then
            write(*,4003) 'Plate',e1,d(14),d(37)
          endif
        endif
      elseif (ietype.eq.5) then
        if(e1.eq.0.0d0 .or. d(37).eq.0.0d0 .or. d(14).eq.0.0d0) then
          write(ilg,4003) 'Shell',e1,d(14),d(37)
          write(iow,4003) 'Shell',e1,d(14),d(37)
          if(ior.lt.0) then
            write(*,4003) 'Shell',e1,d(14),d(37)
          endif
        endif
      elseif (ietype.eq.6) then
        if(d(61).eq.0.0d0) then
          write(ilg,4004) e1,d(61)
          write(iow,4004) e1,d(61)
          if(ior.lt.0) then
            write(*,4004) e1,d(61)
          endif
        endif
      endif

c     Output user history information

      if(uflag) then
        write(iow,5000) nint(d(15)),nh1,nh3
      endif

c     Errors detected in inputs

      if(erflag) call plstop()

c     I/O Formats

2000  format( 5x,'M e c h a n i c a l   P r o p e r t i e s'//
     & 10x,a,' Analysis'//
     & 10x,'Modulus E       ',1p,1e12.5/
     & 10x,'Poisson ratio   ',0p,1f8.5/)

2001  format( 5x,'M e c h a n i c a l   P r o p e r t i e s'//
     & 10x,a,' Analysis'//
     & 10x,'Modulus E-1     ',1p,1e12.5/
     & 10x,'Modulus E-2     ',1p,1e12.5/
     & 10x,'Modulus E-3     ',1p,1e12.5/
     & 10x,'Poisson ratio 12',0p,1f8.5 /
     & 10x,'Poisson ratio 23',0p,1f8.5 /
     & 10x,'Poisson ratio 31',0p,1f8.5 /
     & 10x,'Modulus G-12    ',1p,1e12.5/
     & 10x,'Modulus G-23    ',1p,1e12.5 /
     & 10x,'Modulus G-31    ',1p,1e12.5/
     & 10x,'Angle (psi)   ',0p,1f10.5/)

2002  format(/5x,'T h e r m a l   E x p a n s i o n s'//
     & 10x,'Th. Alpha-1',1p,1e17.5/10x,'Th. Alpha-2',1p,1e17.5/
     & 10x,'Th. Alpha-3',1p,1e17.5/10x,'T_0        ',1p,1e17.5/
     & 10x,'Th. D.O.F. ',i7/)

2003  format(/5x,'M i s e s   P l a s t i c   P a r a m e t e r s'//
     &  10x,'Yield stress short  ',1p,1e15.5/
     &  10x,'Yield stress infin. ',1p,1e15.5/
     &  10x,'Hardening exponent  ',1p,1e15.5/
     &  8x,'Linear hardening parts'/
     &  10x,'Isotropic hardening ',1p,1e15.5/
     &  10x,'Kinematic hardening ',1p,1e15.5/)

2004  format(/5x,'V i s c o e l a s t i c   P a r a m e t e r s'//
     &   (10x,'nu-i :ratio',1p,1e15.5/10x,'lam-i:time ',1p,1e15.5:))

2005  format( 10x,'Formulation : Small deformation.')
2006  format( 10x,'Formulation : Finite deformation.')

2007  format( 10x,'Mass type   : Lumped.')
2008  format( 10x,'Mass type   : Consistent.')
2009  format( 10x,'Mass type   : Interpolated:',1p,1e11.4)

2010  format( 5x,'M e c h a n i c a l   P r o p e r t i e s'//
     &        4x,a,'NeoHookean Stored Energy Function '//
     &       10x,a,' Analysis'//
     &       10x,'Modulus E       ',1p,1e12.5/
     &       10x,'Poisson ratio   ',0p,1f8.5/
     &       10x,'Bulk Modulus    ',1p,1e12.5/
     &       10x,'Shear Modulus   ',1p,1e12.5/1x)

2011  format(10x,'Error estimator ',1p,1e12.5 )

2012  format(10x,'Element type: ',a)

2013  format( 24x,a:/,24x,'Shear factor (k) =',1p,1e12.5)

2014  format(/10x,'Augmented Solution for Inextensible Behavior')

2015  format(/5x,'C r o s s   S e c t i o n   P a r a m e t e r s'//
     &    10x,'Area      ',1p,1e15.5:/10x,'I_xx      ',1p,1e15.5/
     &    10x,'I_yy      ',1p,1e15.5 /10x,'I_xy      ',1p,1e15.5/
     &    10x,'J_zz      ',1p,1e15.5 /10x,'Kappa_x   ',1p,1e15.5/
     &    10x,'Kappa_y   ',1p,1e15.5 / )

2016  format( 10x,'Non-linear analysis')

2017  format( 10x,'Plate/Shell : Kappa ',1p,1e15.5)

2018  format(10x,'Thickness       ',1p,1e12.5)
2019  format(5x,a,'   Q u a d r a t u r e'/
     & 10x,'Quadrature: Arrays',i3    /10x,'Quadrature: Output',i3/)

2020  format(/5x,'T h e r m a l   P r o p e r t i e s'//
     & 10x,a,' Analysis'//
     & 10x,'Cond. K-1    ',1p,1e15.5/ 10x,'Cond. K-2    ',1p,1e15.5/
     & 10x,'Cond. K-3    ',1p,1e15.5/ 10x,'Specific Heat',1p,1e15.5/
     & 10x,'Heat Source  ',1p,1e15.5/:10x,'Density      ',1p,1e15.5)

2021  format( 10x,'Loading - q ',1p,1e16.5)

2022  format( 10x,'Penalty - k ',1p,1e16.5)

2023  format(
     &    /5x,'P r o p o r t i o n a l   B o d y   L o a d i n g s',//
     &    10x,'1-Dir. Factor',1p,1e15.5,': Proportional Load No.',i3/:
     &    10x,'2-Dir. Factor',1p,1e15.5,': Proportional Load No.',i3/:
     &    10x,'3-Dir. Factor',1p,1e15.5,': Proportional Load No.',i3/)

2024  format( 5x,'M e c h a n i c a l   P r o p e r t i e s'//
     &        5x,'Ogden Stored Energy Function '//
     &       10x,'Bulk  Modulus   ',1p,1e12.5//
     &       10x,'Modulus  - C1   ',1p,1e12.5/
     &       10x,'Exponent - n1   ',1p,1e12.5/:/
     &       10x,'Modulus  - C2   ',1p,1e12.5/
     &       10x,'Exponent - n2   ',1p,1e12.5/:/
     &       10x,'Modulus  - C3   ',1p,1e12.5/
     &       10x,'Exponent - n3   ',1p,1e12.5/1x)

2025  format(/5x,'R e f e r e n c e    C o o r d i n a t e s'//
     &        (10x,'X-',i1,' = ',1p,1e12.5:))

2026  format(/5x,'R e f e r e n c e    V e c t o r'//
     &        (10x,'v-',i1,' = ',1p,1e12.5:))

2027  format( 5x,'M e c h a n i c a l   P r o p e r t i e s'//
     &        5x,'Logarithmic Stretch Stored Energy Function '//
     &       10x,'Modulus E       ',1p,1e12.5/
     &       10x,'Poisson ratio   ',0p,1f8.5/
     &       10x,'Bulk Modulus    ',1p,1e12.5/
     &       10x,'Shear Modulus   ',1p,1e12.5/1x)

2028  format(/5x,'D a m a g e    P a r a m e t e r s'//
     &    10x,'Damage Limit',1p,1e13.5/10x,'Decay rate',1p,1e15.5)

2029  format(10x,'Density         ',1p,1e12.5//
     &       10x,'1-Gravity Load  ',1p,1e12.5/
     &       10x,'2-Gravity Load  ',1p,1e12.5/
     &       10x,'3-Gravity Load  ',1p,1e12.5/1x)

2030  format(/10x,'Angular Velocity (radians/time)',1p,1e15.5/)

2031  format(/5x,'S h e a r   C e n t e r   C o o r d i n a t e s'//
     &        (10x,'X-',i1,' = ',1p,1e12.5:))

2032  format( 10x,'Rotational Mass Factor:',1p,1e12.5)

2033  format( 15x,'Follower loading.')
2034  format( 10x,'Global axis loading.')
2035  format( 10x,'Local axis loading.')
2036  format( 10x,'Global axis follower loading.')
2037  format( 10x,'Local axis follower loading.')

2038  format(/5x,'B e a m   L a y e r   D a t a'/
     &       10x,'Layer  Z-location        Width'/(i14,1p2e13.4))

2039  format(10x,'Number of Lobatto points/layer =',i3)

2040  format(/5x,'T u b u l a r   B e a m   D a t a'/
     &       10x,'Radius            =',1p,e13.4/
     &       10x,'Thickness         =',1p,e13.4/
     &       10x,'No. Sectors       =',i6/
     &       10x,'Quad. Pts/Sectors =',i6/)

2041  format(/5x,'H a r d e n i n g    P a r a m e t e r s'//
     &    12x,'Strain e_p  Isotropic Yield-Y  Kinematic Hardening-H'/
     &     (11x,1p,1e12.5,1p,1e19.5,1p,1e22.5))

2042  format(/5x,'R e c t a n g u l a r   B e a m   D a t a'/
     &       10x,'  y-low left  z-low left  y-up right  z-up right',
     &           '  y-quadr  z-quadr'/(10x,1p,4e12.4,2i9))

2043  format(/5x,'P i e z o - E l e c t r i c   D a t a'//
     &       10x,'Voltage D.O.F.    =',i3//
     &       10x,'Permeability -  1 =',1p,e13.4/
     &       10x,'Permeability -  2 =',1p,e13.4/
     &       10x,'Permeability -  3 =',1p,e13.4//
     &       10x,'Coupling(+)  e_12 =',1p,e13.4/
     &       10x,'Coupling(+)  e_22 =',1p,e13.4/
     &       10x,'Coupling(+)  e_45 =',1p,e13.4/:
     &       10x,'Coupling(-)  e_12 =',1p,e13.4/
     &       10x,'Coupling(-)  e_22 =',1p,e13.4/
     &       10x,'Coupling(-)  e_45 =',1p,e13.4/)

2044  format(/5x,'I n i t i a l   S t r e s s   D a t a'//
     &       10x,'11-Stress         =',1p,e13.4/:
     &       10x,'22-Stress         =',1p,e13.4/
     &       10x,'33-Stress         =',1p,e13.4/
     &       10x,'12-Stress         =',1p,e13.4/
     &       10x,'23-Stress         =',1p,e13.4/
     &       10x,'31-Stress         =',1p,e13.4/)

2045  format(/5x,'I n i t i a l   S t r a i n   D a t a'//
     &       10x,'11-Strain         =',1p,e13.4/:
     &       10x,'22-Strain         =',1p,e13.4/
     &       10x,'33-Strain         =',1p,e13.4/
     &       10x,'12-Strain         =',1p,e13.4/
     &       10x,'23-Strain         =',1p,e13.4/
     &       10x,'31-Strain         =',1p,e13.4/)

2046  format( 10x,'Damping     ',1p,1e16.5)

2047  format(/5x,'M i s e s   G e n e r a l i z e d   ',
     &           'P l a s t i c   P a r a m e t e r s'//
     &       10x,'Yield stress        ',1p,1e15.5,' (Sigma_0)'/
     &       10x,'Yield stress infin. ',1p,1e15.5,' (Sigma_inf)'/
     &       10x,'Transition parameter',1p,1e15.5,' (Delta)'/
     &        8x,'Linear hardening parts'/
     &       10x,'Isotropic hardening ',1p,1e15.5,' (H_iso)'/
     &       10x,'Kinematic hardening ',1p,1e15.5,' (H_kin)'/)

2048  format(/5x,'W i d e f l a n g e   B e a m   D a t a'/
     &       10x,'Height             ',1p,1e15.5/
     &       10x,'Top    flange width',1p,1e15.5/
     &       10x,'Bottom flange width',1p,1e15.5/
     &       10x,'Top    flange thick',1p,1e15.5/
     &       10x,'Bottom flange thick',1p,1e15.5/
     &       10x,'Web thickness      ',1p,1e15.5/)

2049  format(/5x,'C h a n n e l   B e a m   D a t a'/
     &       10x,'Height             ',1p,1e15.5/
     &       10x,'Top    flange width',1p,1e15.5/
     &       10x,'Bottom flange width',1p,1e15.5/
     &       10x,'Top    flange thick',1p,1e15.5/
     &       10x,'Bottom flange thick',1p,1e15.5/
     &       10x,'Web thickness      ',1p,1e15.5/)

2050  format(/5x,'A n g l e   B e a m   D a t a'/
     &       10x,'Height             ',1p,1e15.5/
     &       10x,'Width              ',1p,1e15.5/
     &       10x,'Height thickness   ',1p,1e15.5/
     &       10x,'Width  thickness   ',1p,1e15.5/)

2051  format( 10x,'Thickness pts',i12)

2052  format( 10x,'Tension only material')

2053  format( 10x,'Compression only material')

2054  format( 5x,'M e c h a n i c a l   P r o p e r t i e s'//
     & 10x,a,' Analysis'//5x,'Fung Exponential Model parameters:'/
     & 10x,'Modulus C   ',    1p,1e16.5/
     & 10x,'A-11        ',    1p,1e16.5/
     & 10x,'A-22        ',    1p,1e16.5/
     & 10x,'A-12        ',    1p,1e16.5/
     & 10x,'A-44        ',    1p,1e16.5/)

2055  format( 5x,'M e c h a n i c a l   P r o p e r t i e s'//
     & 10x,a,' Analysis'//5x,'Fung Exponential Model parameters:'/
     & 10x,'Modulus C   ',    1p,1e16.5/
     & 10x,'A-11        ',    1p,1e16.5/
     & 10x,'A-22        ',    1p,1e16.5/
     & 10x,'A-33        ',    1p,1e16.5/
     & 10x,'A-12        ',    1p,1e16.5/
     & 10x,'A-23        ',    1p,1e16.5/
     & 10x,'A-31        ',    1p,1e16.5/
     & 10x,'A-44        ',    1p,1e16.5/
     & 10x,'A-55        ',    1p,1e16.5/
     & 10x,'A-66        ',    1p,1e16.5/)

2056  format(/8x,'von Mises Yield function'/
     &       10x,'Yield stress      (initial) :  ',1p,e12.5/
     &       10x,'Yield stress     (infinity) :  ',1p,e12.5/
     &       10x,'Hardening exponent  (delta) :  ',1p,e12.5/
     &        8x,'Linear hardening parts'/
     &       10x,'Isotropic hardening (H_iso) :  ',1p,1e12.5/
     &       10x,'Kinematic hardening (H_kin) :  ',1p,1e12.5/)

2057  format(/8x,'Drucker-Prager Yield function'/
     &       10x,'Yield stress in tension     :  ',1p,e12.5/
     &       10x,'Yield stress in compression :  ',1p,e12.5/
     &       10x,'Yield function radius       :  ',1p,e12.5/
     &       10x,'Alpha parameter             :  ',1p,e12.5/
     &        8x,'Linear hardening parts'/
     &       10x,'Isotropic hardening (H_iso) :  ',1p,1e12.5/
     &       10x,'Kinematic hardening (H_kin) :  ',1p,1e12.5/)

2058  format(/8x,'Drucker-Lode Yield function'/
     &       10x,'Yield stress in tension     :  ',1p,e12.5/
     &       10x,'Yield stress in compression :  ',1p,e12.5/
     &       10x,'Yield function radius       :  ',1p,e12.5/
     &       10x,'Lode angle parameter        :  ',1p,e12.5/
     &        8x,'Linear hardening parts'/
     &       10x,'Isotropic hardening (H_iso) :  ',1p,1e12.5/
     &       10x,'Kinematic hardening (H_kin) :  ',1p,1e12.5/)

2059  format(/5x,'C i r c u l a r   B e a m   D a t a'/
     &       10x,'Radius             ',1p,1e15.5/
     &       10x,'Quadrature order   ',i9/
     &       15x,'1 =  4 point (4-interior)'/
     &       15x,'2 =  5 point (4-perimeter, center)'/
     &       15x,'3 =  9 point (4-perimeter, 4-interior, center)'/
     &       15x,'4 = 17 point (8-perimeter, 8-interior, center)')

2060  format(/8x,'Rayleigh Damping Ratios'/
     &       10x,'Mass  value: a0',1p,1e14.5/
     &       10x,'Stiff value: a1',1p,1e14.5)

2061  format(/8x,'Method Type Values'/
     &       10x,'Type 1: Value  ',1p,1e14.5/
     &       10x,'Type 2: Value  ',1p,1e14.5)

2062  format(10x,'Volume Model    ',i2/
     &       15x,'1: U(J) = 0.25*(J**2 - 1 - 2 * ln J)'/
     &       15x,'2: U(J) = 0.50*(J - 1)**2'/
     &       15x,'3: U(J) = 0.50*(ln J)**2'/1x)

2063  format( 5x,'M e c h a n i c a l   P r o p e r t i e s'//
     &       10x,' Moduli From Full Inputs'//
     &       10x,1p,1e11.3/10x,1p,2e11.3/10x,1p,3e11.3/
     &       10x,1p,4e11.3/10x,1p,5e11.3/10x,1p,6e11.3/)

2064  format(/10x,'Axial Reference Vector')

2065  format(/5x,'P o l a r   A n g l e   O r i g i n'/
     &        10x,'1-Coord   = ',1p,1e12.5/
     &        10x,'2-Coord   = ',1p,1e12.5:/
     &        10x,'3-Coord   = ',1p,1e12.5:)

2066  format(/5x,'C o n s t i t u t i v e    S t a r t    S t a t e'//
     &        (10x,'Value    = ',1p,1e12.5:))

2067  format( 10x,'Incompressible Formulation')

2068  format( 5x,'U s e r   P a r a m e t e r s'//
     &      (10x,'Parameter',i3,'   =',1p,e13.4:))

2069  format( 5x,'E l e m e n t   V a r i a b l e s'//
     &       10x,'Number/element    =',i3/)

2070  format( 10x,'Activation thermal')
2071  format( 10x,'Activation mechanical')

2072  format( 10x,'Transient type =',i5/
     &        15x,' 1 = Euler   Method'/
     &        15x,' 2 = Newmark Method'/
     &        15x,' 3 = User    Method')

2073  format(/5x,'S u r f a c e   C o n v e c t i o n'/
     &       10x,'Convection  (h)    ',1p,1e15.5/
     &       10x,'Temperature (T_inf)',1p,1e15.5)

2074  format( 5x,'M e c h a n i c a l   P r o p e r t i e s'//
     &        4x,a,'Mooney-Rivlin Stored Energy Function '//
     &       10x,a,' Analysis'//
     &       10x,'Modulus E       ',1p,1e12.5/
     &       10x,'Poisson ratio   ',0p,1f8.5/
     &       10x,'Bulk Modulus    ',1p,1e12.5/
     &       10x,'Shear Modulus   ',1p,1e12.5/
     &       10x,'2nd Invariant c ',1p,1e12.5/1x)

2075  format( 8x,'Viscoplastic parameters'/
     &       10x,'Rate time parameter ',1p,1e15.5/
     &       10x,'Yield function power',i8/)

2076  format( 8x,'Absolute Reference Temperature'//
     &       10x,'Temperature (T_ref)',1p,1e15.5)

2077  format( 5x,'M e c h a n i c a l   P r o p e r t i e s'//
     &        4x,a,'Arruda-Boyce Stored Energy Function '//
     &       10x,a,' Analysis'//
     &       10x,'Modulus E     - E ',1p,1e12.5/
     &       10x,'Poisson ratio - v ',0p,1f8.5/
     &       10x,'Reciprocal n  - m ',1p,1e12.5/
     &       10x,'Bulk Modulus  - K ',1p,1e12.5/
     &       10x,'Shear Modulus - G ',1p,1e12.5/1x)

2078  format( 5x,'M e c h a n i c a l   P r o p e r t i e s'//
     &        4x,a,'Yeoh Stored Energy Function '//
     &       10x,a,' Analysis'//
     &       10x,'Modulus E     - E ',1p,1e12.5/
     &       10x,'Poisson ratio - v ',0p,1f8.5/
     &       10x,'Factor 1      - k1',0p,1f8.5/
     &       10x,'Factor 2      - k2',0p,1f8.5/
     &       10x,'Bulk Modulus  - K',1p,1e12.5/
     &       10x,'Shear Modulus - G',1p,1e12.5/1x)

2079  format(/5x,'A x i a l    R e f e r e n c e    V e c t o r'/)

2080  format(/10x,'Material Set is Rigid Body Number',i5)

2081  format( 8x,'Hierarchical formulation')

3000  format(/10x,'Material density is zero.')

4000  format(/' *ERROR* INPT: No elastic properties input.'/)

4001  format(/' *ERROR* INPT: Solid Element: Modulus zero.')

4002  format(/' *ERROR* INPT: ',a,': Modulus  = ',1p,1e12.5/
     &        '                      Area     = ',1p,1e12.5/:
     &        '                      Intertia = ',1p,1e12.5/)

4003  format(/' *ERROR* INPT: ',a,': Modulus  = ',1p,1e12.5/
     &        '                      Thickness= ',1p,1e12.5/
     &        '                      Kappa    = ',1p,1e12.5/)

4004  format(/' *ERROR* INPT: Thermal: Conductivity zero.')

4005  format(/' *ERROR* INPT: Elastic properties not postive definite'/
     &        '         Computed Moduli'/(10x,1p,3e15.5))

4006  format(/' *ERROR* INPT: Incorrect size for ',a,' =',i5/)

4007  format(/' *ERROR* INPT: Incorrect proportional load number'/)

4008  format(/' *ERROR* INPT: Incorrect reference vector or node',
     &        ' specification.'/)

4009  format(/' *ERROR* INPT: Inertia determinant zero or negative'/
     &       15x,'Det = I_xx*I_yy - I_xy*I_xy =',1p,1e12.5/)

4010  format(/' *ERROR* INPT: Too many layer sections: Limit = 13.'/)

4011  format(/' *ERROR* INPT: Too many rectangular sections:',
     &        ' Limit = 5.'/)

4012  format(/' *ERROR* INPT: Too many hardening segments: Limit = 6.'/)

4013  format(/' *WARNING* Non-convex yield function'/)

4014  format(/' *ERROR* INPT: Incorrect elastic model for plasticity',
     &       ' solution: Use ISOTropic'/)

4015  format(/' *ERROR* INPT: Too many viscoelastic terms: Limit = 3.'/)

4016  format(/' *ERROR* INPT: Layers are only for 2-d problems.'/)

c     User function Formats

5000  format(//5x,'U s e r   H i s t o r y   I n f o r m a t i o n'//
     &        10x,'Number of history variables/quad point    ',i8/
     &        10x,'Number of history variables total         ',i8/
     &        10x,'Number of element history variables total ',i8/)

6000  format(/' *ERROR* INPT: Frame element with user model must input'/
     &        '               ELAStic ISOTropic E nu before UCON'/)

      end
