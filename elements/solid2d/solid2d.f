c$Id: solid2d.f,v 1.1 2006/11/20 20:32:56 rlt Exp $
      subroutine solid2d(d,ul,xl,ix,tl,s,p,ndf,ndm,nst,isw)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Plane and axisymmetric linear thermo-elastic element
c               driver

c      Inputs:
c         d(*)  - Element parameters
c         ul(ndf,*) - Current nodal solution parameters
c         xl(ndm,*) - Nodal coordinates
c         ix(*)     - Global nodal connections
c         tl(*)     - Nodal temp vector
c         ndf       - Degree of freedoms/node
c         ndm       - Mesh coordinate dimension
c         nst       - Element array dimension
c         isw       - Solution option switch

c      Outputs:
c         s(nst,*)  - Element array
c         p(*)      - Element vector
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'cdata.h'
      include  'eldata.h'
      include  'eqsym.h'
      include  'evdata.h'
      include  'hdata.h'
      include  'iofile.h'
      include  'part0.h'
      include  'part1.h'
      include  'pmod2d.h'
      include  'strnum.h'
      include  'comblk.h'

      logical   mech,ther,mther, errck
      integer   ndf,ndm,nst,isw, i,tdof,pdof, ix(*)
      real*8    d(*),ul(ndf,16),xl(ndm,*),tl(*),s(nst,nst),p(nst)
      real*8    shp(3,16),th(16)

      save

c     Extract type data

      stype = d(16)
      etype = d(17)
      dtype = d(18)
      hflag = d(67).eq.1.0d0

c     Set nodal temperatures: Can be specified or computed

      if(isw.gt.1) then
        tdof = d(19)
        if(tdof.le.0) then
          do i = 1,nel ! {
            th(i) = tl(i)
          end do ! i     }
        else
          do i = 1,nel ! {
            th(i) = ul(tdof,i)
          end do ! i     }
        endif
      endif

c     Output element type

      if(isw.eq.0 .and. ior.lt.0) then
        write(*,*) '   Elmt  1: 2-d Solid Plane/Axisy Linear Elastic',
     &             '  (Displ. Model)'

c     Input material properties

      elseif(isw.eq.1) then

        write(iow,2001)
        if(ior.lt.0) write(*,2001)
        call inmate(d,tdof,0,1)
        if(etype.eq.2) then
          nh1 = nh1 + 2
        endif

c       Set tdof to zero if 1, 2, or larger than ndf

        if(tdof.gt.ndf) then
          write(iow,3003)
          if(ior.lt.0) write(*,3003)
          tdof = 0
        elseif(tdof.eq.1 .or. tdof.eq.2) then
          write(iow,3004)
          if(ior.lt.0) write(*,3004)
          tdof = 0
        endif

c       Deactivate dof in element for dof > 2

        if(nint(d(16)).eq.8) then
          do i = 4,ndf
            ix(i) = 0
          end do ! i
        else
          do i = 3,ndf
            ix(i) = 0
          end do ! i
        endif

c       If temperature dof is specified activate dof

        if(tdof.gt.0) then
          ix(tdof) = 1
        endif

c       Set plot sequence

        pstyp = 2

c       Set number of projected stresses

        istv = max(istv,9)

c       Set to check in each element

        mech  = .true.
        ther  = .false.
        stype = d(16)
        etype = d(17)
        dtype = d(18)

c       Check dimensions and material for errors

        errck = .false.

        if(ndf.lt.2 .or. nen.lt.3) then
          write(  *,4000) ndf,nen
          write(ilg,4000) ndf,nen
          write(iow,4000) ndf,nen
          errck = .true.
        endif
        if(stype.eq.1 .and. etype.eq.2) then
          write(  *,4001)
          write(ilg,4001)
          write(iow,4001)
          errck = .true.
        endif

        if(dtype.le.0 .and. etype.eq.3) then
          write(  *,4002)
          write(ilg,4002)
          write(iow,4002)
          errck = .true.
        endif

        if(dtype.le.0 .and. etype.eq.5) then
          write(  *,4003)
          write(ilg,4003)
          write(iow,4003)
          errck = .true.
        endif

        if(errck) then
          call plstop()
        endif

c     Check element for jacobian errors in input data

      elseif(isw.eq.2) then
        if(nel.eq.3. .or. nel.eq.6 .or. nel.eq.7) then
          call cktris(ix,xl,shp,ndm)
        else
          call ckisop(ix,xl,shp,ndm)
        endif
        return

c     Compute mass or geometric stiffness matrix

      elseif(isw.eq.5) then

        mech =  npart.eq.ndfp(1)
        if(tdof.gt.ndm .and. hflag) then
          ther = npart.eq.ndfp(tdof)
        else
          ther = .false.
        endif

        if(imtyp.eq.1) then
          if(mech) then
            call mass2d(d,xl,ul,ix,s,p,ndf,ndm,nst)
          endif
          if(ther) then
            call therm2d(d,ul(tdof,1),xl,ix,s(tdof,tdof),p(tdof),
     &                   ndf,ndm,nst,isw)
          endif
        else
c         put call to geometric stiffness routine here
        endif
        return

c     Compute surface tractions

      elseif(isw.eq.7) then
        call surf2d(d,ul,xl,ma,ndf,ndm,nel,mct,nst, p,s)
        return

c     Compute damping matrix

      elseif(isw.eq.9) then

        mech =  npart.eq.ndfp(1)
        if(tdof.gt.ndm .and. hflag) then
          ther = npart.eq.ndfp(tdof)
        else
          ther = .false.
        endif

        if(mech) then
          call damp2d(d,xl,ix,s,ndf,ndm,nst)
        endif
        if(ther) then
          call therm2d(d,ul(tdof,1),xl,ix,s(tdof,tdof),p(tdof),
     &                  ndf,ndm,nst,isw)
        endif
        return

c     Body Force computation

      elseif(isw.eq.15 .or. isw.eq.23) then
        call sbody2d(d,xl,ix, p,ndm,ndf ,isw)
        return

c     Normal to surface computation for projections

      elseif(isw.eq.24) then
        call pnorm2d(p,xl,ndm,nel)
        return

c     Face set data

      elseif(isw.eq.26) then
        call pcorner2d()
        return

c     Compute residuals and tangents for parts

      else
        mech =  npart.eq.ndfp(1) .or. isw.eq.4 .or. isw.eq.8
        if(tdof.gt.ndm .and. hflag) then
          ther = npart.eq.ndfp(tdof) .or. isw.eq.4 .or. isw.eq.8
        else
          ther = .false.
        endif

      endif

c     Check if symmetric or unsymmetric tangent

      if(neqs.eq.1) then
        mther = ther
      else
        mther = .false.
      endif

c     Compute stress-divergence vector (p) and stiffness matrix (s)

      if(mech) then

c       Displacement Model

        if(etype.eq.1) then

          if(dtype.gt.0) then
            call sld2d1(d,ul,xl,ix,th,s,p,ndf,ndm,nst,isw,mther)
          else
            call fld2d1(d,ul,xl,ix,s,p,ndf,ndm,nst,isw)
          endif

c       Mixed Model (B-Bar)

        elseif(etype.eq.2) then

          if(dtype.gt.0) then
            call sld2d2(d,ul,xl,ix,th,s,p,ndf,ndm,nst,isw,mther)
          else
            call fld2d2(d,ul,xl,ix,s,p,ndf,ndm,nst,isw)
          endif

c       Enhanced Strain Model (B-Bar)

        elseif(etype.eq.3) then

          if(dtype.gt.0) then
            call sld2d3(d,ul,xl,ix,th,s,p,ndf,ndm,nst,isw)
          else
c           call fld2d3(d,ul,xl,ix,s,p,ndf,ndm,nst,isw)
          endif

c       Energy Conserving Model

        elseif(etype.eq.4) then

          if(dtype.gt.0) then
            call sld2d1(d,ul,xl,ix,th,s,p,ndf,ndm,nst,isw,mther)
          else
            call fld2d4(d,ul,xl,ix,s,p,ndf,ndm,nst,isw)
          endif

c       Mixed-Enhanced Strain Model

        elseif(etype.eq.5) then

          if(dtype.gt.0) then
            call sld2d5(d,ul,xl,ix,th,s,p,ndf,ndm,nst,isw)
          else
c           call fld2d5(d,ul,xl,ix,s,p,ndf,ndm,nst,isw)
          endif

        endif
      endif

c     Compute thermal vector (p) and matrix (s)

      if(ther) then

c       Prevent accumulation of weight for projections

        if(isw.eq.8) then
          pdof = 1
          do i = 1,nen
            p(i) = 0.0d0
          end do ! i
        else
          pdof = tdof
        endif

        call therm2d(d,ul(tdof,1),xl,ix,s(pdof,pdof),p(pdof),
     &               ndf,ndm,nst,isw)
      endif

c     Formats for input-output

2001  format(
     & /5x,'T w o   D i m e n s i o n a l   S o l i d   E l e m e n t'/)

3003  format(' *WARNING* Thermal d.o.f. > active d.o.f.s : Set to 0')

3004  format(' *WARNING* Thermal d.o.f. can not be 1 or 2: Set to 0')

4000  format(' *ERROR* Problem control record incorrect:'/
     &       '         DOFs/Node (ndf) = ',i4,' - Should be 2 or more'/
     &       '         Nodes/Elm (nen) = ',i4,' - Should be 3 or more')

4001  format(' *ERROR* Plane stress not implemented for mixed type.')

4002  format(' *ERROR* Enhanced element type does not currently exist.')

4003  format(' *ERROR* Mixed-Enhanced element type does not currently',
     &       ' exist.')

      end
