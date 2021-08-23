c$Id: solid3d.f,v 1.1 2006/11/20 20:32:56 rlt Exp $
      subroutine solid3d(d,ul,xl,ix,tl,s,p,ndf,ndm,nst,isw)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Three Dimensional Solid Element Driver

c      Inputs:
c         d(*)  - Element parameters
c         ul(ndf,*) - Current nodal solution parameters
c         xl(ndm,*) - Nodal coordinates
c         ix(*)     - Global nodal connections
c         ndf       - Degree of freedoms/node
c         ndm       - Mesh coordinate dimension
c         nst       - Element array dimension
c         isw       - Solution option switch

c      Outputs:
c         s(nst,*)  - Element array
c         p(ndf,*)  - Element vector
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
      real*8    d(*),ul(ndf,*),xl(ndm,*),tl(*),s(nst,nst),p(nst)
      real*8    shp(4,64),th(64)

      save

c     Extract type data

      stype = d(16)
      etype = d(17)
      dtype = d(18)
      hflag = d(67).eq.1.d0

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

c     Set zero state

c     Output element type

      if(isw.eq.0 .and. ior.lt.0) then
        write(*,*) '   Elmt  1: 3-d Solid Linear/Finite Defm. Element.'

c     Input material properties

      elseif(isw.eq.1) then

        write(iow,2000)
        if(ior.lt.0) write(*,2000)
c       call inmate(d,tdof,ndf*nen,1)
        call inmate(d,tdof,   0   ,1)
        nh1 = nh1 + 2
        if(etype.eq.3) then
          nh1 = nh1 + 12
        endif

c       Set tdof to zero if 1, 2, 3, or larger than ndf

        if(tdof.gt.ndf) then
          write(iow,3000)
          if(ior.lt.0) write(*,3000)
          tdof = 0
        elseif(tdof.ge.1 .and. tdof.le.3) then
          write(iow,3001)
          if(ior.lt.0) write(*,3001)
          tdof = 0
        endif

c       Deactivate dof in element for dof > 3

        do i = 4,ndf
          ix(i) = 0
        end do ! i

c       If temperature dof is specified activate dof

        if(tdof.gt.0) then
          ix(tdof) = 1
        endif

c       Set plot sequence

        pstyp = 3

c       Check for errors in problem size

        errck = .false.
        if(ndf.lt.3 .or. nen.lt.4) then
          write(  *,4000) ndf,nen
          write(ilg,4000) ndf,nen
          write(iow,4000) ndf,nen
          errck = .true.
        endif

        if(etype.eq.5) then
          write(  *,4001)
          write(ilg,4001)
          write(iow,4001)
          errck = .true.
        endif

        if(errck) then
          call plstop()
        endif

c     Check element for errors

      elseif(isw.eq.2) then

        if(nel.eq.4 .or. nel.eq.10) then
          call cktets ( n, ix, xl, ndm, nel, shp )
        else
          call ckbrk8 ( n, ix, xl, ndm, nel, shp )
        endif

c     Compute mass or geometric stiffness matrix

      elseif(isw.eq.5) then
        mech =  npart.eq.ndfp(1)
        if(tdof.gt.ndm .and. hflag) then
          ther = npart.eq.ndfp(tdof)
        else
          ther = .false.
        endif

c       Mass contribution

        if(imtyp.eq.1) then
          if(mech) then
            call mass3d(d,xl,s,p,ndf,ndm,nst)
          elseif(ther) then
            call therm3d(d,ul(tdof,1),xl,ix,s(tdof,tdof),p(tdof),
     &                   ndf,ndm,nst,isw)
          endif

c       Geometric stiffness here

        else
        endif

c     Compute damping matrix

      elseif(isw.eq.9) then
        mech =  npart.eq.ndfp(1)
        if(tdof.gt.ndm .and. hflag) then
          ther = npart.eq.ndfp(tdof)
        else
          ther = .false.
        endif

c       Damping contribution

        if(mech) then
          call damp3d(d,xl,s,ndf,ndm,nst)
        elseif(ther) then
          call therm3d(d,ul(tdof,1),xl,ix,s(tdof,tdof),p(tdof),
     &                 ndf,ndm,nst,isw)
        endif

c     Body Force computation

      elseif(isw.eq.15 .or.isw.eq.23) then

        call sbody3d(d,xl, p,ndm,ndf, isw)

c     Compute face normals

      elseif(isw.eq.24) then

        call pnorm3d(p,xl,ndm,nel)

c     Compute boundary nodes

      elseif(isw.eq.26) then

        call pcorner3d()

c     Compute residuals and tangents for parts

      elseif(isw.eq.3 .or.  isw.eq. 4  .or. isw.eq. 6 .or.
     &      (isw.ge.8 .and. isw.le.14) .or. isw.eq.16 ) then

        mech =  npart.eq.ndfp(1) .or. isw.eq.4 .or.
     &          isw.eq.8 .or. isw.eq.14
        if(tdof.ne.0 .and. tdof.gt.ndm .and. hflag) then
          ther = npart.eq.ndfp(tdof) .or. isw.eq.4 .or.
     &          isw.eq.8 .or. isw.eq.14
        else
          ther = .false.
        endif

c       Check if symmetric or unsymmetric tangent

        if(neqs.eq.1) then
          mther = ther
        else
          mther = .false.
        endif

c       Compute stress-divergence vector (p) and stiffness matrix (s)

        if(mech) then

c         Displacement Model

          if(etype.eq.1) then

            if(dtype.gt.0) then
              call sld3d1(d,ul,xl,th,s,p,ndf,ndm,nst,isw)
            else
              call fld3d1(d,ul,xl,s,p,ndf,ndm,nst,isw)
            endif

c         Mixed Model (B-Bar)

          elseif(etype.eq.2) then

            if(dtype.gt.0) then
              call sld3d2(d,ul,xl,th,s,p,ndf,ndm,nst,isw, mther)
            else
              call fld3d2(d,ul,xl,s,p,ndf,ndm,nst,isw)
            endif

c         Enhanced Strain Model

          elseif(etype.eq.3) then

            if(dtype.gt.0) then
              call sld3d3(d,ul,xl,th,s,p,ndf,ndm,nst,isw)
            else
              call fld3d3(d,ul,xl,s,p,ndf,ndm,nst,isw, .true.)
            endif

c         Energy Conserving Model

          elseif(etype.eq.4) then

            if(dtype.gt.0) then
              call sld3d1(d,ul,xl,th,s,p,ndf,ndm,nst,isw)
            else
              call fld3d4(d,ul,xl,s,p,ndf,ndm,nst,isw)
            endif

c         Mixed-Enhanced Model

          elseif(etype.eq.5) then

          endif

        endif

c       Compute thermal vector (p) and matrix (s)

        if(ther) then

c         Prevent multiple accumulation of projection weights

          if(isw.eq.8) then
            pdof = 1
            do i = 1,8
              p(i) = 0.0d0
            end do ! i
          else
            pdof = tdof
          endif
          call therm3d(d,ul(tdof,1),xl,ix,s(pdof,pdof),p(pdof),
     &                 ndf,ndm,nst,isw)
        endif

      endif

c     Formats for input-output

2000  format(
     & /5x,'T h r e e   D i m e n s i o n a l   S o l i d',
     &     '   E l e m e n t'/)

3000  format(' *WARNING* Thermal d.o.f. > active d.o.f.s : Set to 0')

3001  format(' *WARNING* Thermal d.o.f. can not be 1 to 3: Set to 0')

4000  format(' *ERROR* Problem control record incorrect:'/
     &       '         DOFs/Node (ndf) = ',i4,' - Should be 3 or more'/
     &       '         Nodes/Elm (nen) = ',i4,' - Should be 4 or more')

4001  format(' *ERROR* Mixed-Enhanced element type does not currently',
     &       ' exist.')

      end
