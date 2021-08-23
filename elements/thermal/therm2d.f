c$Id: therm2d.f,v 1.1 2006/11/20 20:32:56 rlt Exp $
      subroutine therm2d(d,ul,xl,ix,s,p,ndf,ndm,nst,isw)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c     Two dimensional (plane/axisymmetric) Linear Thermal Element

c     N.B.  Surface flux loading may be specified using: ELMT08

c-----[--.----+----.----+----.-----------------------------------------]

c        This is a two dimensional element which can analyze plane
c        or axisymmetric geometries.  Set control parameters as
c        follows:

c           ndm - set to 2     (x,y or r,z-coords)
c           ndf - set > or = 1 (nodal temperatures)
c           nel - set > or = 4

c                    A eta
c             4      |      3
c              o-----|-----o
c              |     |     |
c              |     |     |
c              |     +-----|----> xi
c              |           |
c              |           |
c              o-----------o
c             1             2

c               Node numbering
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'bdata.h'
      include  'cdata.h'
      include  'complx.h'
      include  'eldata.h'
      include  'elplot.h'
      include  'eltran.h'
      include  'fdata.h'
      include  'iofile.h'
      include  'mdata.h'
      include  'part0.h'
      include  'pmod2d.h'
      include  'rdata.h'
      include  'comblk.h'

      logical   quad
      integer   ndf,ndm,nst,isw, i,j, i1,j1, l,lint, tdof, ix(*)
      real*8    xx,yy, xsj, a1,a2,a3,a4,shj,tdot,cfac,lfac, hh,tinf
      real*8    d(*),ul(ndf,nen,*),xl(ndm,*),s(nst,*),p(*)
      real*8    temp,gradt(2),flux(2),dd(2,2)
      real*8    shp(3,25),sg(3,25),el(4,7)

      save

c     Set mass factors

      if(d(7).ge.0.0d0) then
        cfac = d(7)
        lfac = 1.d0 - cfac
      else
        cfac = 0.0d0
        lfac = 0.0d0
      endif

c     Input material properties

      if(isw.eq.1) then

        if(ior.lt.0) write(*,2000)
        write(iow,2000)
        call inmate(d,tdof,0,6)

c       Delete unused parameters

        do i = 2,ndf
          ix(i) = 0
        end do ! i

c       Set to preclude sloping boundary transformations

        ea(1,-iel) = 0
        ea(2,-iel) = 0

c       Set plot sequence

        pstyp = 2

c     Check of mesh if desired (chec)

      elseif(isw.eq.2) then

        if(nel.eq.3 .or. nel.eq.6 .or. nel.eq.7) then
          call cktris(ix,xl,shp,ndm)
        else
          call ckisop(ix,xl,shp,ndm)
        endif

c     Compute conductivity (stiffness) matrix

      elseif(isw.eq.3 .or. isw.eq.6) then

        if(nel.eq.6) then
          if(nint(d(182)).gt.0) then
            call tint2dn(nel,lint,el)
          else
            l =  7
            call tint2d(l,lint,el)
          endif
          quad = .false.
        elseif(nel.eq.7) then
          if(nint(d(182)).gt.0) then
            call tint2dn(nel,lint,el)
          else
            l =  7
            call tint2d(l,lint,el)
          endif
          quad = .false.
        else
          quad = .true.
          if(nint(d(182)).gt.0) then
            call int2dn(nel,lint,sg)
          else
            l = nint(d(5))
            call int2d (l,lint,sg)
          endif
        endif

c       Get global dof for thermal variable

        tdof = max(1,nint(d(19)))
        hh   = d(127)
        tinf = d(128)

        do l = 1,lint

          if(quad) then
            call shp2d(sg(1,l),xl,shp,xsj,ndm,nel,ix,.false.)
            xsj = xsj*sg(3,l)*d(14)
          else
            call trishp(el(1,l),xl,ndm,nel-4,xsj,shp)
            xsj = xsj*el(4,l)*d(14)
          endif

c         Compute flux

          call thfx2d(d,xl,ul,xx,yy,shp, temp,gradt,flux,dd,ndm,ndf,nel)

c         Save data for tplot

          j       = 4*(l-1)
          tt(j+1) = flux(1)
          tt(j+2) = flux(2)
          tt(j+3) = gradt(1)
          tt(j+4) = gradt(2)

c         Compute thermal rate

          tdot = 0.0d0
          do j = 1,nel
            tdot = tdot + shp(3,j)*ul(1,j,4)
          end do ! j

          if(stype.eq.3) then
            xsj = xsj*xx
          endif

          j1 = 1
          do j = 1,nel

            a1 = (dd(1,1)*shp(1,j) + dd(1,2)*shp(2,j))*xsj
            a2 = (dd(2,1)*shp(1,j) + dd(2,2)*shp(2,j))*xsj
            a3 = d(4)*d(64)*shp(3,j)*xsj
            a4 = d(127)*shp(3,j)*xsj

c           Compute residual

            p(j1) = p(j1) - a1*gradt(1) - a2*gradt(2)
     &                    - a3*(cfac*tdot + lfac*ul(1,j,4))
     &                    - a4*(temp - d(128))
     &                    + d(66)*shp(3,j)*xsj*dm

c           Compute tangent

            a1 = a1*ctan(1)
            a2 = a2*ctan(1)
            if(shflg) then
              if(cplxfl) then
                a3 = 0.0d0
              else
                a3 = a3*ctan(3)
              endif
            elseif(ndfo(tdof).gt.0) then
              a3 = a3*ctan(2)
            else
              a3 = 0.0d0
            endif
            a4 = a4*ctan(1) + a3*cfac

c           Lumped rate terms

            s(j1,j1) = s(j1,j1) + a3*lfac

c           Consistent rate and conductivity terms

            i1 = 1
            do i = 1,nel
              s(i1,j1) = s(i1,j1) + a1*shp(1,i) + a2*shp(2,i)
     &                            + a4*shp(3,i)
              i1 = i1 + ndf
            end do ! i
            j1 = j1 + ndf
          end do ! j

c         Complex solution part

          if(cplxfl) then
            call ctherm2d(d,xl,ul, shp,xsj,cfac,lfac,
     &                    s(1,nst+1),p(nst+1))
          endif
        end do ! l

c     Output heat flux

      elseif(isw.eq.4) then

        if(nel.eq.6) then
          if(nint(d(182)).gt.0) then
            l =  6
          else
            l =  7
          endif
          quad = .false.
          call tint2d(l,lint,el)
        elseif(nel.eq.7) then
          if(nint(d(182)).gt.0) then
            l = -7
          else
            l =  7
          endif
          quad = .false.
          call tint2d(l,lint,el)
        else
          quad = .true.
          if(nint(d(182)).gt.0) then
            call int2dn(nel,lint,sg)
          else
            l = nint(d(5))
            call int2d (l,lint,sg)
          endif
        endif

        do l=1,lint

        if(quad) then
          call shp2d(sg(1,l),xl,shp,xsj,ndm,nel,ix,.false.)
        else
          call trishp(el(1,l),xl,ndm,nel-4,xsj,shp)
        endif

c       Compute flux and gradients

        call thfx2d(d,xl,ul, xx,yy,shp,temp,gradt,flux,dd, ndm,ndf,nel)

        a4 = -d(127)*(temp - d(128))

        mct = mct - 1
        if(mct.le.0) then
          write(iow,2002) o,head
          if(d(127).gt.0.0d0) then
            write(iow,2004)
          endif
          if(ior.lt.0 .and. pfr) then
            write(*,2002) o,head
            if(d(127).gt.0.0d0) then
              write(*,2004)
            endif
          endif
          mct = 50
        endif
        write(iow,2003) n,ma,xx,yy,flux,gradt
        if(d(127).gt.0.0d0) then
          write(iow,2005) a4
        endif
        if(ior.lt.0 .and. pfr) then
          write(*,2003) n,ma,xx,yy,flux,gradt
          if(d(127).gt.0.0d0) then
            write(*,2005) a4
          endif
        endif

        end do ! l

c     Compute heat capacity (mass) matrix

      elseif(isw.eq.5) then
        if(nel.eq.6) then
          if(nint(d(182)).gt.0) then
            l =  6
          else
            l =  7
          endif
          quad = .false.
          call tint2d(l,lint,el)
        elseif(nel.eq.7) then
          if(nint(d(182)).gt.0) then
            l = -7
          else
            l =  7
          endif
          quad = .false.
          call tint2d(l,lint,el)
        else
          if(nint(d(182)).gt.0) then
            call int2dn(nel,lint,sg)
          else
            l = nint(d(5))
            call int2d (l,lint,sg)
          endif
          quad = .true.
        endif
        do l = 1,lint
          if(quad) then
            call shp2d(sg(1,l),xl,shp,xsj,ndm,nel,ix,.false.)
            xsj = xsj*sg(3,l)*d(14)
          else
            call trishp(el(1,l),xl,ndm,nel-4,xsj,shp)
            xsj = xsj*el(4,l)*d(14)
          endif
          if(stype.eq.3) then
            xx = 0.0d0
              do i = 1,nel
              xx = xx + shp(3,i)*xl(1,i)
            end do ! i
            xsj = xsj*xx
          endif
          j1 = 1
          do j = 1,nel
            shj = d(4)*d(64)*shp(3,j)*xsj

c           Lumped capacity (lmas)

            p(j1) = p(j1) + shj
            i1 = 1

c           Consistent (interpolated ) capacity (mass)

            s(i1,i1) = s(i1,j1) + shj*lfac
            do i = 1,nel
              s(i1,j1) = s(i1,j1) + shj*shp(3,i)*cfac
              i1 = i1 + ndf
            end do ! i
            j1 = j1 + ndf
          end do ! j
        end do ! l

c     Compute surface flux loading (not implemented)

c     elseif(isw.eq.7) then

c     Compute nodal heat flux for output/plots

      elseif(isw.eq.8) then

        call thcn2d(ix,d,xl,ul,shp,p,s,p(nen+1),ndf,ndm,nel)

c     Compute error data for heat flux

      elseif(isw.eq.11) then

        call ther2d(ix,d,xl,ul,shp,s,ndf,ndm,nel,nen)

c     External node determination

      elseif(isw.eq.26) then

        call pcorner2d()

      endif

c     Formats

2000  format(5x,'F o u r i e r   H e a t   C o n d u c t i o n')

2002  format(a1,20a4//5x,'Element Flux'//'  Elmt Mat 1-Coord  2-Coord'
     &            ,'      1-Flux      2-Flux      1-Grad      2-Grad')
2003  format(i6,i4,0p,2f9.3,1p,4e12.3)

2004  format(28x,' Surf. Conv.')
2005  format(28x,1p,1e12.3)

      end

      subroutine thcn2d(ix,d,xl,ul,shp,dt,st,ser,ndf,ndm,nel)

      implicit  none

      include  'iodata.h'
      include  'cdata.h'
      include  'prstrs.h'
      include  'strnum.h'

      logical   quad
      integer   ndf,ndm,nel, j,l,lint, ix(*), ixl(9)
      real*8    xx,yy,xsj,xg,d(*), sg(3,16),el(4,7)
      real*8    dt(*),st(nen,*),ser(*),xl(ndm,*),shp(3,*)
      real*8    temp,gradt(2),flux(2),dd(2,2),ul(ndf,*)

      save

      data      ixl/ 1,2,3,4,5,6,7,8,9 /

c     Lumped projection routine

      if(nel.eq.6) then
        l    =  6
        quad = .false.
        call tint2d(l,lint,el)
      elseif(nel.eq.7) then
        l    = -7
        quad = .false.
        call tint2d(l,lint,el)
      else
        l    =  d(5)
        l    =  max(2,l)
        quad = .true.
        call int2d(l,lint,sg)
        if(nel.eq.8) then
          call meanx(ix,xl,ndm)
        endif
      endif

      do l = 1,lint
        if(quad) then
          if(nel.eq.8) then
            call shp2d(sg(1,l),xl,shp,xsj,ndm, 9,ixl,.false.)
          else
            call shp2d(sg(1,l),xl,shp,xsj,ndm,nel,ix,.false.)
          endif
          xsj = xsj*sg(3,l)*d(14)
        else
          call trishp(el(1,l),xl,ndm,nel-4,xsj,shp)
          xsj = xsj*el(4,l)*d(14)
        endif

        call thfx2d(d,xl,ul, xx,yy,shp,temp,gradt,flux,dd, ndm,ndf,nel)

        temp = -d(127)*(temp - d(128))

c       Compute lumped projection and assemble stress integrals

        do j = 1,nel
          xg      = xsj*shp(3,j)
          dt(j)   = dt(j) + xg
          st(j,7) = st(j,7) + flux(1)*xg
          st(j,8) = st(j,8) + flux(2)*xg
          st(j,9) = st(j,9) + temp   *xg
          ser(j)  = ser(j)  + erav   *xg
        end do ! j
      end do ! l

      iste = 9

      end

      subroutine thfx2d(d,xl,ul, xx,yy,shp, temp,gradt,flux,dd,
     &                  ndm,ndf,nel)

c     Compute thermal gradient and flux

      implicit  none

      integer   ndm,ndf,nel, i
      real*8    d(*),xl(ndm,*),ul(ndf,*), shp(3,*)
      real*8    xx,yy,psi,cs,sn,c2,s2, temp,gradt(2),flux(2),dd(2,2)

      save

      temp     = 0.0d0
      gradt(1) = 0.0d0
      gradt(2) = 0.0d0
      xx       = 0.0d0
      yy       = 0.0d0
      do i = 1,nel
        gradt(1) = gradt(1) + shp(1,i)*ul(1,i)
        gradt(2) = gradt(2) + shp(2,i)*ul(1,i)
        temp     = temp     + shp(3,i)*ul(1,i)
        xx       = xx       + shp(3,i)*xl(1,i)
        yy       = yy       + shp(3,i)*xl(2,i)
      end do ! i

c     Compute thermal flux

      psi = d(31)
      cs  = cos(psi)
      sn  = sin(psi)
      c2  = cs*cs
      s2  = sn*sn
      cs  = cs*sn

      dd(1,1) = c2*d(61) + s2*d(62)
      dd(2,2) = s2*d(61) + c2*d(62)
      dd(1,2) = cs*(d(61) - d(62))
      dd(2,1) = dd(1,2)

      flux(1) = -dd(1,1)*gradt(1) - dd(1,2)*gradt(2)
      flux(2) = -dd(2,1)*gradt(1) - dd(2,2)*gradt(2)

      end

      subroutine ther2d(ix,d,xl,ul,shp,st,ndf,ndm,nel,nen)

      implicit  none

      include  'adapt1.h'
      include  'adapt2.h'
      include  'errind.h'

      integer   ndf,ndm,nel,nen, i,j,ii, ix(*)
      real*8    g,xx,yy,xsj,detd, st(nen,*),xl(ndm,*),shp(3,*)
      real*8    d(*),gradt(2),flux(2),dd(2,2),ul(ndf,*),ss(9),tt(9)
      real*8    temp,gradp(2),fluxp(2),ddp(2,2),sg(2)

      save

      data      ss/-1.d0, 1.d0,1.d0,-1.d0, 0.d0,1.d0,0.d0,-1.d0,0.d0/
      data      tt/-1.d0,-1.d0,1.d0, 1.d0,-1.d0,0.d0,1.d0, 0.d0,0.d0/

c     Simple routine

      vfem   = 0.d0
      vproj  = 0.d0
      verror = 0.d0
      vener  = 0.d0
      venere = 0.d0
      heta   = 0.0d0
      g      = 1.d0/sqrt(3.0d0)
      do ii = 1,4
        sg(1) = ss(ii)*g
        sg(2) = tt(ii)*g
        call shp2d(sg,xl,shp,xsj,ndm,nel,ix,.false.)
        call thfx2d(d,xl,ul, xx,yy,shp,temp,gradt,flux,dd, ndm,ndf,nel)
        do i = 1,2
          fluxp(i) = 0.0d0
        end do ! i
        do i = 1,nel
          do j = 1,2
            fluxp(j) = fluxp(j) + shp(3,i)*st(i,j+6)
          end do ! j
        end do ! i

c       Compute integral of stress squares for error indicator use

        detd     =  1.d0/(dd(1,1)*dd(2,2) - dd(1,2)*dd(2,1))
        ddp(1,1) =  dd(2,2)*detd
        ddp(1,2) = -dd(1,2)*detd
        ddp(2,1) = -dd(2,1)*detd
        ddp(2,2) =  dd(1,1)*detd
        gradp(1) = -ddp(1,1)*fluxp(1) - ddp(1,2)*fluxp(2)
        gradp(2) = -ddp(2,1)*fluxp(1) - ddp(2,2)*fluxp(2)

        heta = heta + xsj
        do i = 1,2
         vfem    = vfem   + flux(i)*flux(i)*xsj
         vproj   = vproj  + fluxp(i)*fluxp(i)*xsj
         verror  = verror + ((fluxp(i)-flux(i))**2)*xsj
         vener   = vener  + flux(i)*gradt(i)*xsj
         venere  = venere + (fluxp(i)-flux(i))*(gradp(i)-gradt(i))*xsj
        end do ! i

      end do ! ii
      arsq  =  arsq  + heta
      efem  =  efem  + vfem
      eproj =  eproj + vproj
      eerror=  eerror+ verror
      eener =  eener + vener
      eenere=  eenere+ venere
      areai = heta

c     Check for triangle

      if(nel.lt.4 .or. ix(1).eq.ix(2) .or. ix(2).eq.ix(3)
     1            .or. ix(3).eq.ix(4) .or. ix(4).eq.ix(1) ) then
        heta = heta*2.d0
      endif
      heta  =  d(50)*sqrt(heta)

      end

      subroutine ctherm2d(d,xl,ul, shp,xsj,cfac,lfac, s,p)

      implicit   none

      include   'cdata.h'
      include   'eldata.h'
      include   'pmod2d.h'
      include   'rdata.h'
      include   'sdata.h'

      real*8     d(*),xl(ndm,*),ul(ndf,nen,*), shp(3,*)
      real*8     s(nst,*),p(*), xsj, cfac,lfac

      integer    i,j, i1,j1
      real*8     xx,yy,tdot,temp, a1,a2,a3,a4, omega
      real*8     gradt(2),flux(2),dd(2,2)

c     Shift frequency

      omega = sqrt(abs(shift))

c     Compute flux

      call thfx2d(d,xl,ul(1,1,8),xx,yy,shp, temp,gradt,flux,dd,
     &            ndm,ndf,nel)

c     Compute thermal rate

      tdot = 0.0d0
      do j = 1,nel
        tdot = tdot + shp(3,j)*ul(1,j,11)
      end do ! j

      if(stype.eq.3) then
        xsj = xsj*xx
      endif

      j1 = 1
      do j = 1,nel

        a1 = (dd(1,1)*shp(1,j) + dd(1,2)*shp(2,j))*xsj
        a2 = (dd(2,1)*shp(1,j) + dd(2,2)*shp(2,j))*xsj
        a3 = d(4)*d(64)*shp(3,j)*xsj*omega
        a4 = d(127)*shp(3,j)*xsj

c       Compute residual

        p(j1) = p(j1) - a1*gradt(1) - a2*gradt(2)
     &                - a3*(cfac*tdot + lfac*ul(1,j,4))
     &                - a4*(temp - d(128))
     &                + d(66)*shp(3,j)*xsj*dm

c       Compute tangent

        a4 = a3*cfac

c       Lumped rate terms

        s(j1,j1) = s(j1,j1) + a3*lfac

c       Consistent rate and conductivity terms

        i1 = 1
        do i = 1,nel
          s(i1,j1) = s(i1,j1) + a4*shp(3,i)
          i1 = i1 + ndf
        end do ! i
        j1 = j1 + ndf
      end do ! j

      end
