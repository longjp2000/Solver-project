c$Id: therm3d.f,v 1.1 2006/11/20 20:32:56 rlt Exp $
      subroutine therm3d(d,ul,xl,ix,s,p,ndf,ndm,nst,isw)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c     Three dimensional Linear Thermal Element

c-----[--.----+----.----+----.-----------------------------------------]

c        This is a three dimensional element which can analyze
c        general geometries.  Set control parameters as
c        follows:

c           ndm - set to 3      (x,y or r,z-coords)
c           ndf - set > or =  1 (nodal temperatures)
c           nel - set > or =  8 for  8-node linear brick
c                     > or = 27 for 27-node quadratic brick
c                     > or = 64 for 64-node cubic brick
c                     > or =  4 for  4-node linear tetrahedron
c                     > or = 10 for 10-node qudratic tetrahedron

c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'bdata.h'
      include  'cdata.h'
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

      integer   ndf,ndm,nst,isw, i,j, i1,j1, l,lint, tdof, ord, ix(*)
      real*8    xsj, a0,a1,a2,a3,shj,tdot,lfac,cfac, sv(5,16)
      real*8    d(*),ul(ndf,nen,*),xl(ndm,*),s(nst,*),p(*)
      real*8    xx(3),gradt(3),flux(3),dd(3,3),shp(4,64),sg(4,64),sg0(3)

      save

      data      sg0 /0.0d0, 0.0d0, 0.0d0/

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

        pstyp = 3

c     Check of mesh if desired (chec)

      elseif(isw.eq.2) then

        call ckbrk8(n,ix,xl,ndm,nel,shp)

c     Compute conductivity (stiffness) matrix

      elseif(isw.eq.3 .or. isw.eq.6) then

        if(nel.eq.4) then
          ord = 1
          if(nint(d(182)).gt.0) then
            call tint3dn(nel,lint,sv)
          else
            l =  2
            call tint3d (l,lint,sv)
          endif
        elseif(nel.eq.10) then
          ord =  2
          l   =  3
          call tint3d(l,lint,sv)
        elseif(nel.eq.11) then
          ord =  2
          if(nint(d(182)).gt.0) then
            call tint3dn(nel,lint,sv)
          else
            l =  4
            call tint3d (l,lint,sv)
          endif
        else
          ord = 0
          if(nint(d(182)).gt.0) then
            call int3dn(nel,lint,sg)
          else
            l = nint(d(5))
            call int3d(l,lint,sg)
          endif
        endif

c       Check for thermal dof

        tdof = max(1,nint(d(19)))

        do l = 1,lint

          if(ord.gt.0) then
            call tetshp(sv(1,l),xl,ndm,ord,xsj,shp)
            xsj = xsj*sv(5,l)
          else
            call shp3d(sg(1,l),xsj,shp,xl,ndm,nel)
            xsj = xsj*sg(4,l)
          endif

c         Compute flux

          call thfx3d(d,xl,ul, xx,shp,gradt,flux,dd, ndm,ndf,nel)

c         Save data for tplot

          j       = 6*(l-1)
          tt(j+1) = flux(1)
          tt(j+2) = flux(2)
          tt(j+3) = flux(3)
          tt(j+4) = gradt(1)
          tt(j+5) = gradt(2)
          tt(j+6) = gradt(3)

c         Compute thermal rate

          tdot = 0.0d0
          do j = 1,nel
            tdot = tdot + shp(4,j)*ul(1,j,4)
          end do ! j

          j1 = 1
          do j = 1,nel

            a1 = (dd(1,1)*shp(1,j) + dd(1,2)*shp(2,j)
     &         +  dd(1,3)*shp(3,j))*xsj
            a2 = (dd(2,1)*shp(1,j) + dd(2,2)*shp(2,j)
     &         +  dd(2,3)*shp(3,j))*xsj
            a3 = (dd(3,1)*shp(1,j) + dd(3,2)*shp(2,j)
     &         +  dd(3,3)*shp(3,j))*xsj

            a0 = d(4)*d(64)*shp(4,j)*xsj

c           Compute residual

            p(j1) = p(j1) - a1*gradt(1) - a2*gradt(2) - a3*gradt(3)
     &                    - a0*(cfac*tdot + lfac*ul(1,j,4))
     &                    + d(66)*shp(4,j)*xsj*dm

c           Compute tangent

            if(shflg) then
              a0 = a0*ctan(3)
            elseif(ndfo(tdof).gt.0) then
              a0 = a0*ctan(2)
            else
              a0 = 0.0d0
            endif
            a1 = a1*ctan(1)
            a2 = a2*ctan(1)
            a3 = a3*ctan(1)

c           Lumped rate terms

            s(j1,j1) = s(j1,j1) + a0*lfac

            i1 = 1
            do i = 1,nel
              s(i1,j1) = s(i1,j1) + a1*shp(1,i) + a2*shp(2,i)
     &                            + a3*shp(3,i) + a0*shp(4,i)*cfac
              i1 = i1 + ndf
            end do ! i
            j1 = j1 + ndf
          end do ! j
        end do ! l

c     Output heat flux

      elseif(isw.eq.4) then

        if(nel.eq.4 .or. nel.eq.10) then
          sv(1,1) = 0.25d0
          sv(2,1) = 0.25d0
          sv(3,1) = 0.25d0
          sv(4,1) = 0.25d0
          if(nel.eq.4) then
            ord = 1
          else
            ord = 2
          endif
          call tetshp(sv,xl,ndm,ord,xsj,shp)
        else
          call shp3d(sg0,xsj,shp,xl,ndm,nel)
        endif

c       Compute flux and gradients

        call thfx3d(d,xl,ul, xx,shp,gradt,flux,dd, ndm,ndf,nel)

        mct = mct - 1
        if(mct.le.0) then
          write(iow,2002) o,head
          if(ior.lt.0 .and. pfr) then
            write(*,2002) o,head
          endif
          mct = 50
        endif
        write(iow,2003) n,ma,xx,flux,gradt
        if(ior.lt.0 .and. pfr) then
          write(*,2003) n,ma,xx,flux,gradt
        endif

c     Compute heat capacity (mass) matrix

      elseif(isw.eq.5) then

        if(nel.eq.4) then
          ord = 1
          if(nint(d(182)).gt.0) then
            call tint3dn(nel,lint,sv)
          else
            l =  2
            call tint3d (l,lint,sv)
          endif
        elseif(nel.eq.10) then
          ord =  2
          l   =  3
          call tint3d(l,lint,sv)
        elseif(nel.eq.11) then
          ord =  2
          if(nint(d(182)).gt.0) then
            call tint3dn(nel,lint,sv)
          else
            l =  4
            call tint3d (l,lint,sv)
          endif
        else
          ord = 0
          if(nint(d(182)).gt.0) then
            call int3dn(nel,lint,sg)
          else
            l = nint(d(5))
            call int3d(l,lint,sg)
          endif
        endif

        do l = 1,lint
          if(ord.gt.0) then
            call tetshp(sv(1,l),xl,ndm,ord,xsj,shp)
            xsj = xsj*sv(5,l)
          else
            call shp3d(sg(1,l),xsj,shp,xl,ndm,nel)
            xsj = xsj*sg(4,l)
          endif
          j1 = 1
          do j = 1,nel
            shj = d(4)*d(64)*shp(4,j)*xsj

c           Lumped capacity (lmas)

            p(j1) = p(j1) + shj
            i1 = 1

c           Consistent capacity (cmas)

            do i = 1,nel
              s(i1,j1) = s(i1,j1) + shj*shp(4,i)
              i1 = i1 + ndf
            end do ! i
            j1 = j1 + ndf
          end do ! j
        end do ! l

c     Compute surface flux loading (not implemented)

c     elseif(isw.eq.7) then

c     Compute nodal heat flux for output/plots

      elseif(isw.eq.8) then

        call thcn3d(d,xl,ul,shp,p,s,p(nen+1),ndf,ndm,nel)

c     Compute error data for heat flux

      elseif(isw.eq.11) then

        call ther3d(d,xl,ul,shp,s,ndf,ndm,nel,nen)

c     External node determination

      elseif(isw.eq.26) then

        call pcorner3d()

      endif

c     Formats

2000  format(5x,'F o u r i e r   H e a t   C o n d u c t i o n')

2002  format(a1,20a4//5x,'Element Fluxes'//' Elmt Matl 1-coord  2-coord'
     1            ,'  3-coord      1-flux      2-flux      3-flux'/
     2         37x,'      1-grad      2-grad      3-grad')

2003  format(2i5,0p,3f9.3,1p,3e12.3/37x,1p,3e12.3)

      end

      subroutine thcn3d(d,xl,ul,shp,dt,st,ser,ndf,ndm,nel)

      implicit  none

      include  'cdata.h'
      include  'iodata.h'
      include  'prstrs.h'
      include  'strnum.h'

      integer   ndf,ndm,nel, j,l,lint,ord
      real*8    xsj,xg, sg(4,64),sv(5,16)
      real*8    d(*),dt(*),st(nen,*),ser(*),xl(ndm,*),shp(4,*)
      real*8    xx(3),gradt(3),flux(3),dd(3,3),ul(ndf,*)

      save

c     Lumped projection routine

      if(nel.eq.4) then
        ord =  1
        call tint3dn(nel,lint,sv)
      elseif(nel.eq.10) then
        ord =  2
        l   =  3
        call tint3d(l,lint,sv)
      elseif(nel.eq.11) then
        ord =  2
        call tint3dn(nel,lint,sv)
      else
        ord = 0
        if(nint(d(182)).gt.0) then
          call int3dn(nel,lint,sg)
        else
          l = nint(d(5))
          call int3d(l,lint,sg)
        endif
      endif

      do l = 1,lint
        if(ord.gt.0) then
          call tetshp(sv(1,l),xl,ndm,ord,xsj,shp)
          xsj = xsj*sv(5,l)
        else
          call shp3d(sg(1,l),xsj,shp,xl,ndm,nel)
          xsj = xsj*sg(4,l)
        endif

        call thfx3d(d,xl,ul, xx,shp,gradt,flux,dd, ndm,ndf,nel)

c       Compute lumped projection and assemble stress integrals

        do j = 1,nel
          xg      = xsj*shp(4,j)
          dt(j)   = dt(j) + xg
          st(j,7) = st(j,7) + flux(1)*xg
          st(j,8) = st(j,8) + flux(2)*xg
          st(j,9) = st(j,9) + flux(3)*xg
          ser(j)  = ser(j)  + erav*xg
        end do ! j
      end do ! l

      iste = 9

      end

      subroutine thfx3d(d,xl,ul, xx,shp,gradt,flux,dd, ndm,ndf,nel)

c     Compute thermal gradient and flux

      implicit  none

      integer   ndm,ndf,nel, i
      real*8    psi,cs,sn,c2,s2,d(*),xl(ndm,*),ul(ndf,*), shp(4,*)
      real*8    xx(3),gradt(3),flux(3),dd(3,3)

      save

      gradt(1) = 0.0d0
      gradt(2) = 0.0d0
      gradt(3) = 0.0d0
      xx(1)    = 0.0d0
      xx(2)    = 0.0d0
      xx(3)    = 0.0d0
      do i = 1,nel
        gradt(1) = gradt(1) + shp(1,i)*ul(1,i)
        gradt(2) = gradt(2) + shp(2,i)*ul(1,i)
        gradt(3) = gradt(3) + shp(3,i)*ul(1,i)
        xx(1)    = xx(1) + shp(4,i)*xl(1,i)
        xx(2)    = xx(2) + shp(4,i)*xl(2,i)
        xx(3)    = xx(3) + shp(4,i)*xl(3,i)
      end do ! i

c     Compute thermal flux

      psi = d(31)
      cs  = cos(psi)
      sn  = sin(psi)
      c2  = cs*cs
      s2  = sn*sn
      cs  = cs*sn

      dd(1,1) = c2*d(61) + s2*d(62)
      dd(1,2) = cs*(d(61) - d(62))
      dd(1,3) = 0.0d0

      dd(2,1) = dd(1,2)
      dd(2,2) = s2*d(61) + c2*d(62)
      dd(2,3) = 0.0d0

      dd(3,1) = 0.0d0
      dd(3,2) = 0.0d0
      dd(3,3) = d(63)

      flux(1) = -dd(1,1)*gradt(1) - dd(1,2)*gradt(2)
      flux(2) = -dd(2,1)*gradt(1) - dd(2,2)*gradt(2)
      flux(3) = -dd(3,3)*gradt(3)

      end

      subroutine ther3d(d,xl,ul,shp,st,ndf,ndm,nel,nen)

      implicit  none

      include  'adapt1.h'
      include  'adapt2.h'
      include  'errind.h'

      integer   ndf,ndm,nel,nen, i,j,l,lint, ord
      real*8    xsj, st(nen,*),xl(ndm,*),shp(4,*)
      real*8    d(*),gradt(3),flux(3),dd(3,3),ul(ndf,*),sg(4,64)
      real*8    xx(3),gradp(3),fluxp(3),sv(5,4)

      save

      vfem   = 0.d0
      vproj  = 0.d0
      verror = 0.d0
      vener  = 0.d0
      venere = 0.d0
      heta   = 0.0d0
      if(nel.eq.4) then
        ord = 1
        if(nint(d(182)).gt.0) then
          call tint3dn(nel,lint,sv)
        else
          l =  2
          call tint3d(l,lint,sv)
        endif
      elseif(nel.eq.10) then
        ord =  2
        l   =  3
        call tint3d(l,lint,sv)
      elseif(nel.eq.11) then
        ord =  2
        if(nint(d(182)).gt.0) then
          call tint3dn(nel,lint,sv)
        else
          l =  4
          call tint3d(l,lint,sv)
        endif
      else
        ord = 0
        if(nint(d(182)).gt.0) then
          call int3dn(nel,lint,sg)
        else
          l = nint(d(5))
          call int3d(l,lint,sg)
        endif
      endif

      do l = 1,lint
        if(ord.gt.0) then
          call tetshp(sv(1,l),xl,ndm,ord,xsj,shp)
        else
          call shp3d(sg(1,l),xsj,shp,xl,ndm,nel)
        endif
        call thfx3d(d,xl,ul, xx,shp,gradt,flux,dd, ndm,ndf,nel)
        do i = 1,3
          fluxp(i) = 0.0d0
        end do ! i
        do i = 1,nel
          do j = 1,3
            fluxp(j) = fluxp(j) + shp(4,i)*st(i,j+6)
          end do ! j
        end do ! i

c       Compute integral of stress squares for error indicator use

        call invert(dd,3,3)

        gradp(1) = -dd(1,1)*fluxp(1) - dd(1,2)*fluxp(2)
        gradp(2) = -dd(2,1)*fluxp(1) - dd(2,2)*fluxp(2)

        heta = heta + xsj
        do i = 1,3
         vfem    = vfem   + flux(i)*flux(i)*xsj
         vproj   = vproj  + fluxp(i)*fluxp(i)*xsj
         verror  = verror + ((fluxp(i)-flux(i))**2)*xsj
         vener   = vener  + flux(i)*gradt(i)*xsj
         venere  = venere + (fluxp(i)-flux(i))*(gradp(i)-gradt(i))*xsj
        end do ! i

      end do ! l
      arsq  =  arsq  + heta
      efem  =  efem  + vfem
      eproj =  eproj + vproj
      eerror=  eerror+ verror
      eener =  eener + vener
      eenere=  eenere+ venere
      areai = heta

      heta  =  d(50)*sqrt(heta)

      end
