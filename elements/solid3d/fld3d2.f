c$Id: fld3d2.f,v 1.1 2006/11/20 20:32:56 rlt Exp $
      subroutine fld3d2(d,ul,xl,s,r,ndf,ndm,nst,isw)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c     1. Use constants from 'pconstant.h'                   14/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Finite deformation mixed model: u-p-theta

c      Inputs:
c         d(*)  - Element parameters
c         ul(ndf,*) - Current nodal solution parameters
c         xl(ndm,*) - Nodal coordinates
c         ndf       - Degree of freedoms/node
c         ndm       - Mesh coordinate dimension
c         nst       - Element array dimension
c         isw       - Solution option switch

c      Outputs:
c         s(nst,*)  - Element array
c         r(ndf,*)  - Element vector

c     MATERIAL PROPERTIES

c          imat  = 1 (for compressible J_2/3 regularized neo-Hookean)
c                  2 (for compressible neo-Hookean)
c                  3 (for Ogden-type; elasto/visco/damage)
c          d(4)  = mass density (rho)

c     MODEL 1. (imat = 1 or 2)

c       Properties for Compressible Neo-Hookean (elastic only)

c          d(21)  = K            - bulk  modulus
c          d(22)  = mu           - shear modulus


c          hr( 1)      :  Augmented Lagrangian parameter (xlam1)
c          hr( 2)      :  Constraint fn. for augmented Lagrangian
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'augdat.h'
      include  'bdata.h'
      include  'cdata.h'
      include  'elcoor.h'
      include  'eldata.h'
      include  'elengy.h'
      include  'elplot.h'
      include  'eltran.h'
      include  'hdata.h'
      include  'iofile.h'
      include  'pconstant.h'
      include  'part0.h'
      include  'pmod2d.h'
      include  'prstrs.h'
      include  'rdata.h'
      include  'comblk.h'

      logical   bflg
      integer   ndf,ndm,nst,isw,i,ii,i1,j,jj,j1,l,lint,nhi,nhv,nn, istrt
      integer   ord,npm

      real*8    augfp, d1, epp, thlog, dtheta,qfact, xsj,ta
      real*8    dsigtr, dpress, mpress, dmass, cmshp,lmshp, bdb
      real*8    d(*),  ul(ndf,nen,*),  xl(ndm,*), s(nst,*),r(ndf,*)

      real*8    sg(4,125),sv(5,11),fi(9,2,125),finv(9,125),df(9,125)
      real*8    bbd(3,7),  bei(6), ad(7,7,5,125), dd(7,7), sigm(10)
      real*8    shp(4,27,125), dvol(125),dvol0(125), detf(2,125),xxm(3)
      real*8    sigl(10,125),shpbar(3,27,125),weng(125),theta(2,125)
      real*8    xr(3,125),ur(3,125)
      real*8    xu(3,27),ru(3,27),acc(3),bpra(3),body(3), bf(3),bt(3,3)
      real*8    hh(4,4),hsig(4),phi(4,125),press(125),x0(3)

      save

      data      nhi   / 2 /

c     TEMPORARY TEMPERATURE

      data      ta    / 0.0d0 /

c     Augmented Lagrangian update for nested iteration

      if(isw.eq.10) then
       d1      = augfp*d(21)
       hr(nh2) = hr(nh2) + d1*hr(nh2+1)

c     Compute tangent stiffness and residual force vector

      elseif(isw.eq.3 .or. isw.eq.4  .or. isw.eq.6 .or.
     &       isw.eq.8 .or. isw.eq.14 .or. isw.eq.16) then
       augfp  = augf
       estore = 0.0d0

c      Compute current geometry

       do i = 1,3
        do j = 1,nel
         xu(i,j) = xl(i,j) + ul(i,j,1)
        end do ! j
       end do ! i

c      Set quadrature and order

       if(nel.eq.4) then
         ord = 1
         npm = 1
         if(nint(d(182)).gt.0) then
           call tint3dn(nel,lint,sv)
         else
           l = 1
           call tint3d (l,lint,sv)
         endif
       elseif(nel.eq.10) then
         ord = 2
         npm = 1
         l   = 3
         call tint3d(l,lint,sv)
       elseif(nel.eq.11) then
         ord = 2
         npm = 4
         if(nint(d(182)).gt.0) then
           call tint3dn(nel,lint,sv)
         else
           l =  4
           call tint3d (l,lint,sv)
         endif
       else
         ord = 0
         npm = 1
         if(nint(d(182)).gt.0) then
           call int3dn(nel,lint,sg)
         else
           l = nint(d(5))
           call int3d(l,lint,sg)
         endif
         if(nel.ge.20) then
           npm = 4
         else
           npm = 1
         endif
       endif

c      Mean coordinate of vertex nodes

       if(npm.gt.1) then
         do i = 1,3
           x0(i) = 0.0d0
           do j = 1,nel
             x0(i) = x0(i) + xl(i,j)
           end do ! j
           x0(i) = x0(i)/dble(nel)
         end do ! i
       endif

c      Get shape functions and derivatives in geometry at time t_n+1

       do l = 1,lint
         if(ord.eq.0) then
           call shp3d(sg(1,l),xsj,shp(1,1,l),xl,ndm,nel)
           dvol0(l) = xsj*sg(4,l)
         else
           call tetshp(sv(1,l),xl,ndm,ord,xsj,shp(1,1,l))
           dvol0(l) = xsj*sv(5,l)
         endif
         do i = 1,3
           xr(i,l) = 0.0d0
           ur(i,l) = 0.0d0
           do j = 1,nel
             xr(i,l) = xr(i,l) + shp(4,j,l)*xl(i,l)
             ur(i,l) = ur(i,l) + shp(4,j,l)*ul(i,l,1)
           end do ! j
c          ur(i,l) = ur(i,l) + xr(i,l)
         end do ! j

         phi(1,l) = 1.d0
         if(npm.gt.1) then
           phi(2,l) = xr(1,l) - x0(1)
           phi(3,l) = xr(2,l) - x0(2)
           phi(4,l) = xr(3,l) - x0(3)
         endif
       end do ! l

c      Set number of history terms / quadradure point

       nhv   = nint(d(15))
       istrt = nint(d(84))

c      MECHANICAL ELEMENT

       if(isw.eq.3 .or. isw.eq. 6 .or. isw.eq.14) then

c       Compute f, finv, df and det(fei) at conf t-n+1

        call kine3m(shp,ul,fi,finv,df,detf,ndf,nel,nen,lint)

c       Compute volume at current state

        do l = 1,lint
          dvol(l) = dvol0(l)*detf(1,l)
        end do ! l

c       Mixed model for volumetric response

        call bbar3m(phi,shp,dvol,detf,lint,nel,npm,hh,theta,shpbar)

c       Compute mixed model deformation gradient

        call fbar3m(fi,detf,theta,lint)

c       Compute Cauchy stresses and spatial tangent tensor at t-n+1

        nn = nhi
        do l = 1,lint

         do i = 1,3
           xref(i) = xr(i,l)
           xcur(i) = xr(i,l) + ur(i,l)
         end do ! i

         call modlfd(d,fi(1,1,l),finv(1,l),df(1,l),theta(1,l),ta,
     &               hr(nn+nh1),hr(nn+nh2),nhv,istrt,ad(1,1,1,l),
     &               sigl(1,l),bei,hr(nh2),hr(nh2+1),.true.,isw)

         nn = nn + nhv
        end do ! l

        if(isw.eq.14) return

c       Compute mixed pressure

        if(isw.eq.3 .or. isw.eq.6) then

         call sbodyf(d, body) ! Compute body force values

         if(npm.eq.1) then

           press(1) = 0.0d0
           do l = 1,lint

c           Modify volume element and integrate pressure
c           over reference volume

            press(1) = press(1) + one3*(sigl(1,l) + sigl(2,l)
     &                                 + sigl(3,l))*dvol0(l)
            dvol(l)  = dvol0(l) * theta(1,l)

           end do ! l

c          Divide pressure by reference volume

           press(1) = press(1) * hh(1,1)
           do l = 2,lint
             press(l) = press(1)
           end do ! l
         else
           do i = 1,npm
             sigm(i) = 0.0d0
           end do ! i

           do l = 1,lint

c            Modify volume element and integrate pressure
c            over reference volume

             mpress   = one3*(sigl(1,l) + sigl(2,l)
     &                       + sigl(3,l))*dvol0(l)
             sigm(1) = sigm(1) + mpress
             do j = 2,npm
               sigm(j) = sigm(j) + mpress*phi(j,l)
             end do ! j
             dvol(l) = dvol0(l) * theta(1,l)

           end do ! l

c          Divide pressure by reference volume

           do i = 1,npm
             hsig(i) = 0.0d0
             do j = 1,npm
               hsig(i) = hsig(i) + hh(i,j)*sigm(j)
             end do ! j
           end do ! i
           do l = 1,lint
             press(l) = hsig(1)
             do j = 2,npm
               press(l) = press(l) + hsig(j)*phi(j,l)
             end do ! j
           end do ! l

         endif
         bflg  = d(4).gt.0.0d0 .and. d(65).gt.0.0d0

c        Compute final residual and tangent arrays

         do l = 1,lint

c         Angular velocity: d(4) = rho; d(65) = omega

          do i = 1,3
            bf(i) = 0.0d0
          end do ! i
          if(bflg) then
            call sbodyw(d(4),d(65),ur(1,l), bf,bt, .true.)
          endif

c         Compute mixed stress and multiply by volume element

          dsigtr  = press(l)*detf(1,l)/theta(1,l)
     &            - (sigl(1,l)+sigl(2,l)+sigl(3,l))*one3

          sigm(1) =  sigl(1,l) + dsigtr
          sigm(2) =  sigl(2,l) + dsigtr
          sigm(3) =  sigl(3,l) + dsigtr
          sigm(4) =  sigl(4,l)
          sigm(5) =  sigl(5,l)
          sigm(6) =  sigl(6,l)

c         Store time history plot data for element

          i = 6*(l-1)
          do j = 1,6
           tt(j+i) = sigm(j)
           sigm(j) = sigm(j)*dvol(l)
          end do ! j

c         Compute acceleration

          if(d(7).ge.0.0 .and. (ndfo(1).gt.0 .or. shflg)) then
           dmass = d(4)*dvol0(l)
          else
           dmass = 0.0d0
          endif ! d(7) test

          cmshp = dmass*d(7)      ! Consistent mass factor
          lmshp = dmass - cmshp   ! Lumped     mass factor
          do i = 1,3
           acc(i) = 0.0d0
           do j = 1,nel
            acc(i) = acc(i) + shp(4,j,l)*ul(i,j,5)
           end do ! j
           acc(i) = acc(i)*cmshp
          end do ! i

c         Compute residual

          do j = 1,nel

            ru(1,j) = shp(1,j,l)*sigm(1)
     &              + shp(2,j,l)*sigm(4)
     &              + shp(3,j,l)*sigm(6)

            ru(2,j) = shp(1,j,l)*sigm(4)
     &              + shp(2,j,l)*sigm(2)
     &              + shp(3,j,l)*sigm(5)
 
            ru(3,j) = shp(1,j,l)*sigm(6)
     &              + shp(2,j,l)*sigm(5)
     &              + shp(3,j,l)*sigm(3)

            do i = 1,3
              r(i,j)  = r(i,j) + shp(4,j,l)*(dvol0(l)*(body(i) + bf(i))
     &                         - (acc(i) + ul(i,j,5)*lmshp))  - ru(i,j)
            end do ! i
          end do ! j

c         Compute mixed tangent stiffness matrix

          if(isw.eq.3) then

c          Part 1: Geometric tangent matrix

           if(gflag) then

             i1 = 0
             do i = 1,nel
               j1 = 0
               do j = 1,nel
                 bdb = (shp(1,i,l)*ru(1,j)
     &               +  shp(2,i,l)*ru(2,j)
     &               +  shp(3,i,l)*ru(3,j))*ctan(1)
                 do jj = 1,3
                   s(i1+jj,j1+jj) = s(i1+jj,j1+jj) + bdb
                 end do ! jj
                 j1 = j1 + ndf
               end do ! j
               i1 = i1 + ndf
             end do ! i

           endif

c          Part 2: Material tangent matrix

c          Modify tangent moduli for stress factors

           mpress = press(l)*detf(1,l)/theta(1,l)
           dpress = one3*(sigl(1,l) + sigl(2,l) + sigl(3,l))

           call dmatdx(ad(1,1,1,l),sigl(1,l),dpress,mpress)

c          Multiply tangent moduli by volume element

           d1 = dvol(l)*ctan(1)
           do i = 1,7
             do j = 1,7
               dd(i,j) = ad(i,j,1,l)*d1
             end do ! j
           end do ! i

c          Compute row terms

           dmass = ctan(3)*dmass
           i1    = 0
           do i = 1,nel

c            Angular velocity tangent

             if(bflg) then
               do jj = 1,3
                 do ii = 1,3
                   bdb = shp(4,i,l)*bt(ii,jj)
                   j1  = 0
                   do j = 1,i
                     s(i1+ii,j1+jj) = s(i1+ii,j1+jj) + bdb*shp(4,j,l)
                     j1             = j1 + ndf
                   end do ! j
                 end do ! ii
               end do ! jj
             endif

c           Compute bmat-t * dd * dvol

            do jj = 1,7

             bbd(1,jj) =    shp(1,i,l)*dd(1,jj)
     &                 +    shp(2,i,l)*dd(4,jj)
     &                 +    shp(3,i,l)*dd(6,jj)
     &                 + shpbar(1,i,l)*dd(7,jj)

             bbd(2,jj) =    shp(2,i,l)*dd(2,jj)
     &                 +    shp(1,i,l)*dd(4,jj)
     &                 +    shp(3,i,l)*dd(5,jj)
     &                 + shpbar(2,i,l)*dd(7,jj)

             bbd(3,jj) =    shp(3,i,l)*dd(3,jj)
     &                 +    shp(2,i,l)*dd(5,jj)
     &                 +    shp(1,i,l)*dd(6,jj)
     &                 + shpbar(3,i,l)*dd(7,jj)
            end do ! jj

            lmshp = shp(4,i,l)*dmass
            cmshp = lmshp*d(7)             ! Consistent mass
            lmshp = lmshp - cmshp          ! Lumped     mass
            do jj = 1,3
              s(i1+jj,i1+jj) = s(i1+jj,i1+jj) + lmshp
            end do ! jj

            j1 = 0
            do j = 1,i

c            Inertial tangent

             do jj = 1,3
              s(i1+jj,j1+jj) = s(i1+jj,j1+jj) + cmshp*shp(4,j,l)
             end do ! jj

c            Compute mechanics part of tangent stiffness

             do jj = 1,3

              s(i1+jj,j1+1) = s(i1+jj,j1+1) + bbd(jj,1)*shp(1,j,l)
     &                                      + bbd(jj,4)*shp(2,j,l)
     &                                      + bbd(jj,6)*shp(3,j,l)
     &                                      + bbd(jj,7)*shpbar(1,j,l)

              s(i1+jj,j1+2) = s(i1+jj,j1+2) + bbd(jj,2)*shp(2,j,l)
     &                                      + bbd(jj,4)*shp(1,j,l)
     &                                      + bbd(jj,5)*shp(3,j,l)
     &                                      + bbd(jj,7)*shpbar(2,j,l)

              s(i1+jj,j1+3) = s(i1+jj,j1+3) + bbd(jj,3)*shp(3,j,l)
     &                                      + bbd(jj,5)*shp(2,j,l)
     &                                      + bbd(jj,6)*shp(1,j,l)
     &                                      + bbd(jj,7)*shpbar(3,j,l)

             end do ! jj

             j1 = j1 + ndf
            end do ! j
            i1 = i1 + ndf
           end do ! i
          endif ! isw = 3
         end do ! l

c        Compute lower part by symmetry

         if(isw .eq. 3) then
          do i = 1,nst
           do j = 1,i
            s(j,i) = s(i,j)
           end do ! j
          end do ! i
         endif
        endif ! isw = 3 or 6

c      Output stresses and computre fracture force

       elseif(isw.eq.4 .or. isw.eq.8 .or. isw.eq.16) then

        do i = 1,10
          sigm(i) = 0.0d0
        end do ! i
        do i = 1,3
          bpra(i) = 0.0d0
          xxm(i)  = 0.0d0
        end do ! i
        epp = 0.0d0
        dtheta = 0.0d0
        qfact  = 1.d0/dble(lint)

c       Compute f, finv, df and det(fei) at conf t-n+1

        call kine3m(shp,ul,fi,finv,df,detf,ndf,nel,nen,lint)

        call bbar3m(phi,shp,dvol,detf,lint,nel,npm,hh,theta,shpbar)

        call fbar3m(fi,detf,theta,lint)

c       Second loop over Gauss points

        nn  = nhi
        do l = 1,lint

c        Compute Cauchy stresses and spatial tangent tensor at t-n+1

         do i = 1,3
           xref(i) = xr(i,l)
           xcur(i) = xr(i,l) + ur(i,l)
         end do ! i

         call modlfd(d,fi(1,1,l),finv(1,l),df(1,l),theta(1,l),ta,
     &               hr(nn+nh1),hr(nn+nh2),nhv,istrt,ad,sigl(1,l),
     &               bei,hr(nh2),hr(nh2+1),.true.,isw)
         weng(l)   = estore

c        Compute principal stretches

         call pstr3d(bei, bpr)

c        Average stresses and stretches for printing

         do i = 1,3
          bpra(i) = bpra(i) + 0.5d0*qfact*log(bpr(i))
          xxm(i)  = xxm(i)  + qfact*xu(i,l)
         end do ! i
         do i = 1,6
          sigm(i) = sigm(i) + qfact*sigl(i,l)
         end do ! i
         sigm(10) = sigm(10) + qfact*sigl(10,l)
         epp      = epp      + qfact*sigl( 9,l)
         dtheta   = dtheta   + qfact*theta(1,l)
         nn = nn + nhv
        end do ! l

c       Output stresses

        if (isw .eq. 4) then

          call pstr3d(sigm,sigm(7))

          mct = mct - 2
          if(mct.le.0) then
           write(iow,2001) o,head
           if(ior.lt.0) write(*,2001) o,head
           mct = 50
          endif

          thlog = log(abs(dtheta))
          write(iow,2002) n,ma,(sigm(i),i=1,9),bpra,
     &                    xxm,thlog,epp,sigm(10)
          if(ior.lt.0) then
            write(*,2002) n,ma,(sigm(i),i=1,9),bpra,
     &                    xxm,thlog,epp,sigm(10)
          endif
        elseif(isw.eq.8) then

c         Project stresses onto nodes

          call slcn3d(sigl,shp,dvol, r,s, lint,nel,27)

c       Compute fracture indices

        elseif(isw.eq.16) then

          call pfrac3f(fi,detf,sigl,weng, shp,dvol, r,
     &               lint,ndf,ndm,3)

        endif ! isw = 4 or 8 or 16

       endif ! isw = 3 or 6 or 4 or 8 or 14 or 16

      endif ! isw tests

c     Formats for input-output

2001  format(a1,20a4//5x,'Element Stresses'//'  Elmt  Matl',
     &   '  11-stress  22-stress  33-stress  12-stress',
     &   '  23-stress  13-stress'/12x,
     &   '   1-stress   2-stress   3-stress',
     &   '  log(lam1)  log(lam2)  log(lam3)'/12x,
     &   '    1-coord    2-coord    3-coord    log-J     eff-ep',
     &   '     Yield')

2002  format(2i6,1p6e11.3/12x,1p6e11.3/12x,0p3f11.5,1p3e11.3/
     &       23x,1p5e11.3,0p,1f11.3/1x)

      end
