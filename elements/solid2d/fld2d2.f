c$Id: fld2d2.f,v 1.1 2006/11/20 20:32:56 rlt Exp $
      subroutine fld2d2(d,ul,xl,ix,s,r,ndf,ndm,nst,isw)

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
c         ix(*)     - Global nodal connections
c         ndf       - Degree of freedoms/node
c         ndm       - Mesh coordinate dimension
c         nst       - Element array dimension
c         isw       - Solution option switch

c      Outputs:
c         s(nst,*)  - Element array
c         r(ndf,*)  - Element vector
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
      include  'part0.h'
      include  'pconstant.h'
      include  'pmod2d.h'
      include  'ptdat6.h'
      include  'rdata.h'
      include  'comblk.h'

      logical   dynflg, quad
      integer   ndf,ndm,nst, isw, i,i1,is,istrt, j,jj,j1,js, k, l,lint
      integer   nhi,nhv,nn, npm,ix(*)
      real*8    augfp,d1,bdb,bd3,epp,thlog,xsj(16)
      real*8    ta,qfact,dsigtr,dpress,mpress,dmass,cmshp,lmshp,dtheta
      real*8    d(*),         ul(ndf,nen,*),   xl(ndm,*),    s(nst,*)
      real*8    sg(3,25),     df(9,25),        fi(9,2,25),   finv(9,25)
      real*8    xr(2,25),     xu(2,16),        ru(3,16),     r(ndf,*)
      real*8    bbd(3,7),     bei(6),          ad(7,7,5,25), dd(7,7)
      real*8    shp(3,16,25), shpbar(2,16,25), dvol0(25),    detf(2,25)
      real*8    sigm(9),      sigl(16,25),     bpra(3),      dvol(25)
      real*8    theta(2,25),  hh(6,6),         el(4,7),      vl(3)
      real*8    press(25),    hsig(6),         weng(25),     body(3)
      real*8    ur(2,25),     phi(6,25),       shpr(16,25),  x0(2)

      save

      data    nhi   / 2 /

c     TEMPORARY SET OF TEMPERATURE

      data    ta    / 0.0d0 /

c     Data input

      if(isw.eq.1) then

c     Augmented Lagrangian update for nested iteration

      elseif(isw.eq.10) then

        d1      = augfp*d(21)
        hr(nh2) = hr(nh2) + d1*hr(nh2+1)

c     Compute tangent stiffness and residual force vector

      elseif(isw.eq. 3 .or. isw.eq. 4 .or. isw.eq. 6 .or.
     &       isw.eq. 8 .or. isw.eq.13 .or. isw.eq.14 .or.
     &       isw.eq.16 .or. isw.eq.25) then

        augfp  = augf
        estore = 0.0d0

c       Compute current geometry

        do j = 1,nel
          do i = 1,2
            xu(i,j) = xl(i,j) + ul(i,j,1)
          end do ! i
        end do ! j

c       Set element quadrature and number of pressure/volume modes

        if(nel.eq.3) then
          if(d(182).gt.0.0d0) then
            call tint2dn(nel,lint,el)
          else
            l =  1
            call tint2d(l,lint,el)
          endif
          npm =  1
          quad = .false.
        elseif(nel.eq.6 .or. nel.eq.7) then
          if(d(182).gt.0.0d0) then
            call tint2dn(nel,lint,el)
          else
            l =  7
            call tint2d(l,lint,el)
          endif
          if(nel.eq.6) then
            npm =  1
          else
            npm =  3
          endif
          quad = .false.
        else
          if(nel.le.4) then
            l   = 2
            npm = 1
          elseif(nel.le.9) then
            l   = 3
            npm = 3
          else
            l   = 4
            npm = 6
          endif
          if(nint(d(182)).gt.0) then
            call int2dn(nel,lint,sg)
          else
            call int2d(l,lint,sg)
          endif
          quad = .true.
        endif

c       Center coordinate

        if(npm.gt.1) then
          x0(1) = 0.0d0
          x0(2) = 0.0d0
          do j = 1,nel
            x0(1) = x0(1) + xl(1,j)
            x0(2) = x0(2) + xl(2,j)
          end do ! j
          x0(1) = x0(1)/dble(nel)
          x0(2) = x0(2)/dble(nel)
        endif

c       Get shape functions and derivatives in geometry at time t_n+1

        do l = 1,lint
          if(quad) then
            call shp2d(sg(1,l),xl,shp(1,1,l),xsj(l),ndm,nel,ix,.false.)
            dvol0(l) = xsj(l)*sg(3,l)
          else
            call trishp(el(1,l),xl,ndm,nel-4,xsj(l),shp(1,1,l))
            dvol0(l) = xsj(l)*el(4,l)
          endif

c         Compute coordinates at gauss points

          xr(1,l) = 0.0d0
          xr(2,l) = 0.0d0
          ur(1,l) = 0.0d0
          ur(2,l) = 0.0d0
          do i = 1,nel
            xr(1,l) = xr(1,l) + xl(1,i)  *shp(3,i,l)
            xr(2,l) = xr(2,l) + xl(2,i)  *shp(3,i,l)
            ur(1,l) = ur(1,l) + ul(1,i,1)*shp(3,i,l)
            ur(2,l) = ur(2,l) + ul(2,i,1)*shp(3,i,l)
          end do ! i

c         Axisymmetric volume

          if(stype.eq.3 .or. stype.eq.8) then
            dvol0(l) = dvol0(l)*xr(1,l)
            do i = 1,nel
              shpr(i,l)  = shp(3,i,l)/(xr(1,l)+ur(1,l))
            end do ! i
          else
            do i = 1,nel
              shpr(i,l)  = 0.0d0
            end do ! i
          endif

          phi(1,l) = 1.0d0
          if(npm.gt.1) then
            phi(2,l) = xr(1,l) - x0(1)
            phi(3,l) = xr(2,l) - x0(2)
            phi(4,l) = phi(2,l)**2
            phi(5,l) = phi(2,l)*phi(3,l)
            phi(6,l) = phi(3,l)**2
          endif

        end do ! l
        xref(3) = 0.0d0
        xcur(3) = 0.0d0

c       Set number of history terms / quadradure point

        nhv   = nint(d(15))
        istrt = nint(d(84))

c       Set for torsion terms

        if(stype.eq.8) then
          is    = 3
          js    = 6
        else
          is    = 2
          js    = 4
          vl(3) = 0.0d0
        endif

c       MECHANICAL ELEMENT

        if(isw.eq.3 .or. isw.eq.6 .or. isw.eq.13 .or. isw.eq.14) then

c         Compute f, finv, df, detf and shp at conf t-n+1

          call kine2d(shp,xl,ul,fi,finv,df,detf,ndm,ndf,nel,nen,lint)

c         Compute volume at current state

          do l = 1,lint
            dvol(l) = dvol0(l)*detf(1,l)
          end do ! l

c         Mixed model for volumetric response

          call bbar2m(phi,shp,shpr,dvol,detf,lint,nel,npm,
     &                hh,theta,shpbar)

c         Compute mixed model deformation gradient

          call fbar2m(fi,detf,theta,lint)

c         Compute Cauchy stresses and spatial tangent tensor at t-n+1

          dynflg = ctan(3).ne.0.0d0
          nn     = nhi
          do l = 1,lint

c           Compute coordinates in reference and current configuration

            do i = 1,2
              xref(i) = xr(i,l)
              xcur(i) = xr(i,l) + ur(i,l)
            end do ! i

            call modlfd(d,fi(1,1,l),finv(1,l),df(1,l),theta(1,l),ta,
     &                  hr(nn+nh1),hr(nn+nh2),nhv,istrt,ad(1,1,1,l),
     &                  sigl(1,l),bei,hr(nh2),hr(nh2+1),.true.,isw)
            nn = nn + nhv

c           Accumulate the energy

            if(isw.eq.13) then

              epl(8) = epl(8) + estore*dvol0(l)

c             Compute velocity at point

              do i = 1,is
                vl(i) = 0.0d0
                do j = 1,nel
                  vl(i) = vl(i) + ul(i,j,4)*shp(3,j,l)
                end do ! j
              end do ! i

              d1 = 0.0d0
              if(stype.eq.8) then
                do i = 1,nel
                  d1 = d1 + (ul(1,i,4)**2 +  ul(2,i,4)**2
     &                    +  ul(3,i,4)**2)*shp(3,i,l)
                end do ! i
              else
                do i = 1,nel
                  d1 = d1 + (ul(1,i,4)**2 +  ul(2,i,4)**2)*shp(3,i,l)
                end do ! i
              endif

c             Accumulate kinetic energy

              epl(7) = epl(7) + 0.5d0*((1.d0-d(7))*d1
     &                        + d(7)*(vl(1)**2
     &                              + vl(2)**2
     &                              + vl(3)**2))*dvol0(l)*d(4)
            endif ! isw.eq.13

          end do ! l

          if(isw.eq.14) return

c         Compute mixed pressure

          if(isw.eq.3 .or. isw.eq.6) then

c           Set body load levels

            call sbodyf(d, body)

            if(npm.eq.1) then

              press(1) = 0.0d0
              do l = 1,lint

c               Modify volume element and integrate pressure
c               over reference volume

                press(1) = press(1) + one3*(sigl(1,l) + sigl(2,l)
     &                              + sigl(3,l))*dvol0(l)
                dvol(l)  = dvol0(l) * theta(1,l)

              end do ! l

c             Divide pressure by reference volume

              press(1) = press(1) * hh(1,1)
              do l = 2,lint
                press(l) = press(1)
              end do ! l

            else

              do i = 1,npm
                sigm(i) = 0.0d0
              end do ! i

              do l = 1,lint

c               Modify volume element and integrate pressure
c               over reference volume

                mpress   = one3*(sigl(1,l) + sigl(2,l)
     &                          + sigl(3,l))*dvol0(l)
                sigm(1) = sigm(1) + mpress
                do k = 2,npm
                  sigm(k) = sigm(k) + mpress*phi(k,l)
                end do ! k
                dvol(l) = dvol0(l) * theta(1,l)

              end do ! l

c             Divide pressure by reference volume

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

c           Set consistent/lumped factor

            if(stype.eq.8) then
              cmshp = 1.d0            ! Axi+tor consistent
            else
              cmshp = d(7)            ! Consistent
            endif
            lmshp   = 1.d0 - cmshp    ! Lumped

c           Compute final residual and tangent arrays

            do l = 1,lint

c             Compute mixed stress and multiply by volume element

              dsigtr  = press(l)*detf(1,l)/theta(1,l)
     &                - (sigl(1,l)+sigl(2,l)+sigl(3,l))*one3
              sigm(1) =  sigl(1,l) + dsigtr
              sigm(2) =  sigl(2,l) + dsigtr
              sigm(3) =  sigl(3,l) + dsigtr
              sigm(4) =  sigl(4,l)
              sigm(5) =  sigl(5,l)
              sigm(6) =  sigl(6,l)

c             Store time history plot data for element

              i = 6*(l-1)
              do j = 1,js
                tt(j+i) = sigm(j)
                sigm(j) = sigm(j)*dvol(l)
              end do ! j

c             Compute mass factor

              if(d(7).ge.0.0 .and. (ndfo(1).gt.0 .or. shflg)) then
                dmass = d(4)*dvol0(l)
              else
                dmass = 0.0d0    ! No inertia effects
              endif

c             Compute residual

              do j = 1,nel

                ru(1,j) = shp(1,j,l)*sigm(1) + shp(2,j,l)*sigm(4)
                ru(2,j) = shp(1,j,l)*sigm(4) + shp(2,j,l)*sigm(2)

                r(1,j)  = r(1,j) - ru(1,j) - shpr(j,l)*sigm(3)
     &                           + body(1)*shp(3,j,l)*dvol0(l)
                r(2,j)  = r(2,j) - ru(2,j)
     &                           + body(2)*shp(3,j,l)*dvol0(l)

              end do ! j

              if(stype.eq.8) then
                do j = 1,nel
                  ru(3,j) = xcur(1)*(shp(1,j,l)*sigm(6)
     &                             + shp(2,j,l)*sigm(5))
                  r(3,j)  = r(3,j) - ru(3,j)
     &                             + body(3)*shp(3,j,l)*dvol0(l)
                  ru(3,j) = ru(3,j)*2.d0
                end do ! j
              endif

c             Compute mixed tangent stiffness matrix

              if(isw.eq.3) then

c               Part 1: Geometric tangent matrix

                if(gflag) then
                  i1 = 0
                  do i = 1,nel
                    bd3 = shpr(i,l)*sigm(3)*ctan(1)
                    j1 = 0
                    do j = 1,nel
                      bdb = (shp(1,i,l)*ru(1,j)
     &                    +  shp(2,i,l)*ru(2,j))*ctan(1)
                      s(i1+1,j1+1) = s(i1+1,j1+1) + bdb + bd3*shpr(j,l)
                      s(i1+2,j1+2) = s(i1+2,j1+2) + bdb

c                     Torsion part

                      if(stype.eq.8) then
                        s(i1+1,j1+3) = s(i1+1,j1+3)
     &                               + shpr(i,l)*ru(3,j)*ctan(1)
                        s(i1+3,j1+1) = s(i1+3,j1+1)
     &                               + shpr(j,l)*ru(3,i)*ctan(1)
                        s(i1+3,j1+3) = s(i1+3,j1+3)+bdb*xcur(1)*xcur(1)
                      endif
                      j1 = j1 + ndf
                    end do ! j

                    i1 = i1 + ndf
                  end do ! i
                endif ! gflag

c               Part 2: Material tangent matrix

c               Modify tangent moduli for stress factors

                mpress = press(l)*detf(1,l)/theta(1,l)
                dpress = one3*(sigl(1,l) + sigl(2,l) + sigl(3,l))

                call dmatdx(ad(1,1,1,l),sigl(1,l),dpress,mpress)

c               Multiply tangent moduli by volume element

                d1 = dvol(l)*ctan(1)
                do i = 1,7
                  do j = 1,7
                    dd(i,j) = ad(i,j,1,l)*d1
                  end do ! j
                end do ! i

c               Compute row terms

                dmass = ctan(3)*dmass
                i1    = 0
                do i = 1,nel

c                 Compute bmat-t * dd * dvol

                  do jj = 1,7

                    bbd(1,jj) =    shp(1,i,l)*dd(1,jj)
     &                        +    shpr( i,l)*dd(3,jj)
     &                        +    shp(2,i,l)*dd(4,jj)
     &                        + shpbar(1,i,l)*dd(7,jj)

                    bbd(2,jj) =    shp(2,i,l)*dd(2,jj)
     &                        +    shp(1,i,l)*dd(4,jj)
     &                        + shpbar(2,i,l)*dd(7,jj)
                  end do ! jj

c                 Torsion part

                  if(stype.eq.8) then
                    do jj = 1,7
                      bbd(3,jj) = xcur(1)*(shp(2,i,l)*dd(5,jj)
     &                                   + shp(1,i,l)*dd(6,jj))
                    end do ! jj
                  endif

                  j1 = 0
                  do j = 1,nel

c                   Compute mechanics part of tangent stiffness

                    do jj = 1,is
                      s(i1+jj,j1+1) = s(i1+jj,j1+1)
     &                              + bbd(jj,1)*shp(1,j,l)
     &                              + bbd(jj,3)*shpr( j,l)
     &                              + bbd(jj,4)*shp(2,j,l)
     &                              + bbd(jj,7)*shpbar(1,j,l)

                      s(i1+jj,j1+2) = s(i1+jj,j1+2)
     &                              + bbd(jj,2)*shp(2,j,l)
     &                              + bbd(jj,4)*shp(1,j,l)
     &                              + bbd(jj,7)*shpbar(2,j,l)
                    end do ! jj

c                   Torsion part

                    if(stype.eq.8) then
                      do jj = 1,is
                        s(i1+jj,j1+3) = s(i1+jj,j1+3)
     &                                + (bbd(jj,5)*shp(2,j,l)
     &                                +  bbd(jj,6)*shp(1,j,l))*xcur(1)
                      end do ! jj
                    endif

                    j1 = j1 + ndf
                  end do ! j
                  i1 = i1 + ndf
                end do ! i
              endif ! isw = 3

c             Add inertia parts

              if(dynflg) then
                call fdyn2d(ul,shp(1,1,l),s,r,is,xcur(1),
     &                      cmshp,lmshp,dmass,isw)
              endif
            end do ! l

c           Multiply by thickness if not unity

            if((isw.eq.3 .or.isw.eq.6) .and. d(14).ne.1.d0) then

              do j = 1,nst
                do i = 1,nst
                  s(i,j) = s(i,j)*d(14)
                end do ! i
              end do ! j
              do j = 1,nel
                do i = 1,ndf
                  r(i,j) = r(i,j)*d(14)
                end do ! i
              end do ! j

            endif

          endif ! isw = 3 or 6

c       Output stresses.

        elseif(isw.eq.4.or.isw.eq.8.or.isw.eq.16.or.isw.eq.25) then

          do i = 1,9
            sigm(i) = 0.0d0
          end do ! i
          do i = 1,3
            bpra(i) = 0.0d0
          end do ! i
          epp    = 0.0d0
          dtheta = 0.0d0
          qfact  = 1.d0/dble(lint)

c         Compute f, finv, df and det(fei) at conf t-n+1

          call kine2d(shp,xl,ul,fi,finv,df,detf,ndm,ndf,nel,nen,lint)

          call bbar2m(phi,shp,shpr,dvol,detf,lint,nel,npm,
     &                hh,theta,shpbar)

          call fbar2m(fi,detf,theta,lint)

c         Second loop over Gauss points

          nn  = nhi
          do l = 1,lint

c           Compute coordinates in reference and current configuration

            do i = 1,2
              xref(i) = xr(i,l)
              xcur(i) = xr(i,l) + ur(i,l)
            end do ! i

c           Compute Cauchy stresses and spatial tangent tensor at t-n+1

            call modlfd(d,fi(1,1,l),finv(1,l),df(1,l),theta(1,l),ta,
     &                  hr(nn+nh1),hr(nn+nh2),nhv,istrt,ad,sigl(1,l),
     &                  bei,hr(nh2),hr(nh2+1),.true.,isw)
            weng(l) = estore

c           Compute principal stretches

            call pstr3d(bei, bpr)

c           Average stresses and stretches for printing

            do i = 1,3
              bpra(i) = bpra(i) + 0.5d0*qfact*log(bpr(i))
            end do ! i
            do i = 1,js
              sigm(i) = sigm(i) + qfact*sigl(i,l)
            end do ! i
            epp      = epp      + qfact*sigl( 9,l)
            dtheta   = dtheta   + qfact*theta(1,l)
            nn = nn + nhv
          end do ! l

c         Output stresses

          if (isw .eq. 4) then

            if(stype.eq.8) then
              call pstr3d(sigm,sigm(7))
            else
              call pstr2d(sigm,sigm(7))
            endif

            mct = mct - 2
            if(mct.le.0) then
              if(stype.eq.8) then
                write(iow,2001) o,head,'3-Stress'
                if(ior.lt.0) write(*,2001) o,head,'3-Stress'
              else
                write(iow,2001) o,head,'   Angle'
                if(ior.lt.0) write(*,2001) o,head,'   Angle'
              endif
              mct = 50
            endif

c           Compute potential damage variable

            thlog = log(abs(dtheta))
            write(iow,2002) n,ma,(sigm(i),i=1,9),bpra,
     &                      xcur(1),xcur(2),thlog,epp
            if(ior.lt.0) then
              write(*,2002) n,ma,(sigm(i),i=1,9),bpra,
     &                      xcur(1),xcur(2),thlog,epp
            endif

c         Project stresses onto nodes

          elseif(isw.eq.8) then

            call slcn2d(ix,sigl,shp,xsj,r,s,r(nen+1,1),lint,nel,16)

c         Compute fracture indices

          elseif(isw.eq.16) then

            call pfrac2f(fi,theta,sigl,weng, shp,dvol, r,
     &                   lint,ndf,ndm,3)

c         Compute Z-Z projections

          elseif(isw.eq.25) then

            call stcn2z(xl,sigl,shp,dvol,lint,ndm,nel,16)

          endif
        endif ! isw = 4 or 8 or 16 or 25

      endif ! isw

c     Formats

2001  format(a1,20a4//5x,'Element Stresses'//'  Elmt  Matl',
     &   '  11-Stress  22-Stress  33-Stress  12-Stress  23-Stress',
     &   '  13-Stress'/12x,'   1-Stress   2-Stress   ',a8,
     &   '  log(lam1)  log(lam2)  log(lam3)'/12x,'    1-Coord',
     &   '    2-Coord    log-J     eff-ep')

2002  format(2i6,1p6e11.3/12x,1p6e11.3/12x,0p2f11.5,1p2e11.3/1x)

      end
