c$Id: sld3d2.f,v 1.2 2006/12/06 04:59:12 rlt Exp $
      subroutine sld3d2(d,ul,xl,tl,s,r,ndf,ndm,nst,isw, ther)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c     1. Use constants from 'pconstant.h'                   14/11/2006
c     2. Add rotational body forces                         05/12/2006
c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose: Three-dimensional u-p-theta small deformation element

c      Inputs:
c         d(*)  - Element parameters
c         ul(ndf,*) - Current nodal solution parameters
c         xl(ndm,*) - Nodal coordinates
c         tl(*)     - Nodal temp vector
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
      include  'elbody.h'
      include  'eldata.h'
      include  'elengy.h'
      include  'elplot.h'
      include  'eltran.h'
      include  'fdata.h'
      include  'hdata.h'
      include  'iofile.h'
      include  'pconstant.h'
      include  'part0.h'
      include  'pmod2d.h'
      include  'ptdat6.h'
      include  'rdata.h'
      include  'comblk.h'

      logical   ther,quad
      integer   ndf,ndm,nst,isw,i,i1,j,jj,j1,l,lint
      integer   ord,nhi,nhv,istrt,nn,npm
      real*8    augfp,dv,dl,d1,xsj,bt(3), aj0,aj1
      real*8    dsigtr, mpress, dvs, dmshp, cfac,lfac, fac
      real*8    d(*),       ul(ndf,nen,*),  xl(ndm,*),   tl(*), s(nst,*)
      real*8    sg(4,64),   r(ndf,*),       xx(3,27) ,   shp(4,27,64)
      real*8    bbd(3,7),   aa(6,6,5,64),   dd(7,7)  ,   dvol(64)
      real*8    sigm(9),    sigl(10,64),    bbar(3,27,64)
      real*8    al(3),      ac(3),          vl(3)    ,   x0(3)
      real*8    phi(4,64),  theta(3,64),    hh(4,4)  ,   el(5,11)
      real*8    press(64),  pbar(64),       hsig(4)  ,   eps(9,3,64)
      real*8    ta(64),     epsm(9),        epsv(64) ,   mom(3,27)
      real*8    dther(7),   body(3),        bf(3)    ,   btan(3,3)

      save

      data    nhi   / 2 /

c     TEMPORARY SET OF TEMPERATURE

      data    ta    / 64*0.0d0 /

c     Data inputs

      if( isw.eq.1 ) then

c     Augmented Lagrangian update for nested iteration

      elseif(isw.eq.10) then

        d1      = augfp*d(21)
        hr(nh2) = hr(nh2) + d1*hr(nh2+1)

c     Compute tangent stiffness and residual force vector

      elseif(isw.eq. 3 .or. isw.eq. 4 .or. isw.eq. 6 .or.
     &       isw.eq. 8 .or. isw.eq.13 .or. isw.eq.14 .or.
     &       isw.eq.16 .or. isw.eq.25) then

c       Integration order set to static

        if(d(7).ge.0.0 .and. (ndfo(1).gt.0 .or. shflg)) then
          cfac = d(7)
          lfac = 1.d0 - cfac
        else
          cfac = 0.0d0
          lfac = 0.0d0
        endif

c       Initialize momenta

        if(isw.eq.13) then
          do j = 1,nel
            do i = 1,3
              mom(i,j) = 0.0d0
            end do ! i
          end do ! j
        endif

c       Proportional body forces

        call sbodyf(d, body)

        augfp  = augf
        estore = 0.0d0

c       Set element quadrature order

        if(nel.eq.4) then
          if(d(182).gt.0.0d0) then
            call tint3dn(nel,lint,el)
          else
            l = 1
            call tint3d(l,lint,el)
          endif
          ord = 1
          npm = 1
          quad = .false.
        elseif(nel.eq.10 .or. nel.eq.11) then
          if(d(182).gt.0.0d0) then
            call tint3dn(nel,lint,el)
          else
            l = 4
            call tint3d(l,lint,el)
          endif
          if(nel.eq.10) then
            npm = 1
          else
            npm = 4
          endif
          ord = 2
          quad = .false.
        else
          if(nel.le.8) then
            npm = 1
            l   = 2
          else
            npm = 4
            l   = 3
          endif
          if(nint(d(182)).gt.0) then
            call int3dn(nel,lint,sg)
          else
            call int3d(l,lint,sg)
          endif
          quad = .true.
          ord = 0
        endif

c       Set number of history terms / quadradure point

        nhv    = nint(d(15))
        istrt  = nint(d(84))

c       Center estimate

        if(npm.gt.1) then
          do i = 1,2
            x0(i) = 0.0d0
            do j = 1,nel
              x0(i) = x0(i) + xl(i,j)
            end do ! j
            x0(i) = x0(i)/dble(nel)
          end do ! i
        endif

c       MECHANICAL ELEMENT

        do l = 1,lint

c         Shape functions and derivatives

          if(ord.gt.0) then
            call tetshp(el(1,l),xl,ndm,ord,xsj,shp(1,1,l))
            dvol(l) = xsj*el(5,l)
          else
            call shp3d(sg(1,l),xsj,shp(1,1,l),xl,ndm,nel)
            dvol(l) = xsj*sg(4,l)
          endif

c         Mixed volume effect and temperature projection

          do i = 1,3
            theta(i,l) = 0.0d0
          end do ! i

          ta(l) = -d(9)

c         Compute coordinates

          xx(1,l) = 0.0d0
          xx(2,l) = 0.0d0
          xx(3,l) = 0.0d0
          do i = 1,nel
            xx(1,l) = xx(1,l) + shp(4,i,l)*xl(1,i)
            xx(2,l) = xx(2,l) + shp(4,i,l)*xl(2,i)
            xx(3,l) = xx(3,l) + shp(4,i,l)*xl(3,i)
          end do ! i

c         Compute volumetric strain from displacements

          do i = 1,nel
            theta(1,l) = theta(1,l) + shp(1,i,l) * ul(1,i,1)
     &                              + shp(2,i,l) * ul(2,i,1)
     &                              + shp(3,i,l) * ul(3,i,1)
            theta(2,l) = theta(2,l) + shp(1,i,l) * ul(1,i,2)
     &                              + shp(2,i,l) * ul(2,i,2)
     &                              + shp(3,i,l) * ul(3,i,2)
            theta(3,l) = theta(3,l) + shp(1,i,l) * ul(1,i,3)
     &                              + shp(2,i,l) * ul(2,i,3)
     &                              + shp(3,i,l) * ul(3,i,3)
            ta(l)      = ta(l)      + shp(4,i,l) * tl(i)
          end do ! i

c         Set the pressure functions

          phi(1,l) = 1.d0
          if(npm.gt.1) then
            phi(2,l) = xx(1,l) - x0(1)
            phi(3,l) = xx(2,l) - x0(2)
            phi(4,l) = xx(3,l) - x0(3)
          endif
        end do ! l

c       Mixed model for volumetric and temperature response

        call bbar3s(phi,shp,dvol,lint,nel,npm,hh,theta,bbar)

c       Compute strains and stresses at quadrature points

        nn = nhi
        do l = 1,lint
          call strn3m(shp(1,1,l),xl,ul,theta(1,l),
     &                ndm,ndf,nel,nen,eps(1,1,l))

          epsv(l) = theta(1,l)

c         Set rotation angle

          if(nint(d(85)).eq.1) then
            psil = atan2(xx(2,l)-d(87),xx(1,l)-d(86))
          else
            psil = 0.0d0
          endif

          call modlsd(d,ta(l),eps(1,1,l),hr(nn+nh1),hr(nn+nh2),nhv,
     &                istrt,aa(1,1,1,l),sigl(1,l),isw)

c         Volumetric stress

          pbar(l) = one3*(sigl(1,l) + sigl(2,l) + sigl(3,l))

c         Compute energy

          if(isw.eq.13) then
            aj0 = d(14)*dvol(l)
            aj1 = d( 4)*aj0
            call sengy(ul,shp, aj0,aj1, lfac,cfac, 3,4, mom)
          endif

          nn = nn + nhv
        end do ! l

c       Accumulate energy

        if(isw.eq.13) then
          do j = 1,nel
            do i = 1,3
              epl(i) = epl(i) + mom(i,j)
            end do ! i
            epl(4) = epl(4) + xl(2,j)*mom(3,j) - xl(3,j)*mom(2,j)
            epl(5) = epl(5) + xl(3,j)*mom(1,j) - xl(1,j)*mom(3,j)
            epl(6) = epl(6) + xl(1,j)*mom(2,j) - xl(2,j)*mom(1,j)
          end do ! j
          return
  	endif

c       Integrate constant pressure over volume

        if(npm.eq.1) then

          mpress = 0.0d0
          do l = 1,lint
            mpress  = mpress  + pbar(l)*dvol(l)
          end do ! l

c         Divide pressure by volume

          press(1) = mpress * hh(1,1)
          do l = 2,lint
            press(l) = press(1)
          end do ! l

c       Higher order element pressures

        else

          do i = 1,npm
            sigm(i) = 0.0d0
          end do ! i
          do l = 1,lint
            mpress  = pbar(l)*dvol(l)
            sigm(1) = sigm(1) + mpress
            do i = 2,npm
              sigm(i) = sigm(i) + mpress*phi(i,l)
            end do ! i
          end do ! l

c         Divide pressure by reference volume

          do i = 1,npm
            hsig(i) = 0.0d0
            do j = 1,npm
              hsig(i) = hsig(i) + hh(i,j)*sigm(j)
            end do ! j
          end do ! i

          do l = 1,lint
            press(l) = hsig(1)
            do i = 2,npm
              press(l) = press(l) + hsig(i)*phi(i,l)
            end do ! i
          end do ! l

        endif

c       Compute mixed stress

        do l = 1,lint
          dsigtr    =  press(l)  - pbar(l)
          sigl(1,l) =  sigl(1,l) + dsigtr
          sigl(2,l) =  sigl(2,l) + dsigtr
          sigl(3,l) =  sigl(3,l) + dsigtr
        end do ! l

c       Tangent and residual computations

        if(isw.eq.3 .or. isw.eq.6 .or. isw.eq.14) then

c         Compute mixed pressure

          if(isw.eq.3 .or. isw.eq.6) then

            do l = 1,lint

c             Angular velocity body force: d(4) = density, d(65) = omega

              if(d(4).gt.0.0d0 .and. d(65).gt.0.0d0) then
                do i = 1,3
                  vl(i) = 0.0d0
                  bf(i) = 0.0d0
                  do j = 1,nel
                    vl(i) = vl(i) + shp(4,j,l)*xl(i,j)
                  end do ! j
                end do ! i
                call sbodyw(d(4),d(65),vl, bf, btan, .false.)
              else
                do i = 1,3
                  bf(i) = 0.0d0
                end do ! i
              endif

c             Store time history plot data for element

              i = 6*(l-1)
              do j = 1,6
                tt(j+i) = sigl(j,l)
                sigm(j) = sigl(j,l)*dvol(l)
              end do ! j

c             Compute acceleration

              do i = 1,3
                al(i) = 0.0d0
                do j = 1,nel
                  al(i) = al(i) + shp(4,j,l)*ul(i,j,5)
                end do ! j
                al(i) = al(i)*cfac
              end do ! i

c             Rayleigh damping

              if(d(77).ne.0.0d0) then
                do i = 1,3
                  vl(i) = 0.0d0
                  do j = 1,nel
                    vl(i) = vl(i) + shp(4,j,l)*ul(i,j,4)
                  end do ! j
                  vl(i)   = vl(i)*cfac
                end do ! i

c               Compute mass damping residual

                do i = 1,nel
                  fac    = shp(4,i,l)*dvol(l)*d(77)*d(4)
                  do j = 1,3
                    r(j,i) = r(j,i) - (vl(j) + lfac*ul(j,i,4))*fac
                  end do ! j
                end do ! i
              endif

              if(d(78).ne.0.0d0) then
                do i = 1,6
                  epsm(i) = 0.0d0
                end do ! i
                do j = 1,nel
                  epsm(1) = epsm(1) + shp(1,j,l)*ul(1,j,4)
                  epsm(2) = epsm(2) + shp(2,j,l)*ul(2,j,4)
                  epsm(3) = epsm(3) + shp(3,j,l)*ul(1,j,4)
                  epsm(4) = epsm(4) + shp(2,j,l)*ul(1,j,4)
     &                              + shp(1,j,l)*ul(2,j,4)
                  epsm(5) = epsm(5) + shp(3,j,l)*ul(2,j,4)
     &                              + shp(2,j,l)*ul(3,j,4)
                  epsm(4) = epsm(4) + shp(1,j,l)*ul(3,j,4)
     &                              + shp(3,j,l)*ul(1,j,4)
                end do ! j
                fac = one3*(theta(2,l) - epsm(1) - epsm(2) - epsm(3))
                epsm(1) = (epsm(1) + fac)*d(78)*dvol(l)
                epsm(2) = (epsm(2) + fac)*d(78)*dvol(l)
                epsm(3) = (epsm(3) + fac)*d(78)*dvol(l)
                do j = 4,6
                  epsm(j) =  epsm(j)*d(78)*dvol(l)
                end do ! j

                do i = 1,6
                  do j = 1,6
                    sigm(j) = sigm(j) + aa(j,i,5,l)*epsm(i)
                  end do ! j
                end do ! i
              endif

c             Compute residual

              do j = 1,nel

                do i = 1,3
                  ac(i)  = d(4)*(al(i) + lfac*ul(i,j,5))
                end do ! i

                dvs    = shp(4,j,l)*dvol(l)
                r(1,j) = r(1,j) + (body(1) + bf(1) - ac(1))*dvs
     &                          - shp(1,j,l)*sigm(1)
     &                          - shp(2,j,l)*sigm(4)
     &                          - shp(3,j,l)*sigm(6)
                r(2,j) = r(2,j) + (body(2) + bf(2) - ac(2))*dvs
     &                          - shp(1,j,l)*sigm(4)
     &                          - shp(2,j,l)*sigm(2)
     &                          - shp(3,j,l)*sigm(5)
                r(3,j) = r(3,j) + (body(3) + bf(3) - ac(3))*dvs
     &                          - shp(1,j,l)*sigm(6)
     &                          - shp(2,j,l)*sigm(5)
     &                          - shp(3,j,l)*sigm(3)
              end do ! j

c             Compute mixed tangent stiffness matrix

              if(isw.eq.3) then

c               Multiply tangent moduli by volume element

                call dmatmx( aa(1,1,1,l), dd )
                d1 = dvol(l)*ctan(1)
                do i = 1,7
                  do j = 1,7
                    dd(i,j) = dd(i,j)*d1
                  end do ! j
c                 Thermo-mechanical coupling
                  dther(i) = aa(i,1,2,l)*d1
                end do ! i

                if(d(78).ne.0.0d0) then
                  d1 = dvol(l)*d(78)*ctan(2)
                  do i = 1,6
                    do j = 1,6
                      dd(i,j) = dd(i,j) + aa(i,j,5,l)*d1
                    end do ! j
                  end do ! i
                endif

c               Mass factors

                dv = (ctan(3) + d(77)*ctan(2))*dvol(l)*d(4)*cfac
                dl = (ctan(3) + d(77)*ctan(2))*dvol(l)*d(4)*lfac

c               Compute row terms

                i1 = 0
                do i = 1,nel

c                 Compute bmat-t * dd * dvol

                  do jj = 1,7

                    bbd(1,jj) =  shp(1,i,l)*dd(1,jj)
     &                        +  shp(2,i,l)*dd(4,jj)
     &                        +  shp(3,i,l)*dd(6,jj)
     &                        + bbar(1,i,l)*dd(7,jj)

                    bbd(2,jj) =  shp(2,i,l)*dd(2,jj)
     &                        +  shp(1,i,l)*dd(4,jj)
     &                        +  shp(3,i,l)*dd(5,jj)
     &                        + bbar(2,i,l)*dd(7,jj)

                    bbd(3,jj) =  shp(3,i,l)*dd(3,jj)
     &                        +  shp(2,i,l)*dd(5,jj)
     &                        +  shp(1,i,l)*dd(6,jj)
     &                        + bbar(3,i,l)*dd(7,jj)
                  end do ! jj
c                          _
c                 Compute: B_trans * D * alpha * j * w

                  if(ther) then
                    fac = (dther(1) + dther(2) + dther(3))*one3
                    bt(1) =  shp(1,i,l)*dther(1) + shp(2,i,l)*dther(4)
     &                    +  shp(3,i,l)*dther(6)
     &                    + (bbar(1,i,l) - shp(1,i,l))*fac

                    bt(2) =  shp(2,i,l)*dther(2) + shp(1,i,l)*dther(4)
     &                    +  shp(3,i,l)*dther(5)
     &                    + (bbar(2,i,l) - shp(2,i,l))*fac

                    bt(3) =  shp(3,i,l)*dther(3) + shp(2,i,l)*dther(5)
     &                    +  shp(1,i,l)*dther(6)
     &                    + (bbar(3,i,l) - shp(3,i,l))*fac
                  endif

c                 Compute tangent stiffness

                  fac = shp(4,i,l)*dl
                  do jj = 1,3
                    s(i1+jj,i1+jj) = s(i1+jj,i1+jj) + fac
                  end do ! jj

                  dmshp = shp(4,i,l)*dv

                  j1 = 0
                  do j = 1,nel

c                   Inertial tangent

                    do jj = 1,3
                      s(i1+jj,j1+jj) = s(i1+jj,j1+jj) + dmshp*shp(4,j,l)
                    end do ! jj

c                   Compute mechanics part of tangent stiffness

                    do jj = 1,3
                      s(i1+jj,j1+1) = s(i1+jj,j1+1)
     &                              + bbd(jj,1)*shp(1,j,l)
     &                              + bbd(jj,4)*shp(2,j,l)
     &                              + bbd(jj,6)*shp(3,j,l)
     &                              + bbd(jj,7)*bbar(1,j,l)

                      s(i1+jj,j1+2) = s(i1+jj,j1+2)
     &                              + bbd(jj,2)*shp(2,j,l)
     &                              + bbd(jj,4)*shp(1,j,l)
     &                              + bbd(jj,5)*shp(3,j,l)
     &                              + bbd(jj,7)*bbar(2,j,l)

                      s(i1+jj,j1+3) = s(i1+jj,j1+3)
     &                              + bbd(jj,3)*shp(3,j,l)
     &                              + bbd(jj,5)*shp(2,j,l)
     &                              + bbd(jj,6)*shp(1,j,l)
     &                              + bbd(jj,7)*bbar(3,j,l)
                    end do ! jj

c                   Thermo-mechanical coupling matrix

                    if(ther) then
                      do jj = 1,3
                      s(i1+jj,j1+3)  = s(i1+jj,j1+3) + bt(jj)*shp(4,j,l)
                      end do ! jj
                    endif

                    j1 = j1 + ndf
                  end do ! j
                  i1 = i1 + ndf
                end do ! i
              endif ! isw = 3
            end do ! l

          endif ! isw = 3 or 6

c       Output stresses.

        elseif(isw.eq.4 .or. isw.eq.8) then

c         Output stresses

          if (isw .eq. 4) then

            do l = 1,lint
              do i = 1,6
                sigm(i) = sigl(i,l)
                epsm(i) = eps (i,1,l)
              end do ! i

              call pstr3d(sigm,sigm(7))
              call pstr3d(epsm,epsm(7))

              mct = mct - 5
              if(mct.le.0) then
                write(iow,2001) o,head
                if(ior.lt.0) write(*,2001) o,head
                mct = 50
              endif
              write(iow,2002) n,ma,(xx(i,l),i=1,3),
     &                        (sigm(i),i=1,6),(epsm(i),i=1,6),
     &                        (sigm(i),i=7,9),(epsm(i),i=7,9)
              if(ior.lt.0) then
                write(*,2002) n,ma,(xx(i,l),i=1,3),
     &                        (sigm(i),i=1,6),(epsm(i),i=1,6),
     &                        (sigm(i),i=7,9),(epsm(i),i=7,9)
              endif
            end do ! l

c         Project stresses onto nodes

          else
            call slcn3d(sigl,shp,dvol, r,s, lint,nel,10)
          endif

c       Compute J-integrals and material forces

        elseif(isw.eq.16) then

          call pjint3d(ul,tl,shp,dvol,epsv,sigl,r,ndf,ndm,lint,10)

c       Compute Z-Z projections

c       elseif(isw.eq.25) then

c         call stcn3z(xl,sigl,shp,dvol,lint,ndm,nel,16)

        endif ! isw = 4 or 8

      endif ! isw = 3 or 4 or 6 or 8 or 14 or 16

c     Formats

2001  format(a1,20a4//5x,'Element Stresses'//'     Elmt Mat',
     &   '    1-coord    2-coord    3-coord'/12x,
     &   '  11-stress  22-stress  33-stress  12-stress',
     &   '  23-stress  31-stress'/12x,
     &   '  11-strain  22-strain  33-strain  12-strain',
     &   '  23-strain  31-strain'/12x,
     &   '   1-stress   2-stress   3-stress',
     &   '   1-strain   2-strain   3-strain')

2002  format(/i9,i4,0p,3f11.3,1p,3e11.3/(12x,1p,6e11.3))

      end
