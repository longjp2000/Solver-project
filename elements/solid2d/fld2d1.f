c$Id: fld2d1.f,v 1.1 2006/11/20 20:32:56 rlt Exp $
      subroutine fld2d1(d,ul,xl,ix,s,r,ndf,ndm,nst,isw)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose:  2-D Finite Deformation Elasticity Routine
c                Remark: This is a standard displacement model

c      Inputs:
c         d(*)      - Material set parameters
c         ul(ndf,*) - Nodal solution parameters for element
c         xl(ndm,*) - Nodal coordinates for element
c         ix(*)     - Element nodal connection list
c         ndf       - Number dof/node
c         ndm       - Spatial dimension of mesh
c         nst       - Dimension of element arrays
c         isw       - Switch to control action

c      Outputs:
c         s(nst,*)  - Element matrix
c         p(nst)    - Element vector
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

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
      include  'pmod2d.h'
      include  'ptdat6.h'
      include  'rdata.h'
      include  'comblk.h'

      logical   dynflg,quad
      integer   ndf,ndm,nst,isw,i,i1,is,j,jj,j1,js,l,lint,ni,nn,nhv
      integer   istrt
      integer   ix(*)
      real*8    bdb,bd3,dmas0
      real*8    cfac,lfac,xx1,xx2,xx3,yy1, tempi, xlamd, ha, ta
      real*8    d(*),ul(ndf,nen,*),xl(ndm,*),s(nst,*),body(3)
      real*8    r(ndf,*),r1(3,16),vl(3),weng(25)
      real*8    df(3,3,25),f(9,2,25),finv(3,3,25),detf(2,25),xsj(25)
      real*8    ds(6,6,5),dd(6,6),sigv(9),shp(3,16,25),sigl(16,25)
      real*8    sg(3,25),bbd(6,3),bb(6),shpr(16),dvol0(25),dvol(25)
      real*8    xr(2,25),ur(2,25),el(4,7)

      save

      data      lint / 0 /
      data      xlamd,ha /  2*0.0d0 /
      data      f    /450*0.0d0 /
      data      finv /225*0.0d0 /
      data      ta   /    0.0d0 /

c     Check process isw = 1

      if(isw.eq.1) then
        return
      endif

c     Set quadrature points and weights

      if(nel.eq.3) then
        if(d(182).gt.0.0d0) then
          call tint2dn(nel,lint,el)
        else
          l = 1
          call tint2d(l,lint,el)
        endif
        quad = .false.
      elseif(nel.eq.6 .or. nel.eq.7 ) then
        if(d(182).gt.0.0d0) then
          call tint2dn(nel,lint,el)
        else
          l = 7
          call tint2d(l,lint,el)
        endif
        quad = .false.
      else
        if(nint(d(182)).gt.0) then
          call int2dn(nel,lint,sg)
        else
          l = nint(d(5))
          call int2d(l,lint,sg)
        endif
        quad = .true.
      endif

c     COMPUTE TANGENT STIFFNESS AND RESIDUAL FORCE VECTOR

c     Compute shape functions and derivatives in reference configuration

      do l = 1,lint
        if(quad) then
          call shp2d(sg(1,l),xl,shp(1,1,l),xsj(l),ndm,nel,ix,.false.)
          dvol0(l) = xsj(l)*sg(3,l)
        else
          call trishp(el(1,l),xl,ndm,nel-4,xsj(l),shp(1,1,l))
          dvol0(l) = xsj(l)*el(4,l)
        endif

c       Compute coordinates at gauss points

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

c       Axisymmetric volume

        if(stype.eq.3 .or. stype.eq.8) then
          dvol0(l) = dvol0(l)*xr(1,l)
        endif
      end do ! l
      xref(3) = 0.0d0
      xcur(3) = 0.0d0

c     Compute deformation gradient and determinant; transform shape
c     functions to current configuration.

      call kine2d(shp,xl,ul,f,finv,df,detf,ndm,ndf,nel,nen,lint)

      nhv   = nint(d(15))
      istrt = nint(d(84))
      ni    = 0

c     Set the loop limits and consistent/lumped mass factor
      if(stype.eq.8) then
        is = 3
        js = 6
        cfac = 1.0d0 ! Axisymmetric + torsion must be consistent mass
      else
        is = 2
        js = 4
        cfac = d(7)
        vl(3) = 0.0d0
        if(stype.eq.1) then ! plane stress modifications to F, det F
          call kineps(f,finv,df,detf,hr(nh1),hr(nh2), lint)
          ni = lint
        endif
      endif
      lfac = 1.0d0 - cfac

c     Transfer for output related values

      if(isw.eq.4 .or. isw.eq.8 .or. isw.eq.16 .or. isw.eq.25 ) go to 4

c     Compute body forces

      call sbodyf(d, body)

c     LOOP OVER GAUSS POINTS

      dynflg = ctan(3).ne.0.0d0
      nn     = ni
      do l = 1,lint

c       Set reference coordinates

        xref(1) = xr(1,l)
        xref(2) = xr(2,l)
        xcur(1) = xr(1,l) + ur(1,l)
        xcur(2) = xr(2,l) + ur(2,l)

c       Check for axisymmetry

        if(stype.eq.3 .or. stype.eq.8) then
          do i = 1,nel
            shpr(i) = shp(3,i,l)/xcur(1)
          end do ! i
        else
          do i = 1,nel
            shpr(i) = 0.0d0
          end do ! i
        end if

c       Compute Cauchy stresses and spatial tangent tensor

        call modlfd(d,f(1,1,l),finv(1,1,l),df(1,1,l),detf(1,l),ta,
     &             hr(nn+nh1),hr(nn+nh2),nhv,istrt, ds,sigv,bb,
     &             xlamd,ha,.false.,isw)

c       Save plane stress thickness deformation gradient (minus unity)

        if(stype.eq.1) then
          hr(nh2-1+l) = f(9,1,l) - 1.0d0
        endif

c       Compute volume and mass factor

        dvol(l) = dvol0(l)*detf(1,l)
        if(d(7).ge.0.0 .and. (ndfo(1).gt.0 .or. shflg)) then
          dmas0 = dvol0(l)*d(4)
        else
          dmas0 = 0.0d0   ! No inertia effects
        endif

        if(isw.eq.13) then

          epl(8) = epl(8) + estore*dvol0(l)

c         Compute velocity at point

          do i = 1,is
            vl(i) = 0.0d0
            do j = 1,nel
              vl(i) = vl(i) + ul(i,j,4)*shp(3,j,l)
            end do ! j
          end do ! i

          tempi = 0.0d0
          if(stype.eq.8) then
            do i = 1,nel
              tempi = tempi
     &          + (ul(1,i,4)**2+ul(2,i,4)**2+ul(3,i,4)**2)*shp(3,i,l)
            end do ! i
          else
            do i = 1,nel
              tempi = tempi
     &              + (ul(1,i,4)**2 + ul(2,i,4)**2)*shp(3,i,l)
            end do ! i
          endif

c         Accumulate kinetic energy

          epl(7) = epl(7) + 0.5d0*(lfac*tempi
     &                    + cfac*(vl(1)**2 + vl(2)**2 + vl(3)**2))*dmas0

        elseif(isw.ne.14) then

c         Store stress values for tplot

          j1 = 6*(l-1)
          do j = 1,js
            tt(j+j1) = sigv(j)
          end do ! j

c         Multiply tangent moduli and stresses by volume element.

          do i = 1,js
            sigv(i) = sigv(i)*dvol(l)
            do j = 1,js
              dd(i,j) = ds(i,j,1)*dvol(l)*ctan(1)
            end do ! j
          end do ! i

c         COMPUTE STRESS DIVERGENCE AND INERTIA TERMS

          do i = 1,nel

c           Stress divergence term (used in geometric stiffness)

            r1(1,i) = shp(1,i,l)*sigv(1) + shp(2,i,l)*sigv(4)
            r1(2,i) = shp(1,i,l)*sigv(4) + shp(2,i,l)*sigv(2)

c           Element residual

            r(1,i) = r(1,i) - r1(1,i) - shpr(i)*sigv(3)
     &                      + body(1)*dvol0(l)*shp(3,i,l)
            r(2,i) = r(2,i) - r1(2,i)
     &                      + body(2)*dvol0(l)*shp(3,i,l)
          end do ! i

c         Torsion residual

          if(stype.eq.8) then
            do i = 1,nel
              r1(3,i) = xcur(1)*(shp(1,i,l)*sigv(6)+shp(2,i,l)*sigv(5))
              r(3,i)  = r(3,i) - r1(3,i) + body(3)*shp(3,i,l)*dvol0(l)
              r1(3,i) = r1(3,i)*2.d0 ! Geometric term only
            end do ! i
          endif

c         COMPUTE K (s(nst,nst) = K)

          if(isw.eq.3) then

c           PART 1. - Geometric part.

            if(gflag) then
              i1  = 0
              do i = 1,nel
                bd3 = shpr(i)*sigv(3)*ctan(1)
                j1  = 0
                do j = 1,nel
                  bdb          = (r1(1,i)*shp(1,j,l)
     &                         +  r1(2,i)*shp(2,j,l))*ctan(1)
                  s(i1+1,j1+1) = s(i1+1,j1+1) + bdb + bd3*shpr(j)
                  s(i1+2,j1+2) = s(i1+2,j1+2) + bdb
                  if(stype.eq.8) then
                    s(i1+1,j1+3) = s(i1+1,j1+3)+shpr(i)*r1(3,j)*ctan(1)
                    s(i1+3,j1+1) = s(i1+3,j1+1)+shpr(j)*r1(3,i)*ctan(1)
                    s(i1+3,j1+3) = s(i1+3,j1+3)+xcur(1)*xcur(1)*bdb
                  endif
                  j1 = j1 + ndf
                end do ! j
                i1 = i1 + ndf
              end do ! i
            endif ! gflag

c           PART 2. - Tangent modulus part (based upon dd-array)

            i1 = 0
            do i  = 1,nel

c             Compute bmat-t * dd * dvol

              do jj = 1,js
                bbd(jj,1) = shp(1,i,l)*dd(1,jj)
     &                    + shpr( i  )*dd(3,jj)
     &                    + shp(2,i,l)*dd(4,jj)

                bbd(jj,2) = shp(1,i,l)*dd(4,jj)
     &                    + shp(2,i,l)*dd(2,jj)
              end do ! jj

              if(stype.eq.8) then
                do jj = 1,js
                  bbd(jj,3) = xcur(1)*(shp(2,i,l)*dd(5,jj)
     &                               + shp(1,i,l)*dd(6,jj))
                end do ! jj
              endif

c             Compute tangent stiffness

              j1 = 0
              do j  = 1,nel

                do jj = 1,is
                  s(i1+jj,j1+1) = s(i1+jj,j1+1) + bbd(1,jj)*shp(1,j,l)
     &                                          + bbd(3,jj)*shpr( j  )
     &                                          + bbd(4,jj)*shp(2,j,l)

                  s(i1+jj,j1+2) = s(i1+jj,j1+2) + bbd(4,jj)*shp(1,j,l)
     &                                          + bbd(2,jj)*shp(2,j,l)
                end do ! jj

c               Torsion part

                if(stype.eq.8) then
                  do jj = 1,3
                    s(i1+jj,j1+3) = s(i1+jj,j1+3)
     &                            + (bbd(5,jj)*shp(2,j,l)
     &                            +  bbd(6,jj)*shp(1,j,l))*xcur(1)
                  end do ! jj

                endif

                j1 = j1 + ndf
              end do ! j

              i1 = i1 + ndf
            end  do ! i

          endif ! end of tangent

c         Add inertia parts

          if(dynflg) then
            call fdyn2d(ul,shp(1,1,l),s,r,is,xcur(1),
     &                  cfac,lfac,dmas0,isw)
          endif

        endif ! end of isw options

        nn = nn + nhv

      end do ! l

c     Multiply by thickness if not unity

      if((isw.eq.3 .or. isw.eq.6) .and. d(14).ne.1.d0) then

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

      return

c     OUTPUT STRESSES

   4  xx1  = 0.d0
      xx2  = 0.d0
      xx3  = 0.d0
      do i = 1,6
        sigv(i) = 0.0d0
      end do ! i

c     LOOP OVER GAUSS POINTS

      nn = ni

      do l = 1,lint

c       Compute Cauchy stresses and spatial tangent tensor at t-n+1

        call modlfd(d,f(1,1,l),finv(1,1,l),df(1,1,l),detf(1,l),ta,
     &             hr(nn+nh1),hr(nn+nh2),nhv,istrt,ds,sigl(1,l),bb,
     &             xlamd,ha,.false.,isw)
        weng(l) = estore
        dvol(l) = dvol0(l)*detf(1,l)

        yy1 = 1.d0/dble(nel)
        do i=1,nel
          tempi = yy1 * shp(3,i,l)
          xx1   = xx1 + tempi*xl(1,i)
          xx2   = xx2 + tempi*xl(2,i)
        end do ! i

c       Compute average stresses and jacobian for printing

        do i = 1,js
          sigv(i) = sigv(i) + yy1 * sigl(i,l)
        end do ! i

        nn = nn + nhv

      end do ! l

c     Output stresses

      if(isw.eq.4) then
        mct = mct - 2
        if(mct.le.0) then
          write(iow,2001) o,head
          if(ior.lt.0) write(*,2001) o,head
          mct = 50
        endif

        write(iow,2002) n,ma,(sigv(jj),jj=1,6),xx1,xx2,xx3
        if(ior.lt.0) then
          write(*,2002) n,ma,(sigv(jj),jj=1,6),xx1,xx2,xx3
        end if

c     Project stress values to nodes

      elseif(isw.eq.8) then

        call slcn2d(ix,sigl,shp,xsj,r,s,r(nen+1,1),lint,nel,16)

c     Compute fracture indices

      elseif(isw.eq.16) then

        call pfrac2f(f,detf,sigl,weng, shp,dvol, r, lint,ndf,ndm,3)

c     Compute Z-Z projections

      elseif(isw.eq.25) then

        call stcn2z(xl,sigl,shp,dvol,lint,ndm,nel,16)

      end if

c     Format statements

2001  format(a1,20a4//5x,'Element Stresses'//'   Elmt Matl',
     1   '  11-stress  22-stress  33-stress  12-stress',
     2   '  23-stress  13-stress'/16x,'1-coord    2-coord    3-coord ')

2002  format(i8,i4,1p,6e11.3/12x,1p,6e11.3/12x,0p,3f11.5)

      end
