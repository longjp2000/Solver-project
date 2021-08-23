c$Id: fld3d1.f,v 1.1 2006/11/20 20:32:56 rlt Exp $
      subroutine fld3d1(d,ul,xl,s,p,ndf,ndm,nst,isw)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: 3-D finite deformation displacement element
c      Remark: This a completely standard mechanical element
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'bdata.h'
      include  'cdata.h'
      include  'cdamag.h'
      include  'debugs.h'
      include  'elcoor.h'
      include  'eldata.h'
      include  'elengy.h'
      include  'elplot.h'
      include  'eltran.h'
      include  'hdata.h'
      include  'iofile.h'
      include  'ndata.h'
      include  'part0.h'
      include  'plast3f.h'
      include  'pmod2d.h'
      include  'prstrs.h'
      include  'ptdat6.h'
      include  'rdata.h'
      include  'tdata.h'
      include  'comblk.h'

      logical   bflg
      integer   ndf,ndm,nst,isw, i,ii,i1, j,jj,j1,  k
      integer   l,lint, nhv,nn,nhi,istrt, ord
      real*8    bdb, rho,rhol,rhom,pote, xipa,wena, tb4, ta
      real*8    d(*),ul(ndf,nen,*),xl(ndm,*),s(nst,*),p(ndf,*)
      real*8    sg(4,125),sv(5,16),xsj(125),xx(3),bf(3),bt(3,3)
      real*8    f(9,2,125),finv(3,3,125),df(3,3,125),detfi(2,125)
      real*8    shp(4,27,125),cmom(3),be(6), body(3)
      real*8    aa(6,6,5),dd(6,6,125), xu(3,125)
      real*8    sigv(13),sigl(10,125), sigp(10,125)
      real*8    bbd(3,6),r1(3,125), tbi(3)
      real*8    dvol(125),dvol0(125),weng(125),xr(3,125),ur(3,125)

      save

c     MATERIAL DATA

      if(isw.eq.1) then

        ta = 0.0d0

c     AUGMENTED LAGRANGIAN UPDATE FOR NESTED ITERATION

      elseif(isw.eq.10) then

        if(hr(nh2).ne.0.0d0) then
          hr(nh2+1) = hr(nh2+1) + d(1)*log(abs(hr(nh2)))
        endif


c     COMPUTE TANGENT STIFFNESS/RESIDUAL FORCE AND OUTPUTS

      else

c       Get quadrature information

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

c       Compute shape functions

        do l = 1,lint
          if(ord.gt.0) then
            call tetshp(sv(1,l),xl,ndm,ord,xsj(l),shp(1,1,l))
            dvol0(l) = xsj(l)*sv(5,l)
          else
            call shp3d(sg(1,l),xsj(l),shp(1,1,l),xl,ndm,nel)
            dvol0(l) = xsj(l)*sg(4,l)
          endif
        end do ! l

c       Compute coordinates

        do l = 1,lint
          do i = 1,3
            xr(i,l) = 0.0d0
            ur(i,l) = 0.0d0
            do j = 1,nel
              xr(i,l) = xr(i,l) + shp(4,j,l)*xl(i,j)
              ur(i,l) = ur(i,l) + shp(4,j,l)*ul(i,j,1)
            end do ! j
            ur(i,l) = ur(i,l) + xr(i,l)
          end do ! i
        end do ! l

c       Degrees of freedom with enhanced modes.

        nhv   = nint(d(15))
        istrt = nint(d(84))

c       Initialize history variables

        if(isw.eq.14) then

          do l = 1,lint
            do i = 1,9
              f(i,1,l)    = 0.0d0
              f(i,2,l)    = 0.0d0
              finv(i,1,l) = 0.0d0
            end do ! i
            detfi(1,l)  = 1.d0
            detfi(2,l)  = 1.d0
            do i = 1,9,4
              f(i,1,l)    = 1.0d0
              f(i,2,l)    = 1.0d0
            end do ! i
            finv(1,1,l) = 1.0d0
            finv(2,2,l) = 1.0d0
            finv(3,3,l) = 1.0d0
          end do ! l

        else

c         Compute deformation gradient, inverse and determinant.

          do l = 1,lint
            call kine3df(shp(1,1,l),ul,f(1,1,l),finv(1,1,l),df(1,1,l),
     &                   detfi(1,l),ndf,nel,nen)
            dvol(l) = dvol0(l)*detfi(1,l)
          end do ! l

        endif

c       Compute Cauchy stresses and spatial tangent tensor at t-n+1

        nhi = 2
        nn  = nhi
        do l = 1,lint

          estore = 0.0d0
          do i = 1,3
            xref(i) = xr(i,l)
            xcur(i) = xr(i,l) + ur(i,l)
          end do ! i
          call modlfd(d,f(1,1,l),finv(1,1,l),df(1,1,l),detfi(1,l),ta,
     &                hr(nh1+nn),hr(nh2+nn),nhv,istrt,aa,sigv,be,
     &                hr(nh2),hr(nh2+1),.false.,isw)
          weng(l) = estore

c         Multiply tangent moduli and stresses by volume element.
c         Store time history plot data for element

          k = 6*(l-1)
          do i = 1,6
            tt(i+k)   = sigv(i)
            sigp(i,l) = sigv(i)
            sigl(i,l) = sigv(i)*dvol(l)
            do j = 1,6
              dd(j,i,l) = aa(j,i,1)*dvol(l)*ctan(1)
            end do ! j
          end do ! i
          nn = nn + nhv

        end do ! l

c       ENERGY COMPUTATIONS

        if(isw.eq.13) then

          rhom = 0.0d0
          rhol = 0.0d0
          pote = 0.0d0
          do l = 1,lint

            rho = d(4)*dvol0(l)
            do i = 1,3
              tbi(i) = 0.0d0
            end do ! i
            do j = 1,nel
              rhol = rhol + (ul(1,j,4)**2 +  ul(2,j,4)**2
     &                    +  ul(3,j,4)**2)*shp(4,j,l)*rho
              do i = 1,3
                tbi(i) = tbi(i) + shp(4,j,l)*ul(i,j,4)
              end do ! i
            end do ! j
            rhom = rhom + (tbi(1)**2 + tbi(2)**2 + tbi(3)**2)*rho
            pote = pote +  weng(l)*dvol0(l)
          end do ! l

c         Accumulate energies

          epl(7) = epl(7) + (rhol + d(7)*(rhom - rhol))*0.5d0
          epl(8) = epl(8) +  pote

c       STIFFNESS AND RESIDUAL

        elseif(isw.eq.3 .or. isw.eq.6) then

c         Compute body forces values

          call sbodyf(d, body)
          bflg = d(4).gt.0.0d0 .and. d(65).gt.0.0d0  ! angular velocity

          do l = 1,lint

c           Angular velocity: d(4) = rho; d(65) = omega

            do i = 1,3
              bf(i) = 0.0d0
            end do ! i
            if(bflg) then
              call sbodyw(d(4),d(65),ur(1,l), bf,bt, .true.)
              do ii = 1,3
                do jj = 1,3
                  bt(jj,ii) = bt(jj,ii)*dvol0(l)
                end do ! jj
              end do ! ii
            endif

c           Compute change in momentum

            if(d(7).ge.0.0d0) then
              rho  = d(4)*dvol0(l)
              rhom = rho * d(7)
              rhol = rho - rhom
            else
              rhom = 0.0d0
              rhol = 0.0d0
            endif
            do i = 1,3
              cmom(i) = 0.0d0
              do j = 1,nel
                cmom(i) = cmom(i) + shp(4,j,l)*ul(i,j,5)
              end do ! j
              cmom(i) = rhom*cmom(i)
            end do ! i

c           COMPUTE STRESS DIVERGENCE TERM

c           Compatible internal force.

            do i = 1,nel
              r1(1,i) = shp(1,i,l)*sigl(1,l)
     &                + shp(2,i,l)*sigl(4,l)
     &                + shp(3,i,l)*sigl(6,l)
              r1(2,i) = shp(1,i,l)*sigl(4,l)
     &                + shp(2,i,l)*sigl(2,l)
     &                + shp(3,i,l)*sigl(5,l)
              r1(3,i) = shp(1,i,l)*sigl(6,l)
     &                + shp(2,i,l)*sigl(5,l)
     &                + shp(3,i,l)*sigl(3,l)

              do j = 1,3
                p(j,i)  = p(j,i) + shp(4,i,l)*(body(j) + bf(j))*dvol0(l)
     &                  - shp(4,i,l)*(cmom(j)+rhol*ul(j,i,5)) - r1(j,i)
              end do ! j
            end do ! i

c           COMPUTE K11

            if(isw.eq.3) then

c             Inertia factor

              if(d(7).ge.0.0 .and. (ndfo(1).gt.0 .or. shflg)) then
                rhom  = rho*ctan(3)
              else
                rhom = 0.0d0
              endif

              i1 = 0
              do i = 1,nel

c               PART 1. - geometric and inertial parts.

                tb4   = shp(4,i,l)*rhom
                rhol  = tb4*(1.d0 - d(7))
                do jj = 1,3
                  s(i1+jj,i1+jj) = s(i1+jj,i1+jj) + rhol
                end do ! jj

                tb4 = tb4*d(7)
                j1 = 0
                if(gflag) then
                  do j = 1,nel

c                   Accumulate geometric factor with consistent mass

                    bdb = (r1(1,i)*shp(1,j,l)
     &                  +  r1(2,i)*shp(2,j,l)
     &                  +  r1(3,i)*shp(3,j,l))*ctan(1)
     &                  +  tb4*shp(4,j,l)

                    do jj = 1,3
                      s(i1+jj,j1+jj) = s(i1+jj,j1+jj) + bdb
                    end do ! jj
                    j1 = j1 + ndf
                  end do ! j
                else
                  do j = 1,nel
                    do jj = 1,3
                      s(i1+jj,j1+jj) = s(i1+jj,j1+jj)
     &                               + tb4*shp(4,j,l)
                    end do ! jj
                    j1 = j1 + ndf
                  end do ! j
                endif

c               Angular velocity tangent

                if(bflg) then
                  do jj = 1,3
                    do ii = 1,3
                      bdb = shp(4,i,l)*bt(ii,jj)
                      j1  = 0
                      do j = 1,nel
                        s(i1+ii,j1+jj) = s(i1+ii,j1+jj) + bdb*shp(4,j,l)
                        j1             = j1 + ndf
                      end do ! j
                    end do ! ii
                  end do ! jj
                endif

c               PART 2. - tangent modulus part (based upon aa-array)

                do jj = 1,6
                  bbd(1,jj) = shp(1,i,l)*dd(1,jj,l)
     &                      + shp(2,i,l)*dd(4,jj,l)
     &                      + shp(3,i,l)*dd(6,jj,l)
                  bbd(2,jj) = shp(1,i,l)*dd(4,jj,l)
     &                      + shp(2,i,l)*dd(2,jj,l)
     &                      + shp(3,i,l)*dd(5,jj,l)
                  bbd(3,jj) = shp(1,i,l)*dd(6,jj,l)
     &                      + shp(2,i,l)*dd(5,jj,l)
     &                      + shp(3,i,l)*dd(3,jj,l)
                end do ! jj

c               Compute tangent stiffness

                j1 = 0
                do j  = 1,nel
                  s(i1+1,j1+1) = s(i1+1,j1+1)
     &                          + bbd(1,1)*shp(1,j,l)
     &                          + bbd(1,4)*shp(2,j,l)
     &                          + bbd(1,6)*shp(3,j,l)
                  s(i1+1,j1+2) = s(i1+1,j1+2)
     &                          + bbd(1,4)*shp(1,j,l)
     &                          + bbd(1,2)*shp(2,j,l)
     &                          + bbd(1,5)*shp(3,j,l)
                  s(i1+1,j1+3) = s(i1+1,j1+3)
     &                          + bbd(1,6)*shp(1,j,l)
     &                          + bbd(1,5)*shp(2,j,l)
     &                          + bbd(1,3)*shp(3,j,l)
                  s(i1+2,j1+1) = s(i1+2,j1+1)
     &                          + bbd(2,1)*shp(1,j,l)
     &                          + bbd(2,4)*shp(2,j,l)
     &                          + bbd(2,6)*shp(3,j,l)
                  s(i1+2,j1+2) = s(i1+2,j1+2)
     &                          + bbd(2,4)*shp(1,j,l)
     &                          + bbd(2,2)*shp(2,j,l)
     &                          + bbd(2,5)*shp(3,j,l)
                  s(i1+2,j1+3) = s(i1+2,j1+3)
     &                          + bbd(2,6)*shp(1,j,l)
     &                          + bbd(2,5)*shp(2,j,l)
     &                          + bbd(2,3)*shp(3,j,l)
                  s(i1+3,j1+1) = s(i1+3,j1+1)
     &                          + bbd(3,1)*shp(1,j,l)
     &                          + bbd(3,4)*shp(2,j,l)
     &                          + bbd(3,6)*shp(3,j,l)
                  s(i1+3,j1+2) = s(i1+3,j1+2)
     &                          + bbd(3,4)*shp(1,j,l)
     &                          + bbd(3,2)*shp(2,j,l)
     &                          + bbd(3,5)*shp(3,j,l)
                  s(i1+3,j1+3) = s(i1+3,j1+3)
     &                          + bbd(3,6)*shp(1,j,l)
     &                          + bbd(3,5)*shp(2,j,l)
     &                          + bbd(3,3)*shp(3,j,l)
                  j1 = j1 + ndf
                end do ! j
                i1 = i1 + ndf
              end  do ! i

            endif

          end do ! l

c       OUTPUT STRESSES

        elseif(isw.eq.4 .or. isw.eq.8 .or. isw.eq.16) then

          xipa = 0.0d0
          wena = 0.0d0
          do i = 1,3
            xx(i) = 0.0d0
          end do ! i
          do i = 1,13
            sigv(i) = 0.0d0
          end do ! i

          do l = 1,lint

            do j = 1,3
              do i=1,nel
                xx(j)   = xx(j) + 0.125*shp(4,i,l)*xl(j,i)
              end do ! i
            end do ! j

c           Move stresses and jacobian for printing

            xipa = xipa + 0.125d0*xipr
            wena = wena + 0.125d0*wengy
            do i = 1,6
              sigv(i) = sigv(i) + 0.125d0*sigp(i,l)
            end do ! i

c           Store arrays for plotting.

            sigp(10,l) = 0.0d0
            sigp(9,l)  = 0.0d0

          end do ! l

c         OUTPUT STRESSES

          if (isw .eq. 4) then

            call pstr3d(sigv,sigv(7))

            mct = mct - 2
            if(mct.le.0) then
              write(iow,2001) o,head
              if(ior.lt.0) write(*,2001) o,head
              mct = 50
            endif

            write(iow,2002) n,(sigv(ii),ii=1,6),
     &                     ma,(sigv(ii),ii=7,9),xx
            if(ior.lt.0) then
              write(*,2002) n,(sigv(ii),ii=1,6),
     &                     ma,(sigv(ii),ii=7,9),xx
            end if

c         PROJECT STRESSES ONTO THE NODES FOR PLOTTING

          elseif(isw.eq.8) then

c           Compute current geometry

            do i = 1,ndm
              do j = 1,nel
                xu(i,j) = xl(i,j) + ul(i,j,1)
              end do ! j
            end do ! i

            call slcn3d(sigp,shp,xsj, p,s, lint,nel,27)

c         COMPUTE FRACTURE INDICES

          elseif(isw.eq.16) then

            call pfrac3f(f,detfi,sigp,weng, shp,dvol, p,
     &                   lint,ndf,ndm,3)
          endif
        endif
      endif

c     FORMAT STATEMENTS

2001  format(a1,20a4//5x,'Element Stresses'//
     &       '   Elem.   11-Stress   22-Stress   33-Stress   12-Stress',
     &   '   23-Stress   13-Stress'/
     &       '   Matl.    1-Stress    2-Stress    3-Stress',
     &   '     1-Coord     2-Coord     3-Coord ')

2002  format(i8,1p,6e12.4/i8,1p,6e12.4/1x)

      end
