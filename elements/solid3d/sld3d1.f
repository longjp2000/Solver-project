c$Id: sld3d1.f,v 1.1 2006/11/20 20:32:56 rlt Exp $
      subroutine sld3d1(d,ul,xl,tl,s,r,ndf,ndm,nst,isw)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: 3-D linear elastic displacment element for feap

c     Output records:

c     Prints in element: sig-11, sig-22, sig-33, sig-12, sig-23, sig-31
c                        eps-11, eps-22, eps-33, eps-12, eps-23, eps-31

c     Prints at nodes:   1=sig-11, 2=sig-22, 3=sig-33,
c                        4=sig-12  5=sig-23, 6=sig-31
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'bdata.h'
      include  'cdata.h'
      include  'complx.h'
      include  'eldata.h'
      include  'elplot.h'
      include  'eltran.h'
      include  'hdata.h'
      include  'iofile.h'
      include  'part0.h'
      include  'prstrs.h'
      include  'ptdat6.h'
      include  'rdata.h'
      include  'comblk.h'

      integer   i,j,l,nn,i1,j1,ndf,ndm,nst,isw,lint,nhv,istrt,ord,tdof
      real*8    shp(4,27,125),sg(4,125),sig(10,125),eps(9,3),dd(6,6,5)
      real*8    d(*),xl(ndm,*),ul(ndf,nen,*),tl(*),s(nst,*),r(*)
      real*8    th(125),sv(5,16),epsv(125),xsj(125),dv(125)
      real*8    mom(3,27)
      real*8    dvk, dvm, xn, yn, zn, ta, rr, lfac, cfac, am, mass
      real*8    a11, a12, a13, a21, a22, a23, a31, a32, a33
      real*8    a41, a42, a43, a51, a52, a53, a61, a62, a63

      save

c     Set nodal temperatures: Can be specified or computed

      if(isw.gt.1) then
        tdof = d(19)
        if(tdof.le.0) then
          do i = 1,nel ! {
            th(i) = tl(i)
          end do ! i     }
        else
          do i = 1,nel ! {
            th(i) = ul(tdof,i,1)
          end do ! i     }
        endif
      endif

c     Compute element tangent array

      nhv   = nint(d(15))
      istrt = nint(d(84))

c     Get quadrature information

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

c     Compute shape functions

      do l = 1,lint
        if(ord.gt.0) then
          call tetshp(sv(1,l),xl,ndm,ord,xsj(l),shp(1,1,l))
          dv(l) = xsj(l)*sv(5,l)
        else
          call shp3d(sg(1,l),xsj(l),shp(1,1,l),xl,ndm,nel)
          dv(l) = xsj(l)*sg(4,l)
        endif
      end do ! l

c     Compute the residual and tangent arrays or energy

      if(isw.eq.3 .or. isw.eq.6 .or. isw.eq.13) then

c       Set mass factors

        rr   = d(4)
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

c       Compute residual and tangent arrays

        nn = 0
        do l = 1,lint

c         Compute strain at point

          call strn3d(d,xl,ul,th,shp(1,1,l),3,ndf,nel, eps,ta)

c         Compute stress at point

          call modlsd(d,ta,eps,hr(nh1+nn),hr(nh2+nn),nhv,istrt,
     &                dd,sig(1,l),isw)

c         Residual and tangent computation

          if(isw.eq.3 .or. isw.eq.6) then

c           Residual computations

c           Store time history plot data for element

            i = 6*(l-1)
            do j = 1,6
              tt(j+i) = sig(j,l)
            end do ! j

c           Add stiffness part of Rayleigh damping to stress

            if(d(78).ne.0.0d0) then
              call rays3d(d,shp(1,1,l),shp(1,1,l),sig(1,l),dd,ul(1,1,4),
     &                    ndf,nel,.false.)
            endif

c           Form residual

            call resid3d(dv(l),shp(1,1,l),sig(1,l),d,xl,
     &                   ul(1,1,4),ul(1,1,5),r,3,ndf)

c           Stiffness computations

            if(isw.eq.3) then

              dvm   = rr*(ctan(3) + d(77)*ctan(2))*dv(l)
              dvk   =    (ctan(1) + d(78)*ctan(2))*dv(l)

              j1 = 1
              do j = 1,nel

c               Compute d * b matrix = a

                xn  = shp(1,j,l)*dvk
                yn  = shp(2,j,l)*dvk
                zn  = shp(3,j,l)*dvk
                a11 = dd(1,1,1)*xn + dd(1,4,1)*yn + dd(1,6,1)*zn
                a21 = dd(2,1,1)*xn + dd(2,4,1)*yn + dd(2,6,1)*zn
                a31 = dd(3,1,1)*xn + dd(3,4,1)*yn + dd(3,6,1)*zn
                a41 = dd(4,1,1)*xn + dd(4,4,1)*yn + dd(4,6,1)*zn
                a51 = dd(5,1,1)*xn + dd(5,4,1)*yn + dd(5,6,1)*zn
                a61 = dd(6,1,1)*xn + dd(6,4,1)*yn + dd(6,6,1)*zn
                a12 = dd(1,2,1)*yn + dd(1,4,1)*xn + dd(1,5,1)*zn
                a22 = dd(2,2,1)*yn + dd(2,4,1)*xn + dd(2,5,1)*zn
                a32 = dd(3,2,1)*yn + dd(3,4,1)*xn + dd(3,5,1)*zn
                a42 = dd(4,2,1)*yn + dd(4,4,1)*xn + dd(4,5,1)*zn
                a52 = dd(5,2,1)*yn + dd(5,4,1)*xn + dd(5,5,1)*zn
                a62 = dd(6,2,1)*yn + dd(6,4,1)*xn + dd(6,5,1)*zn
                a13 = dd(1,3,1)*zn + dd(1,5,1)*yn + dd(1,6,1)*xn
                a23 = dd(2,3,1)*zn + dd(2,5,1)*yn + dd(2,6,1)*xn
                a33 = dd(3,3,1)*zn + dd(3,5,1)*yn + dd(3,6,1)*xn
                a43 = dd(4,3,1)*zn + dd(4,5,1)*yn + dd(4,6,1)*xn
                a53 = dd(5,3,1)*zn + dd(5,5,1)*yn + dd(5,6,1)*xn
                a63 = dd(6,3,1)*zn + dd(6,5,1)*yn + dd(6,6,1)*xn

c               Add diagonal mass effects

                am           = shp(4,j,l)*dvm
                s(j1  ,j1  ) = s(j1  ,j1  ) + am*lfac
                s(j1+1,j1+1) = s(j1+1,j1+1) + am*lfac
                s(j1+2,j1+2) = s(j1+2,j1+2) + am*lfac

                i1 = 1
                do i = 1,nel

c                 Compute consistent mass matrix

                  xn   = shp(1,i,l)
                  yn   = shp(2,i,l)
                  zn   = shp(3,i,l)
                  mass = shp(4,i,l)*am*cfac

                  s(i1  ,j1  ) = s(i1  ,j1  ) + xn*a11 + yn*a41 + zn*a61
     &                                        + mass
                  s(i1  ,j1+1) = s(i1  ,j1+1) + xn*a12 + yn*a42 + zn*a62
                  s(i1  ,j1+2) = s(i1  ,j1+2) + xn*a13 + yn*a43 + zn*a63
                  s(i1+1,j1  ) = s(i1+1,j1  ) + yn*a21 + xn*a41 + zn*a51
                  s(i1+1,j1+1) = s(i1+1,j1+1) + yn*a22 + xn*a42 + zn*a52
     &                                        + mass
                  s(i1+1,j1+2) = s(i1+1,j1+2) + yn*a23 + xn*a43 + zn*a53
                  s(i1+2,j1  ) = s(i1+2,j1  ) + zn*a31 + yn*a51 + xn*a61
                  s(i1+2,j1+1) = s(i1+2,j1+1) + zn*a32 + yn*a52 + xn*a62
                  s(i1+2,j1+2) = s(i1+2,j1+2) + zn*a33 + yn*a53 + xn*a63
     &                                        + mass
                  i1 = i1 + ndf
                end do ! i
                j1 = j1 + ndf
              end do ! j
            endif
            if(cplxfl) then
              call cst3d1(shp(1,1,l),dv(l),eps,ul(1,1,8),
     &                    dd,dd(1,1,3),s(1,nst+1),r,r(nst+1))
            endif

c         Compute energy

          elseif(isw.eq.13) then

            dvm = rr*dv(l)
            call sengy(ul,shp(1,1,l), dv(l),dvm, lfac,cfac, 3,4, mom)

          endif ! isw
          nn = nn + nhv
        end do ! l

c       Accumulate momentum

        if(isw.eq.13) then
          do j = 1,nel
            do i = 1,3
              epl(i) = epl(i) + mom(i,j)
            end do ! i
            epl(4) = epl(4) + xl(2,j)*mom(3,j) - xl(3,j)*mom(2,j)
            epl(5) = epl(5) + xl(3,j)*mom(1,j) - xl(1,j)*mom(3,j)
            epl(6) = epl(6) + xl(1,j)*mom(2,j) - xl(2,j)*mom(1,j)
          end do ! j
        endif

c     Compute and output element variables

      elseif(isw.eq.4 .or. isw.eq.8 .or. isw.eq.16) then

c       Set initial counter for history terms in stress/strain

        nn = 0
        do l = 1,lint

c         Compute strain at point

          call strn3d(d,xl,ul,th,shp(1,1,l),ndm,ndf,nel, eps,ta)
          epsv(l) = eps(1,1) + eps(2,1) + eps(3,1)

c         Compute stress at point

          call modlsd(d,ta,eps,hr(nh1+nn),hr(nh2+nn),nhv,istrt,
     &                dd,sig(1,l),isw)

C         Output values

          if(isw.eq.4) then

c           Compute coordinates

            xn = 0.0
            yn = 0.0
            zn = 0.0
            do j = 1,nel
              xn = xn + shp(4,j,l)*xl(1,j)
              yn = yn + shp(4,j,l)*xl(2,j)
              zn = zn + shp(4,j,l)*xl(3,j)
            end do ! j

c           Compute principal stress values

            mct = mct - 3
            if(mct.le.0) then
              write(iow,2010) o,head
              if(ior.lt.0) write(*,2010) o,head
              mct = 50
            endif
            write(iow,2011) n,xn,(sig(i,l),i=1,6),ma,yn,(eps(i,1),i=1,6)
            if(ior.lt.0) then
              write(*,2011) n,xn,(sig(i,l),i=1,6),ma,yn,(eps(i,1),i=1,6)
            end if
          endif
          nn = nn + nhv
        end do ! l

c       Plot stress values

        if(isw.eq.8) then
          call slcn3d(sig,shp,xsj, r,s, lint,nel,27)

c       Compute J-integrals and material forces

        elseif(isw.eq.16) then
          call pjint3d(ul,tl,shp,xsj,epsv,sig,r,ndf,ndm,lint,27)
        endif
      endif

c     Formats

2010  format(a1,20a4//5x,'Element Stresses'//' Elmt 1-coord',
     1    2x,'11-stress  22-stress  33-stress  12-stress',
     2    2x,'23-stress  31-stress'/' matl 2-coord  11-strain',
     3    2x,'22-strain  33-strain  12-strain  23-strain',
     4    2x,'31-strain'/39(' -'))

2011  format(i4,0p1f9.3,1p6e11.3/i4,0p1f9.3,1p6e11.3/)

      end

      subroutine cst3d1(shp,xsj,epsr,ul,dr,di, si,pr,pi)

      implicit   none
      include   'cdata.h'
      include   'eldata.h'
      include   'eltran.h'
      include   'pmod2d.h'
      include   'sdata.h'

      integer    j,jj,j1,k,kk,k1
      real*8     shp(4,*),epsr(*),epsi(6),ul(ndf,*)
      real*8     dr(6,6),di(6,6), sigr(6),sigi(6),si(nst,*),pr(*),pi(*)
      real*8     xsj, aj1,aj2,aj3, bd(3,6)

c     Compute imaginary strains

      do j = 1,6
        epsi(j) = 0.0d0
      end do ! j
      do j = 1,nel
        epsi(1) = epsi(1) + shp(1,j)*ul(1,j)
        epsi(2) = epsi(2) + shp(2,j)*ul(2,j)
        epsi(3) = epsi(3) + shp(3,j)*ul(3,j)
        epsi(4) = epsi(4) + shp(1,j)*ul(2,j) + shp(2,j)*ul(1,j)
        epsi(5) = epsi(5) + shp(2,j)*ul(3,j) + shp(3,j)*ul(2,j)
        epsi(6) = epsi(6) + shp(3,j)*ul(1,j) + shp(1,j)*ul(3,j)
      end do ! j

c     Compute imaginary stress

      do j = 1,6
        sigr(j) = 0.0d0
        sigi(j) = 0.0d0
        do k = 1,6
          sigr(j) = sigr(j) + di(j,k)*epsi(k)
          sigi(j) = sigi(j) + dr(j,k)*epsi(k) + di(j,k)*epsr(k)
        end do ! k
      end do ! j

      j1 = 0
      do jj = 1,nel

        aj1 = shp(1,jj)*xsj
        aj2 = shp(2,jj)*xsj
        aj3 = shp(3,jj)*xsj

c       Compute B_trans * D * j * w

        do k = 1,6
          bd(1,k) = (aj1*di(1,k) + aj2*di(4,k) + aj3*di(6,k))*ctan(1)
          bd(2,k) = (aj2*di(2,k) + aj1*di(4,k) + aj3*di(5,k))*ctan(1)
          bd(3,k) = (aj3*di(3,k) + aj2*di(5,k) + aj1*di(6,k))*ctan(1)
        end do ! k

c       Loop over columns (symmetry noted)

        k1 = 0
        do kk = 1,nel
          do j = 1,3
            si(j1+j,k1+1) = si(j1+j,k1+1) + bd(j,1)*shp(1,kk)
     &                                    + bd(j,4)*shp(2,kk)
     &                                    + bd(j,6)*shp(3,kk)

            si(j1+j,k1+2) = si(j1+j,k1+2) + bd(j,4)*shp(1,kk)
     &                                    + bd(j,2)*shp(2,kk)
     &                                    + bd(j,5)*shp(3,kk)

            si(j1+j,k1+3) = si(j1+j,k1+3) + bd(j,6)*shp(1,kk)
     &                                    + bd(j,5)*shp(2,kk)
     &                                    + bd(j,3)*shp(3,kk)
          end do ! j
          k1 = k1 + ndf
        end do ! kk

c       Residual for added real part

        pr(j1+1) = pr(j1+1) + aj1*sigr(1) + aj2*sigr(4) + aj3*sigr(6)
        pr(j1+2) = pr(j1+2) + aj2*sigr(2) + aj1*sigr(4) + aj3*sigr(5)
        pr(j1+3) = pr(j1+3) + aj3*sigr(3) + aj2*sigr(5) + aj1*sigr(6)

c       Residual for imaginary part

        pi(j1+1) = pi(j1+1) - aj1*sigi(1) - aj2*sigi(4) - aj3*sigi(6)
        pi(j1+2) = pi(j1+2) - aj2*sigi(2) - aj1*sigi(4) - aj3*sigi(5)
        pi(j1+3) = pi(j1+3) - aj3*sigi(3) - aj2*sigi(5) - aj1*sigi(6)

        j1 = j1 + ndf
      end do ! jj

      end
