c$Id: sld1d3.f,v 1.1 2006/11/20 20:32:56 rlt Exp $
      subroutine sld1d3(d,ul,xl,ix,tl,s,p,ndf,ndm,nst,isw)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose:  1-d enhanced strain element for FEAP

c      Output records:
c      Prints in element: sig-11, sig-22, sig-33
c                         eps-11, eps-22, eps-33
c                         N.B. These are by default principal values

c      Prints at nodes:   1=sig-11, 2=sig-22, 3=sig-33
c                         psig-1  , psig-2    (computed by FEAP)

c      History Variable Storage (relative to nh1 or nh2)

c      Start           Description             Variable  Length
c      hr(0)           Enhanced displacement      ui(*,1)   3
c      hr(3)           Stress history at point-1    -      nhv
c      hr(3+  nhv)     Stress history at point-2    -      nhv

c      Total number of words / element is 3 + 2*nhv
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'bdata.h'
      include  'cdata.h'
      include  'eldata.h'
      include  'elplot.h'
      include  'eltran.h'
      include  'hdata.h'
      include  'incshp.h'
      include  'iofile.h'
      include  'part0.h'
      include  'pmod2d.h'
      include  'rdata.h'
      include  'comblk.h'

      logical   noconv
      integer   ndm,ndf,nst,isw, ix(*)
      real*8    a11,a31,ta,tol1,tolu
      real*8    cfac,lfac,sfac,dmas,lms,xsj0 ,rr
      integer   i,j,l,lint,nhv,nu1,nu2,ni,nn,nenit,i1,j1, istrt
      real*8    d(*),xl(ndm,*),ul(ndf,nen,*),tl(*),s(nst,*),p(ndf,*)
      real*8    shp(2,4,4),sg(2,4),sig(9,4),eps(9,3),epsv(4)
      real*8    gg(3,4),hh(3,3),bb(3),hg(3,4),dui(3),dd(6,6,5,4)
      real*8    ss(8,8),shpr(4,4),xsj(4),sigl(4),aa(4,4,4)
      real*8    dvol(4),r0(4)

      save

      data      ni /3/, eps / 27*0.0d0 /

c     Recover enhanced modes (saved in last iteration)

      nhv   = nint(d(15))
      istrt = nint(d(84))
      nu1   = nh1 - 1
      nu2   = nh2 - 1
      do i = 1,3
        ui(i,1)   = hr(nu2+i)
        ui(i,2)   = 0.0d0
      end do ! i

c     Compute quadrature and weights

      if(nint(d(182)).eq.1) then
        lint = nel
        call int1dn(lint,sg)
      else
        lint = nint(d(5))
        call int1d(lint,sg)
      endif

c     Initialize history variables only

      if(isw.eq.14) then
        nn = ni
        do l = 1,lint
          call modlsd(d,ta,eps,hr(nn+nh1),hr(nn+nh2),nhv,istrt,
     &                dd(1,1,1,l),sig(1,l),isw)
          nn = nn + nhv
        end do ! l
        return
      endif ! isw.eq.14

c     Compute shape functions

      do l = 1,lint
        call shp2d(sg(1,l),xl,shp(1,1,l),xsj(l),ndm,nel,ix,.false.)

c       Axisymmetry

        if(stype.eq.3) then
          r0(l) = 0.0d0
          do j = 1,nel
            r0(l) = r0(l) + shp(2,j,l)*xl(1,j)
          end do ! j
          dvol(l) = xsj(l)*sg(2,l)*r0(l)
          do j = 1,nel
            shpr(j,l) = shp(2,j,l)/r0(l)
          end do ! j

c       Plane

        else
          dvol(l) = xsj(l)*sg(2,l)
          do j = 1,nel
            shpr(j,l) = 0.0d0
          end do ! j
        endif
      end do ! l

c     Compute enhanced modes by local iteration

      tolu   =  1.d-03*tol*rnmax/dble(numel)
      nenit  =  0
      noconv = .true.
      do while(noconv)

        do j = 1,5
          bb(j) = 0.0d0
          do i = 1,5
            hh(i,j) = 0.0d0
          end do ! i
        end do ! j

c       Set initial counter for history terms in stress/strain

        nn = ni
        do l = 1,lint

c         Compute enhanced strain shape functions

          call shpi1d(sg(1,l),xl,ndm)

c         Compute stress and strain at point

          call stra1d(d,xl,ul,tl,shp(1,1,l),ndf,ndm,nel,rr,ta,eps)

          call modlsd(d,ta,eps,hr(nn+nh1),hr(nn+nh2),nhv,istrt,
     &                dd(1,1,1,l),sig(1,l),isw)

c         Scale moduli and stresses

          do i = 1,3
            do j = 1,3
              aa(j,i,l) = dd(j,i,1,l)*dvol(l)*ctan(1)
            end do ! j
            sigl(i) = sig(i,l)*dvol(l)
          end do ! i

c         Store time history plot data for element

          i = 12*(l-1)
          do j = 1,3
            tt(j+i)   = sig(j,l)
            tt(j+i+6) = eps(j,1)
          end do ! j

c         Enhanced residual computations

          do j = 1,2
            bb(j) = bb(j) - sigl(1)*shpi(1,j)
          end do ! j

          if(stype.eq.3) then
            shpi(2,3) = shpi(2,3)/r0(l)
            bb(3)     = bb(3) - sigl(3)*shpi(2,3)
            do j = 1,2
              hh(3,j) = hh(3,j) + shpi(2,3)*aa(3,1,l)*shpi(1,j)
              hh(j,3) = hh(j,3) + shpi(1,j)*aa(1,3,l)*shpi(2,3)
            end do ! j
            hh(3,3) = hh(3,3) + shpi(2,3)*aa(3,3,l)*shpi(2,3)
          else
            shpi(2,3) = 0.0d0
            hh(3,3)   = 1.0d0
          endif

c         Stiffness computations

          do j = 1,2
            a11 = aa(1,1,l)*shpi(1,j)
            do i = 1,2
              hh(i,j) = hh(i,j) + shpi(1,i)*a11
            end do ! i
          end do ! j
          nn = nn + nhv
        end do ! l

        call invert(hh,3,3)

c       Compute incremental enhanced displacements enhanced modes

        do i = 1,3
          dui(i)  = hh(i,1)*bb(1) + hh(i,2)*bb(2) + hh(i,3)*bb(3)
          ui(i,1) = ui(i,1) + dui(i)
          ui(i,2) = ui(i,2) + dui(i)
        end do ! i

c       Check convergence

        tol1 = abs(bb(1)*dui(1) + bb(2)*dui(2) + bb(3)*dui(3))

        if(tol1.le.tolu .and. nenit.ge.1) then
          noconv = .false.
        endif
        nenit = nenit + 1
        if(nenit.ge.3 .or. tolu.eq.0.0d0) then
          noconv = .false.
        endif
      end do ! while

c     Save enhanced modes

      do i = 1,3
        hr(nu2+i) = ui(i,1)
      end do ! i

c     Set initial counter for history terms in stress/strain

      if(isw.eq.3. or. isw.eq.6) then

c       Time integration order set to static or dynamic

        if(d(7).ge.0.0 .and. (ndfo(1).gt.0 .or. shflg)) then
          cfac = d(7)
          lfac = 1.0d0 - cfac
        else
          cfac = 0.0d0
          lfac = 0.0d0
        endif

        do i = 1,nel
          do j = 1,3
            gg(j,i) = 0.0d0
          end do ! j
        end do ! i

        nn = ni
        do l = 1,lint

          call shpi1d(sg(1,l),xl,ndm)

          if(stype.eq.3) then
            shpi(2,3) = shpi(2,3)/r0(l)
            xsj0      = xsj(l)*sg(2,l)
          else
            shpi(2,3) = 0.0d0
            xsj0      = 0.0d0
          endif

c         Rayleigh damping effects

          if(d(78).ne.0.0d0) then
            call rays1d(d,shp(1,1,l),sig(1,l),dd(1,1,5,l),ul(1,1,4),
     &                  xl,ndf,ndm,nel)
            sfac = d(78)*ctan(2)
          else
            sfac = 0.0d0
          endif

c         Residual computations

          call resid1d(cfac,lfac,dvol(l),xsj0,shp(1,1,l),sig(1,l),d,
     &                 ul(1,1,4),ul(1,1,5), p,ndf)

c         Stiffness computations

          if(isw.eq.3) then

            dmas = d(4)*(ctan(3) + d(77)*ctan(2))*dvol(l)

            do j = 1,3
              do i = 1,3
                aa(i,j,l) = aa(i,j,l) + dd(i,j,5,l)*sfac
              end do ! i
            end do ! j

            j1 = 1
            do j = 1,nel

c             Compute d * b matrix = a

              a11 = aa(1,1,l)*shp(1,j,l) + aa(1,3,l)*shpr(j,l)
              a31 = aa(3,1,l)*shp(1,j,l) + aa(3,3,l)*shpr(j,l)

c             Lumped mass effects

              lms          = shp(2,j,l)*dmas
              s(j1  ,j1  ) = s(j1  ,j1  ) + lms*lfac

              i1 = 1
              do i = 1,nel
                s(i1  ,j1  ) = s(i1  ,j1  ) + shp(1,i,l)*a11
     &                                      + shpr(i,l)*a31
     &                                      + lms*shp(2,i,l)*cfac
                i1 = i1 + ndf
              end do ! i
              do i = 1,2
                gg(i,j) = gg(i,j) + shpi(1,i)*a11
              end do ! i
              gg(3,j) = gg(3,j) + shpi(3,3)*a31
              j1      = j1 + ndf
            end do ! j
          endif
          nn = nn + nhv
        end do ! l

c       Eliminate enhanced modes

        do i = 1,3
          do j = 1,nel
            hg(i,j) = hh(1,i)*gg(1,j) + hh(2,i)*gg(2,j)
     &              + hh(3,i)*gg(3,j)
          end do ! j
        end do ! i

        if(isw.eq.3) then
          do j = 1,nel
            do i = 1,nel
              ss(i,j) = gg(1,i)*hg(1,j) + gg(2,i)*hg(2,j)
     &                + gg(3,i)*hg(3,j)
            end do ! i
          end do ! j

c         Construct static condensation

          j1 = 0
          do j = 1,nel
            i1 = 0
            do i = 1,nel
              s(i1,j1) = s(i1,j1) - ss(i,j)
              i1           = i1 + ndf
            end do ! i
            j1 = j1 + ndf
          end do ! j

c         Compute reduced residual

          do j = 1,nel
            p(1,j) = p(1,j) - hg(1,j)*bb(1) - hg(2,j)*bb(2)
     &                      - hg(3,j)*bb(3)
          end do ! j
        endif

c     Compute and output element variables

      elseif(isw.eq.4 .or. isw.eq.8 .or. isw.eq.16 .or. isw.eq.25) then

        if(nint(d(182)).eq.1) then
          lint = nel
          call int1dn(lint,sg)
        else
          lint = nint(d(5))
          call int1d(lint,sg)
        endif

c       Set initial counter for history terms in stress/strain

        nn = ni
        do l = 1,lint
          call shp1d(sg(1,l),xl,shp(1,1,l),ndm,nel,xsj(l))
          call shpi1d(sg(1,l),xl,ndm)

c         Compute stress and strain at point

          call stra1d(d,xl,ul,tl,shp(1,1,l),ndf,ndm,nel,rr,ta,eps)

          epsv(l) = eps(1,1) + eps(2,1) + eps(3,1)

          call modlsd(d,ta,eps,hr(nn+nh1),hr(nn+nh2),nhv,istrt,
     &                dd,sig(1,l),isw)

c         Compute principal stress values

          if(isw.eq.4) then
            mct = mct - 3
            if(mct.le.0) then
              write(iow,2001) o,head
              if(ior.lt.0) write(*,2001) o,head
              mct = 50
            endif
            write(iow,2002) n,ma,rr,(sig(i,l),i=1,3),
     &                              (eps(i,1),i=1,3)
            if(ior.lt.0) then
              write(*,2002) n,ma,rr,(sig(i,l),i=1,3),
     &                              (eps(i,1),i=1,3)
            endif
          endif
          nn = nn + nhv
        end do ! l

c       Plot stress values

        if(isw.eq.8) then

          call slcn1d(sig,shp,xsj,p,s,p(nen+1,1),lint,nel,9)

c       Compute J-integrals and material forces

        elseif(isw.eq.16) then

          call pjint1d(d,ul,tl,shp,xsj,epsv,sig,p,ndf,ndm,lint)

c       Compute Z-Z projections

        elseif(isw.eq.25) then

          call stcn1z(xl,sig,shp,xsj,lint,ndm,nel,9)

        endif

      endif

c     Formats

2001  format(a1,20a4//5x,'Element Stresses'//'    Elmt Mat',
     &     '    1-coord   11-stress   22-stress   33-stress'/
     &            12x,'   11-strain   22-strain   33-strain'/30(' -'))
2002  format(i8,i4,1p,4e12.4/1p,3e12.4)

      end
