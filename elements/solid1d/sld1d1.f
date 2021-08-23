c$Id: sld1d1.f,v 1.1 2006/11/20 20:32:56 rlt Exp $
      subroutine sld1d1(d,ul,xl,tl,s,p,ndf,ndm,nst,isw, ther)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Plane and axisymmetric linear elastic element routine

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
c         p(ndf,*)  - Element vector
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'bdata.h'
      include  'cdata.h'
      include  'complx.h'
      include  'elbody.h'
      include  'eldata.h'
      include  'elplot.h'
      include  'eltran.h'
      include  'fdata.h'
      include  'hdata.h'
      include  'iofile.h'
      include  'part0.h'
      include  'pmod2d.h'
      include  'strnum.h'
      include  'rdata.h'
      include  'comblk.h'

      logical   ther
      integer   ndf,ndm,nst,isw,j,k,l,ii,j1,k1,lint,nhv,nn,istrt
      real*8    aj1,aj2,aj0,xx,lfac,cfac,sfac,xsj0,dv,ta
      real*8    bd11,bd13, bt1
      real*8    d(*),ul(ndf,nen,*),xl(ndm,*),tl(*),s(nst,*),p(ndf,*)
      real*8    shp(2,4,4),sg(2,4),sig(10,4),eps(9,3),epsv(4)
      real*8    xsj(4),dd(6,6,5),shpr(4)

      save

      data      eps/ 27*0.0d0 /

c     Compute stress-divergence vector (p) and stiffness matrix (s)

      nhv   = nint(d(15))
      istrt = nint(d(84))

      if(isw.eq.3  .or. isw.eq.6  .or. isw.eq.14) then

c       Integration order set to static

        if(d(7).ge.0.0 .and. (ndfo(1).gt.0 .or. shflg)) then
          cfac = d(7)
          lfac = 1.d0 - cfac
        else
          cfac = 0.0d0
          lfac = 0.0d0
        endif

c       Compute Gauss quadrature points and weights

        if(nint(d(182)).eq.1) then
          lint = nel
          call int1dn(lint,sg)
        else
          lint = nint(d(5))
          call int1d(lint,sg)
        endif

c       Zero shpr matrix

        do j = 1,3
          shpr(j) = 0.0d0
        end do ! j

c       Numerical integration loop

        nn = 0
        do l = 1,lint
          call shp1d(sg(1,l),xl,shp(1,1,l),ndm,nel,xsj(l))
          xsj(l) = xsj(l)*sg(2,l)

c         Compute stresses and strains

          call stra1d(d,xl,ul,tl,shp,ndf,ndm,nel, xx,ta,eps)
          call modlsd(d,ta,eps,hr(nh1+nn),hr(nh2+nn),nhv,istrt,
     &                dd,sig(1,l),isw)

          if(isw.eq.3 .or. isw.eq.6) then

c           Multiply jacobian by radius for axisymmetry

            if(stype.eq.3) then
              xsj0   = xsj(l)
              xsj(l) = xsj(l)*xx
              do j = 1,nel
                shpr(j) = shp(2,j,1)/xx
              end do ! j
            else
              xsj0 = 0.0d0
            end if

c           Store time history plot data for element

            k = 10*(l-1)
            do j = 1,3
              tt(j+k) = sig(j,l)
            end do ! j
            k = k + 6
            do j = 1,3
              tt(j+k) = eps(j,1)
            end do ! j

c           Rayleigh damping effects

            dv = d(4)*(ctan(3) + d(77)*ctan(2))*xsj(l)

            if(d(78).ne.0.0d0) then
              call rays1d(d,shp,sig(1,l),dd(1,1,5),ul(1,1,4),xl,
     &                    ndf,ndm,nel)
              sfac = d(78)*ctan(2)
            else
              sfac = 0.0d0
            endif

c           Compute gravity, thermal, inertia, and stress contributions

            call resid1d(cfac,lfac,xsj(l),xsj0,shp,sig(1,l),d,
     &                   ul(1,1,4),ul(1,1,5),p,ndf)

c           Tangent stiffness computation

            if(isw.eq.3) then

c             Modify tangent for stiffness rayleigh damping

              do j = 1,3
                do k = 1,3
                  dd(k,j,1) = dd(k,j,1)*ctan(1) + dd(k,j,5)*sfac
                end do ! k
c               Thermo-mechanical coupling
                dd(k,1,2) = dd(k,1,2)*ctan(1)
              end do ! j

              j1 = 1

              do j = 1,nel

                aj1 = shp(1,j,1)*xsj(l)
                aj2 = shp(2,j,1)*xsj0

c               Compute B_trans * D * j * w

                bd11 = aj1*dd(1,1,1) + aj2*dd(3,1,1)
                bd13 = aj1*dd(1,3,1) + aj2*dd(3,3,1)

c               Compute B_trans * D * alpha * j * w

                bt1  = aj1*dd(1,1,2) + aj2*dd(3,1,2)

c               Compute lumped mass matrix

                aj0          = shp(2,j,1)*dv
                s(j1  ,j1  ) = s(j1  ,j1  ) + aj0*lfac

c               Loop over columns (symmetry noted)

                k1 = 1
                do k = 1,nel
                  s(j1  ,k1  ) = s(j1  ,k1  ) + bd11*shp(1,k,1)
     &                                        + bd13*shpr(k)
     &                                        + cfac*aj0*shp(2,k,1)

c                 Add thermo-mechanical coupling term

                  if(ther) then
                    s(j1  ,k1+2) = s(j1  ,k1+2) + bt1*shp(2,k,1)
                  endif

                  k1 = k1 + ndf
                end do ! k

                j1 = j1 + ndf
              end do ! j
            end if
            if(cplxfl) then
              call cst1d1(shp,shpr,xsj(l),xsj0,xx,eps,ul(1,1,8),
     &                    dd,dd(1,1,3),s(1,nst+1),p,p(1,nen+1))
            endif
          end if
          nn = nn + nhv
        end do ! l

c     Output of element quantities

      elseif(isw.eq.4 .or. isw.eq.8 .or. isw.eq.16 .or. isw.eq.25) then

        if(nint(d(182)).eq.1) then
          lint = nel
          call int1dn(lint,sg)
        else
          lint = nint(d(5))
          call int1d(lint,sg)
        endif

c       Compute element stresses, strains, and forces

        nn = 0
        do l = 1,lint

c         Compute element shape functions

          call shp1d(sg(1,l),xl,shp(1,1,l),ndm,nel,xsj(l))
          xsj(l) = xsj(l)*sg(2,l)

c         Compute strains and coordinates

          call stra1d(d,xl,ul,tl,shp(1,1,l),ndf,ndm,nel,
     &                xx,ta,eps)
          call modlsd(d,ta,eps,hr(nh1+nn),hr(nh2+nn),nhv,istrt,
     &                dd,sig(1,l),isw)

c         Compute volumetric strain

          epsv(l) = eps(1,1) + eps(2,1) + eps(3,1)

          if(isw.eq.4) then

c           Output stresses and strains

            mct = mct - 2
            if(mct.le.0) then
            write(iow,2001) o,head
              if(ior.lt.0 .and. pfr) then
                write(*,2001) o,head
              endif
              mct = 50
            endif
            write(iow,2002)  n,ma,(sig(ii,l),ii=1,3),
     &                         xx,(eps(ii,1),ii=1,3)
            if(ior.lt.0 .and. pfr) then
              write(*,2002)  n,ma,(sig(ii,l),ii=1,3),
     &                         xx,(eps(ii,1),ii=1,3)
            endif
          endif
          nn = nn + nhv
        end do ! l

c       Compute nodal stress values

        if(isw.eq.8) then

          call slcn1d(sig,shp,xsj,p,s,p(nen+1,1),lint,nel,10)

c       Compute J-integrals and material forces

        elseif(isw.eq.16) then

          call pjint1d(d,ul,tl,shp,xsj,epsv,sig,p,ndf,ndm,lint)

c       Compute Z-Z projections

        elseif(isw.eq.25) then

          call stcn1z(xl,sig,shp,xsj,lint,ndm,nel,16)

        endif

c     Compute nodal stress error values

      elseif(isw.eq.11) then

        call ster1d(d,xl,ul,tl,shp,s,ndf,ndm,nel,nen)

      endif

c     Formats for input-output

2001  format(a1,20a4//5x,'Element Stresses'//'    Elmt Mat',
     &   '   11-stress   22-stress   33-stress'/
     &   '     1-coord   11-strain   22-strain   33-strain')
2002  format(i8,i4,1p,3e12.3/1p,4e12.3/1x)

      end
