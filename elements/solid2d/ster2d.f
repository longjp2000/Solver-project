c$Id: ster2d.f,v 1.1 2006/11/20 20:32:56 rlt Exp $
      subroutine ster2d(ix,d,xl,ul,th,shp,st,ndf,ndm,nel,nen)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose:

c      Inputs:

c      Outputs:
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'adapt1.h'
      include  'adapt2.h'
      include  'errind.h'
      include  'hdata.h'
      include  'comblk.h'

      logical   quad
      integer   ndf,ndm,nel,nen, i,j,ii,lint, ix(*)
      real*8    st(nen,*),xl(ndm,*),shp(3,4),sig(6), epsp(4)
      real*8    sigp(4),d(*),eps(9,3),ul(ndf,*),th(*)
      real*8    dd(6,6),dr(6,6),sg(3,16),el(4,7), xsj, xx,yy, ta

      save

c     Error computation routine

      vfem   = 0.0d0
      vproj  = 0.0d0
      verror = 0.0d0
      vener  = 0.0d0
      venere = 0.0d0
      heta   = 0.0d0

c     Set quadrature formula

      if(nel.eq.3 .or. nel.eq.6 .or. nel.eq.7) then
        quad = .false.
        call tint2dn(nel,lint,el)
      else
        quad = .true.
        if(nint(d(182)).gt.0) then
          call int2dn(nel,lint,sg)
        else
          i    = nint(d(5))
          call int2d(i,lint,sg)
        endif
      endif
      do ii = 1,lint
        if(quad) then
          call shp2d(sg(1,ii),xl,shp,xsj,ndm,nel,ix,.false.)
        else
          call trishp(el(1,ii),xl,ndm,nel-4,xsj,shp)
        endif

c       Compute stresses

        call strn2d(d,xl,ul,th,shp,ndf,ndm,nel,xx,yy,ta,eps)
        call estrsd(d,ta,eps,sig,dd,dr)

        do i = 1,4
          sigp(i)= 0.0d0
        end do ! i

        do i = 1,nel
          do j = 1,4
            sigp(j) = sigp(j) + shp(3,i)*st(i,j)
          end do ! j
        end do ! i

c       Compute projected strains

        call invert(dd,4,6)

        do i = 1,3
c         epsp(i) = ta
          epsp(i) = 0.0d0
        end do ! i
        epsp(4) = 0.0d0

        do i = 1,4
          do j = 1,4
            epsp(i) = epsp(i) + dd(i,j)*sigp(j)
          end do ! j
        end do ! i

c       Compute integral of stress squares for error indicator use

        heta = heta + xsj

        do i = 1,4
          vfem   = vfem   + sig(i)*sig(i)*xsj
          vproj  = vproj  + sigp(i)*sigp(i)*xsj
          verror = verror + ((sigp(i)-sig(i))**2)*xsj
          vener  = vener  + sig(i)*eps(i,1)*xsj
          venere = venere + (sigp(i)-sig(i))*(epsp(i)-eps(i,1))*xsj
        end do ! i
      end do ! ii

c     Set error indicators

      arsq   = arsq  + heta
      efem   = efem  + vfem
      eproj  = eproj + vproj
      eerror = eerror+ verror
      eener  = eener + vener
      eenere = eenere+ venere

      areai  = heta

c     Check for triangles

      if(nel.eq.3 .or. nel.eq.6 .or. nel.eq.7) then
        heta = heta + heta
      endif

      heta  =  d(50)*sqrt(heta)

      end
