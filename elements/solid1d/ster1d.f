c$Id: ster1d.f,v 1.1 2006/11/20 20:32:56 rlt Exp $
      subroutine ster1d(d,xl,ul,th,shp,st,ndf,ndm,nel,nen)

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

      integer   ndf,ndm,nel,nen, i,j,ii,lint
      real*8    st(nen,*),xl(ndm,*),shp(2,4),sig(6), epsp(3)
      real*8    sigp(3),d(*),eps(9,3),ul(ndf,*),th(*)
      real*8    dd(6,6),dr(6,6),sg(2,4), xsj, xx,ta

      save

c     Error computation routine

      vfem   = 0.0d0
      vproj  = 0.0d0
      verror = 0.0d0
      vener  = 0.0d0
      venere = 0.0d0
      heta   = 0.0d0

c     Set quadrature formula

      if(nint(d(182)).eq.1) then
        lint = nel
        call int1dn(lint,sg)
      else
        lint = nint(d(5))
        call int1d(lint,sg)
      endif

      do ii = 1,lint
        call shp1d(sg(1,ii),xl,shp,ndm,nel,xsj)

c       Compute stresses

        call stra1d(d,xl,ul,th,shp,ndf,ndm,nel,xx,ta,eps)
        call estrsd(d,ta,eps,sig,dd,dr)

        do i = 1,3
          sigp(i)= 0.0d0
        end do ! i

        do i = 1,nel
          do j = 1,3
            sigp(j) = sigp(j) + shp(2,i)*st(i,j)
          end do ! j
        end do ! i

c       Compute projected strains

        call invert(dd,3,6)

        do i = 1,3
c         epsp(i) = ta
          epsp(i) = 0.0d0
        end do ! i

        do i = 1,3
          do j = 1,3
            epsp(i) = epsp(i) + dd(i,j)*sigp(j)
          end do ! j
        end do ! i

c       Compute integral of stress squares for error indicator use

        heta = heta + xsj

        do i = 1,3
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
      heta  =  d(50)*sqrt(heta)

      end
