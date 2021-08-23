c$Id: resid3d.f,v 1.1 2006/11/20 20:32:56 rlt Exp $
      subroutine resid3d(xsj,shp,sig,d,xl,vl,al,r,ndm,ndf)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: 3-D residual routine

c      Inputs:

c      Outputs:
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'eldata.h'
      include  'eltran.h'

      integer   ndm,ndf, i,j
      real*8    aj1,aj2,aj3,aj0,rr,lfac,cfac,xsj
      real*8    d(*),xl(ndm,*),vl(ndf,*),al(ndf,*),r(ndf,*)
      real*8    shp(4,*),sig(*),xx(3),ac(3),vc(3), bf(3),bt(3,3)

      save

c     Compute stress-divergence vector (p)

      do i = 1,ndm
        bf(i) = 0.0d0
      end do ! i
      call sbodyf(d, bf)

c     Angular velocity body force: d(4) = density, d(65) = omega

      rr = d(4)
      if(rr.gt.0.0d0 .and. d(65).gt.0.0d0) then
        do i = 1,3
          xx(i) = 0.0d0
          do j = 1,nel
            xx(i) = xx(i) + shp(4,j)*xl(i,j)
          end do ! j
        end do ! i
        call sbodyw(rr,d(65),xx, bf,bt, .false.)
      endif

      if(d(7).ge.0.0d0) then
        cfac = d(7)
        lfac = 1.d0 - cfac
      else
        cfac = 0.0d0
        lfac = 0.0d0
      endif

c     Compute accelerations

      do i = 1,3
        ac(i) = 0.0d0
        do j = 1,nel
          ac(i) = ac(i) + shp(4,j)*al(i,j)
        end do ! j
        ac(i)   = rr*ac(i)*cfac
      end do ! i

c     For Rayleigh Mass Damping: Compute velocity

      if(d(77).ne.0.0d0) then
        do i = 1,3
          vc(i) = 0.0d0
        end do ! i
        do j = 1,nel
          do i = 1,ndf
            vc(i) = vc(i) + shp(4,j)*vl(i,j)
          end do ! i
          aj0   = shp(4,j)*xsj*rr*d(77)
          do i = 1,3
            r(i,j) = r(i,j) - (cfac*vc(i) + lfac*vl(i,j))*aj0
          end do ! i
        end do ! j
      endif

c     Loop over rows

      do j = 1,nel
        aj1 = shp(1,j)*xsj
        aj2 = shp(2,j)*xsj
        aj3 = shp(3,j)*xsj
        aj0 = lfac*rr

c       Compute gravity, thermal, inertia, and stress contributions

        r(1,j) = r(1,j) + (bf(1) - ac(1) - aj0*al(1,j))*shp(4,j)*xsj
     &                  - aj1*sig(1)  - aj2*sig(4)  - aj3*sig(6)
        r(2,j) = r(2,j) + (bf(2) - ac(2) - aj0*al(2,j))*shp(4,j)*xsj
     &                  - aj1*sig(4)  - aj2*sig(2)  - aj3*sig(5)
        r(3,j) = r(3,j) + (bf(3) - ac(3) - aj0*al(3,j))*shp(4,j)*xsj
     &                  - aj1*sig(6)  - aj2*sig(5)  - aj3*sig(3)
      end do ! j

      end
