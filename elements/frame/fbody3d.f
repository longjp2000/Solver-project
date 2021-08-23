c$Id: fbody3d.f,v 1.1 2006/11/20 20:32:56 rlt Exp $
      subroutine fbody3d(d,xl, r, ndm,ndf, isw)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Compute body force loads
c-----[--.----+----.----+----.-----------------------------------------]

      implicit   none

      include   'elbody.h'
      include   'eldata.h'

      integer    ndm,ndf,isw, ii
      real*8     d(*),xl(ndm,*), r(ndf,*), body(3), le

c     Set body loading factors

      if(isw.eq.15) then
        do ii = 1,3
          body(ii) = bodyf(ii)
        end do ! ii
      else
        call sbodyf(d, body)
      endif

c     Add body forces

      le       = sqrt((xl(1,2)-xl(1,1))**2 + (xl(2,2)-xl(2,1))**2
     &              + (xl(3,2)-xl(3,1))**2)*0.5d0
      do ii = 1,3
        r(ii,1)   = r(ii,1) + body(ii)*le
        r(ii,2)   = r(ii,2) + body(ii)*le
      end do ! ii

      end
