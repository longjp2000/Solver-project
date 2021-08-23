c$Id: lclip.f,v 1.1 2006/11/20 20:33:12 rlt Exp $
      logical function lclip( ix, nen, x, ndm )

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Clip plot at specified planes

c      Inputs:
c         ix(*)     - List of nodes on element
c         nen       - Number of nodes on element
c         x(ndm,*)  - Nodal coordinates of element
c         ndm       - Spatial dimension of mesh

c      Outputs:
c         lclip     - Flag, true if element is within clip region
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'plclip.h'

      integer   i, n, nen,ndm
      integer   ix(nen)
      real*8    x(ndm,*),x0(3)

      save

      do i = 1, ndm
        x0(i) = 0.0d0
        do n = 1,nen
          if(ix(n).gt.0) x0(i) = x0(i) + x(i,ix(n))
        end do ! n
        x0(i) = x0(i)/nen
      end do ! i
      if(ndm.eq.1) then
        lclip = (x0(1).ge.cmin(1)) .and. (x0(1).le.cmax(1))
      elseif(ndm.eq.2) then
        lclip = (x0(1).ge.cmin(1)) .and. (x0(1).le.cmax(1))
     &                             .and.
     &          (x0(2).ge.cmin(2)) .and. (x0(2).le.cmax(2))
      elseif(ndm.eq.3) then
        lclip = (x0(1).ge.cmin(1)) .and. (x0(1).le.cmax(1))
     &                             .and.
     &          (x0(2).ge.cmin(2)) .and. (x0(2).le.cmax(2))
     &                             .and.
     &          (x0(3).ge.cmin(3)) .and. (x0(3).le.cmax(3))
      else
        lclip = .false.
      endif

      end
