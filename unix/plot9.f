c$Id: plot9.f,v 1.1 2006/11/20 20:33:21 rlt Exp $
      subroutine plot9(iel,ix,x,ndm,nel,isp)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Plot 2-d 3 to 9 node elements

c      Inputs:
c         iel       - Element type number
c         ix(*)     - Node list for element
c         x(ndm,*)  - Nodal coordinates for element
c         nel       - Number of nodes on element
c         isp       - Indicator on line type for outline of panels

c      Outputs:
c         none      - Plot outputs to screen/file
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'pdata2.h'
      include  'plclip.h'

      integer   iel, ndm, nel, isp, i, j, is, iu, ix(*),iplt(30)
      real*8    x(ndm,*),v

      save

      call pltord(ix,iel, iu,iplt)

c     Set for fill panel or line

      if(nel.gt.2) then
        is = abs(isp)
      else
        is      = 3
        iplt(1) = 1
        iplt(2) = 2
        iu      = 2
      endif

      v = 0.0d0
      if(ndm.eq.3) v = x(3,iplt(1))
      call plotl(x(1,iplt(1)),x(2,iplt(1)),v,is)
      do i = 2,iu
        j = iplt(i)
        if((j.le.nel).and.(j.ne.0).and.(ix(j).ne.0)) then
          if(ndm.eq.3) v = x(3,j)
          call plotl(x(1,j),x(2,j),v,2)
        endif
      end do ! i

c     Fill in corners for a clip and close panel

      if(is.eq.1) then

        if(fwin) then
          if(ndm.eq.3) v = x(3,iplt(1))
          call plotl(x(1,iplt(1)),x(2,iplt(1)),v,2)
          if(ienter.gt.0) then
            call pfclip(iexit,ienter)
          endif
        endif

        call clpan

c       Outline part in black or white

        if(isp.gt.0) then
          is = icolr

          if(icolr.ne.0) then
            call pppcol(0,1)
          else
            call pppcol(1,1)
          endif

          v = 0.0d0
          if(ndm.eq.3) v = x(3,iplt(1))
          call plotl(x(1,iplt(1)),x(2,iplt(1)),v,3)
          do i = 2,iu
            j = iplt(i)
            if((j.le.nel).and.(j.ne.0).and.(ix(j).ne.0)) then
              if(ndm.eq.3) v = x(3,j)
              call plotl(x(1,j),x(2,j),v,2)
            endif
          end do ! i
          icolr = is
        endif

      endif

      end
