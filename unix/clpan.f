c$Id: clpan.f,v 1.1 2006/11/20 20:33:21 rlt Exp $
      subroutine clpan()

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Close panel plots

c      Inputs:
c         none

c      Outputs:
c         none
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'pdatap.h'
      include  'pdatps.h'
      include  'pdataq.h'
      include  'plflag.h'
      include  'psdat3.h'
      include  'x11f.h'

      integer   ii

      save

c     Close panel for filled plots

c     Fill panel using X11

      if(ipan.ge.3) then
        do ii=1, ipan
          x11(ii) = xp(ii)*xx(2)
          y11(ii) = yp(ii)*xx(3)*1.28
        end do ! ii

        if(screfl) call gdx11(ipan+1024,x11,y11)

c       Fill panel for PostScript

        if (hdcpy) call fppspl(ipan,xp,yp)

      endif

c     Reinitialize panel counter

      ipan = 0

      end
