c$Id: dplot.f,v 1.1 2006/11/20 20:33:21 rlt Exp $
      subroutine dplot(x,y,ipin)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose:  Line drawing routine for device types

c      Inputs:
c         x,y       - Screen coordinates for move/draw
c         ipin      - Switch: = 1 for panel initiation; = 2 for draw;
c                             = 3 for move (no draw).

c      Outputs:
c         none      - Outputs are graphics lines/moves
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'plclip.h'
      include  'pdata2.h'
      include  'pdatap.h'
      include  'pdataq.h'
      include  'pdatxt.h'
      include  'pdatps.h'
      include  'plflag.h'
      include  'x11f.h'

      integer   ipen,ipin
      real*4    x11, y11, xtem, ytem
      real*8    x, y

      save

c     Pen command motions  (ipen = 1, begin panel)
c                          (ipen = 2, drawline from to x,y)
c                          (ipen = 3, move to position x,y,
c                                     but do not draw line)
      ipen = abs(ipin)

      if(fwin .and. clchk) then
        call pclip(x,y,ipin)
      else

c       Set flag to indicate inside clip region (or not checking)

        clip = .true.

        xtem = 0.78125d0*x
        ytem = 0.78125d0*y

c       X11 device

        if (ipen .eq. 1) then
          xp(1) = xtem
          yp(1) = ytem
          ipan  = 1
        elseif(ipen .eq. 2) then
          if(ipan.eq.0) then
            xp(1) = xp(2)
            yp(1) = yp(2)
            xp(2) = xtem
            yp(2) = ytem
            x11 = xp(1)*xx(2)
            y11 = yp(1)*xx(3)*1.28
            if(screfl) then
              call gdx11(3,x11,y11)
            endif
            x11 = xp(2)*xx(2)
            y11 = yp(2)*xx(3)*1.28
            if(screfl) then
              call gdx11(4,x11,y11)
            endif
          elseif(ipan.ge.1) then
            ipan = ipan + 1
            xp(ipan) = xtem
            yp(ipan) = ytem
          endif
        elseif(ipen .eq. 3) then
          xp(1) = xp(2)
          yp(1) = yp(2)
          xp(2) = xtem
          yp(2) = ytem
        endif

c       Postcript

        if(hdcpy) then

c         Plot line segment

          if(ipan.eq.0 .and. ipen.eq.2) then
            call fpplps(2, xp, yp)
          endif

        endif

      endif

      end
