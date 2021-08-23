c$Id: plstrt.f,v 1.1 2006/11/20 20:33:21 rlt Exp $
      subroutine plstrt()

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Open a graphics device to receive plot data

c      Inputs:
c         none

c      Outputs:
c         none      - Graphics window should appear after call
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'debugs.h'
      include  'pdata2.h'
      include  'plflag.h'
      include  'x11f.h'

      integer   i

      save

c     Open devices for Graphics output

      everon = .true.

c     X11 device


      idx = 440
      idy = 550

      if(screfl) call gdx11(1,xx,yy)
      if(screfl) call gdx11(7,xx,yy)
      if(debug) then
        write(*,2000) (xx(i),i = 2,6)
      endif

c     Format

2000  format('  X Length in cm.     ',f10.5/
     &       '  Y Length in cm.     ',f10.5/
     &       '  X Pixels in cm.     ',f10.5/
     &       '  Y Pixels in cm.     ',f10.5/
     &       '  No. Forground colors',f10.5/)

      end
